import numpy, math
import engine, io_system, cell, upper_T
from utils.depend import *
from utils import units

class ensemble(object): 
   def __init__(self,syst):
      self.syst=syst
      self.econs=depend(name='econs',func=self.get_econs, deplist=[ self.syst.pot, self.syst.kin])
      
   def step(self): pass

   def get_econs(self):
      return self.syst.pot.get()+self.syst.kin.get()

class nve_ensemble(ensemble):
   def __init__(self, syst, dt=1.0):
      super(nve_ensemble,self).__init__(syst = syst)
      self.dt = depend(name='dt',value=dt)
   
   def step(self):
      p=self.syst.p.get();  q=self.syst.q.get();  dt=self.dt.get()
      
      p[:] += self.syst.f.get() * (dt*0.5);  # self.syst.p.taint(taintme=False) this taint step is likely to be unnecessary
      q[:] += p[:] * dt;                      self.syst.q.taint(taintme=False)            
      p[:] += self.syst.f.get() * (dt*0.5);   self.syst.p.taint(taintme=False) 
      
      self.syst.p.taint(taintme=False)
      
class nvt_ensemble(nve_ensemble):
   def __init__(self, syst, thermo, temp=1.0, dt=1.0):
      super(nvt_ensemble,self).__init__(syst=syst, dt=dt)
      self.temp=temp
      
      #hooks to system, thermostat, etc
      self.thermo = thermo
      self.thermo.bind(self.syst.atoms, self.syst.p, self.syst.cell)
      self.thermo.dt.set(self.dt.get()*0.5)   # maybe make thermo.dt a dependant of dt?
      self.thermo.temp.set(self.temp)
      self.econs=depend(name='econs',func=self.get_econs, deplist=[ self.syst.pot, self.syst.kin, self.thermo.econs])
      
   def step(self):
      self.thermo.step()
      nve_ensemble.step(self)
      self.thermo.step()
      
   def get_econs(self):
      return nve_ensemble.get_econs(self)+self.thermo.econs.get()


#TODO rework to be a derived class of nvt_ensemble, trying to keep the integrator as transparent as possible
class nst_ensemble(nvt_ensemble):
   def __init__(self, syst, thermo, cell_thermo, temp=1.0, dt=1.0, pext=numpy.zeros((3,3),float)):
      super(nst_ensemble,self).__init__(syst=syst, thermo=thermo, temp=temp, dt=dt)
      self.cell_thermo=cell_thermo
      self.cell_thermo.bind(self.syst.atoms, self.syst.p, self.syst.cell)
      self.cell_thermo.dt.set(self.dt.get()*0.5)   # maybe make thermo.dt a dependant of dt?
      self.cell_thermo.temp.set(self.temp)
      self.syst.cell.pext.set(pext) 
      self.econs=depend(name='econs',func=self.get_econs, deplist=[ self.syst.pot, self.syst.kin, self.thermo.econs, self.cell_thermo.econs, self.syst.cell.kin, self.syst.cell.pot])
      
   def exp_p(self):
      p=self.syst.cell.p.get()
      dist_mat = p*self.dt.get()/self.syst.cell.w.get()
      eig, eigvals = upper_T.compute_eigp(dist_mat)
      i_eig = upper_T.compute_ih(eig)
   
      exp_mat = numpy.zeros((3,3), float)
      neg_exp_mat = numpy.zeros((3,3), float)
      for i in range(3):
         exp_mat[i,i] = math.exp(eigvals[i])
         neg_exp_mat[i,i] = math.exp(-eigvals[i])
      
      exp_mat = numpy.dot(eig, exp_mat)
      exp_mat = numpy.dot(exp_mat, i_eig)
      
      neg_exp_mat = numpy.dot(eig, neg_exp_mat)
      neg_exp_mat = numpy.dot(neg_exp_mat, i_eig)

      return exp_mat, neg_exp_mat

   def pstep(self):
      #equivalent to P-step in paper

      p = self.syst.p.get(); f = self.syst.f.get(); pc=self.syst.cell.p.get()
      V=self.syst.cell.V.get(); dthalf=self.dt.get()*0.5
      L = numpy.zeros((3,3), float)
      for i in range(3):
         L[i,i] = 3.0-i

      pc[:] = pc[:] + dthalf*(V*( self.syst.stress.get() - self.syst.cell.piext.get()) + 2.0*units.kb*self.temp*L)
            
      for i in range(len(self.syst.atoms)):
         atom_i = self.syst.atoms[i]
         pc[:] += dthalf**2/(2.0*atom_i.mass.get())*(numpy.outer(atom_i.f.get(), atom_i.p.get()) + numpy.outer(atom_i.p.get(), atom_i.f.get()))
         pc[:] += dthalf**3/(3.0*atom_i.mass.get())*numpy.outer(atom_i.f.get(), atom_i.f.get())

      self.syst.cell.p.taint(taintme=False)      
      p[:] += f[:] * dthalf;    self.syst.p.taint(taintme=False) 

   def rstep(self):
      """Takes the atom positions, velocities and forces and integrates the 
         equations of motion forward by a step dt"""

      #equivalent to R-step in paper

      exp_mat, neg_exp_mat = self.exp_p()
      sinh_mat = 0.5*(exp_mat - neg_exp_mat)
      ip_mat = cell.ut_inverse(self.syst.cell.p.get()/self.syst.cell.w.get())

      for i in range(len(self.syst.atoms)):
         atom_i = self.syst.atoms[i]
         (atom_i.q.get())[:] = numpy.dot(exp_mat, atom_i.q.get()) + numpy.dot(ip_mat, numpy.dot(sinh_mat, atom_i.p.get()/atom_i.mass.get()))
         (atom_i.p.get())[:] = numpy.dot(neg_exp_mat, atom_i.p.get())

      self.syst.q.taint(taintme=False)         
      (self.syst.cell.h.get())[:] = numpy.dot(exp_mat, self.syst.cell.h.get())
      self.syst.cell.h.taint(taintme=False)
      
      
   def step(self):
      self.cell_thermo.cell_step()
      self.thermo.step()
      self.pstep()
      
      self.rstep()
      
      self.pstep()
      self.thermo.step()
      self.cell_thermo.cell_step()
   
   def get_econs(self):
      return nvt_ensemble.get_econs(self)+self.cell_thermo.econs.get()+self.syst.cell.pot.get()+self.syst.cell.kin.get()
      
        
   
   def thermo_step(self):
      self.thermo.step()
      self.thermo.cell_step(self.syst.cell)


#   def apply_pbc(self):
#      """Takes the system and applies periodic boundary conditions to fold the
#         particle positions back into the unit cell"""

#      for i in range(self.syst.natoms):
#         q = self.syst.cell.apply_pbc(self.syst.atoms[i])
#         for j in range(3):
#            self.syst.atoms[i].q.set(q)

   def R_update(self):
      self.pot_func.force_update()
      self.syst.get_kinetic()
      self.syst.stress_update()
      self.syst.cell_update()
      self.syst.tot_E_update()

   def TP_update(self):
      self.syst.get_kinetic()
      self.syst.stress_update()
      self.syst.cell_update()
      self.syst.tot_E_update()

   def simulation(self, maxcount = 5):
      self.R_update()
      print self.syst
#      f = open("./pdboutput.pdb","a")
      for i in range(maxcount):
         self.thermo_step()
         self.TP_update()
         self.vel_step()
         self.TP_update()
         self.pos_step()
         self.R_update()

#         if (i%10 == 0):
#            io_system.print_pdb(self.syst.atoms, self.syst.cell, f)

         self.vel_step()
         self.TP_update()
         self.thermo_step()
         self.TP_update()
         print self.syst
      self.apply_pbc()
      #print self.syst
      io_system.print_pdb(self.syst.atoms, self.syst.cell)
