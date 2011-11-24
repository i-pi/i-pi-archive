import numpy, math
import engine, io_system, cell, upper_T, dynamics
from utils.depend import *
from utils import units
      
class rp_nve_ensemble(dynamics.nve_ensemble):
   def __init__(self, syst, dt=1.0):
      super(rp_nve_ensemble,self).__init__(syst=syst, dt=dt)

   def step(self):

      dt = self.dt.get()

      self.syst.p.get_array()[:] += self.syst.f.get_array()*dt/2.0

      q_tilde = numpy.dot(self.syst.trans_mat.get_array(), self.syst.q.get_array())
      p_tilde = numpy.dot(self.syst.trans_mat.get_array(), self.syst.p.get_array())

      nbeads = len(self.syst.systems)
      natoms = len(self.syst.atoms)
      for i in range(natoms):
         cos_mat = numpy.zeros(nbeads)
         sin_mat = numpy.zeros(nbeads)
         neg_sin_mat = numpy.zeros(nbeads)
         mass = self.syst.atoms[i].mass.get()

         for j in range(nbeads):
            omega_j = self.syst.n_frequencies.get_array()[j]
            cos_mat[j] = math.cos(dt*omega_j)
            sin_mat[j] = math.sin(dt*omega_j)/(mass*omega_j)
            neg_sin_mat[j] = -mass*omega_j*math.sin(dt*omega_j)
         sin_mat[0] = dt/mass

         p_tilde_i = numpy.array(p_tilde[:,3*i:3*(i+1)])
         q_tilde_i = numpy.array(q_tilde[:,3*i:3*(i+1)])


         for j in range(nbeads):
            p_tilde[j,3*i:3*(i+1)] = cos_mat[j]*p_tilde_i[j,:] + neg_sin_mat[j]*q_tilde_i[j,:]
            q_tilde[j,3*i:3*(i+1)] = cos_mat[j]*q_tilde_i[j,:] + sin_mat[j]*p_tilde_i[j,:]

      self.syst.p.get_array()[:] = numpy.dot(numpy.transpose(self.syst.trans_mat.get_array()), p_tilde)
      self.syst.q.get_array()[:] = numpy.dot(numpy.transpose(self.syst.trans_mat.get_array()), q_tilde)



      self.syst.q.taint(taintme=False)
      self.syst.p.get_array()[:] += self.syst.f.get_array()*dt/2.0
      self.syst.p.taint(taintme=False)

class rp_nvt_ensemble(dynamics.nvt_ensemble,rp_nve_ensemble):
   """NVT ensemble object, with velocity Verlet time integrator and thermostat
      Contains: syst = System object, containing the atom and cell coordinates 
      econs = conserved energy quantity, dt = time step, temp = temperature,
      thermo = thermostat object for the atoms.
      Initialised by: ens = ensemble(system, thermo, temp, dt)
      system is a System object, containing the atom and cell coordinates
      thermo is a thermostat object, which maintains constant temperature
      temp = temperature, default = 1.0
      dt = time step, default = 1.0"""

   def __init__(self, syst, thermo, temp=1.0, dt=1.0):
      super(rp_nvt_ensemble,self).__init__(syst, thermo, temp=temp, dt=dt)

      self.thermo.temp.set(self.temp*len(syst.systems))

   def step(self):
      self.thermo.step()

      rp_nve_ensemble.step(self)
      #dynamics.nve_ensemble.step(self)

      self.thermo.step()

#class npt_ensemble(nvt_ensemble):
##TODO rework this entirely, so that we separate the NPT and NST implementations
#   """NPT ensemble object, with Bussi time integrator and cell dynamics, 
#      with independent thermostating for the cell and atoms.
#      Contains: syst = System object, containing the atom and cell coordinates 
#      econs = conserved energy quantity, dt = time step, temp = temperature,
#      thermo = thermostat object for the atoms, cell_thermo = thermostat
#      object for the cell.
#      Initialised by: 
#      ens = ensemble(system, thermo, cell_thermo, temp, dt, pext)
#      system is a System object, containing the atom and cell coordinates
#      thermo is a thermostat object, which maintains constant temperature for 
#      the atoms
#      cell_thermo maintains constant temperature for the cell
#      temp = temperature, default = 1.0
#      dt = time step, default = 1.0
#      pext = external pressure, default = 0"""
#
#   def __init__(self, syst, thermo, cell_thermo, temp=1.0, dt=1.0, pext=0.0):
#      super(npt_ensemble,self).__init__(syst=syst, thermo=thermo, temp=temp, dt=dt)
#      self.syst.init_cell_velocities(temp = temp)
#
#      self.cell_thermo=cell_thermo
#      self.cell_thermo.bind(self.syst.atoms, self.syst.p, self.syst.cell)
#      self.cell_thermo.dt.set(self.dt.get()*0.5)   # maybe make thermo.dt a dependant of dt?
#      self.cell_thermo.temp.set(self.temp)
#      self.syst.cell.pext.set(pext*numpy.identity(3,float)) 
#      self.econs=depend(name='econs',func=self.get_econs, deplist=[ self.syst.pot, self.syst.kin, self.thermo.econs, self.cell_thermo.econs, self.syst.cell.kin, self.syst.cell.pot])
#
#   def pstep(self):
#      """Evolves the atom and cell momenta forward in time by a step dt/2"""
#
#      p = self.syst.p.get(); f = self.syst.f.get(); pc=(self.syst.cell.pc.get())
#      V=self.syst.cell.V.get(); dthalf=self.dt.get()*0.5
#      press=self.syst.press.get();  pext=numpy.trace(self.syst.cell.pext.get())/3.0
#
#      pc += dthalf * 3.0* (V* ( press-pext ) + 2.0*units.kb*self.temp)
#            
#      for i in range(len(self.syst.atoms)):
#         atom_i = self.syst.atoms[i]
#         pc += dthalf**2/atom_i.mass.get()*numpy.inner(atom_i.f.get(), atom_i.p.get())
#         pc += dthalf**3/(3.0*atom_i.mass.get())*numpy.inner(atom_i.f.get(), atom_i.f.get())
#
#      self.syst.cell.pc.set(pc)
#      p[:] += f[:] * dthalf;   
#      self.syst.p.taint(taintme=False)
#
#   def rstep(self):
#      """Takes the atom positions, velocities and forces and integrates the 
#         equations of motion forward by a step dt"""
#
#      expeta = math.exp(self.syst.cell.pc.get()*self.dt.get()/self.syst.cell.w.get())
#      sinheta=0.5*(expeta-1.0/expeta)
#      for i in range(len(self.syst.atoms)):
#         atom_i = self.syst.atoms[i]
#         (atom_i.q.get())[:] *= expeta 
#         (atom_i.q.get())[:] += (atom_i.p.get())[:]/atom_i.mass.get() *sinheta/(self.syst.cell.pc.get()/self.syst.cell.w.get())
#         (atom_i.p.get())[:] /= expeta
#
#      self.syst.q.taint(taintme=False);  self.syst.p.taint(taintme=False)    
#      (self.syst.cell.h.get())[:] *= expeta
#      self.syst.cell.h.taint(taintme=False)
#      
#      
#   def step(self):
#      """NPT time step, with appropriate thermostatting steps"""
#
#      self.cell_thermo.NPT_cell_step()
#      self.thermo.step()
#      self.pstep()
#      self.rstep()
#      self.pstep()
#      self.thermo.step()
#      self.cell_thermo.NPT_cell_step()
#   
#   def get_econs(self): 
#      """Calculates the conserved energy quantity for the NPT ensemble"""
#
#      return nvt_ensemble.get_econs(self)-2.0*units.kb*self.temp*math.log(self.syst.cell.V.get())+self.syst.cell.pext.get()[0,0]*self.syst.cell.V.get()+self.cell_thermo.econs.get()+self.syst.cell.kin.get()
#
#class nst_ensemble(nvt_ensemble):
#   """NST ensemble object, with Bussi time integrator and cell dynamics, 
#      with independent thermostating for the cell and atoms.
#      Contains: syst = System object, containing the atom and cell coordinates 
#      econs = conserved energy quantity, dt = time step, temp = temperature,
#      thermo = thermostat object for the atoms, cell_thermo = thermostat
#      object for the cell.
#      Initialised by: 
#      ens = ensemble(system, thermo, cell_thermo, temp, dt, pext)
#      system is a System object, containing the atom and cell coordinates
#      thermo is a thermostat object, which maintains constant temperature for 
#      the atoms
#      cell_thermo maintains constant temperature for the cell
#      temp = temperature, default = 1.0
#      dt = time step, default = 1.0
#      pext = external pressure tensor"""
#
#   def __init__(self, syst, thermo, cell_thermo, temp=1.0, dt=1.0, pext=numpy.zeros((3,3),float)):
#      super(nst_ensemble,self).__init__(syst=syst, thermo=thermo, temp=temp, dt=dt)
#      self.syst.init_cell_velocities(temp = temp)
#
#      self.cell_thermo=cell_thermo
#      self.cell_thermo.bind(self.syst.atoms, self.syst.p, self.syst.cell)
#      self.cell_thermo.dt.set(self.dt.get()*0.5)   # maybe make thermo.dt a dependant of dt?
#      self.cell_thermo.temp.set(self.temp)
#      self.syst.cell.pext.set(pext) 
#      self.econs=depend(name='econs',func=self.get_econs, deplist=[ self.syst.pot, self.syst.kin, self.thermo.econs, self.cell_thermo.econs, self.syst.cell.kin, self.syst.cell.pot])
#      
##   def exp_p(self):
##      """Exponentiates the displacement matrix p*dt/w, as required for rstep"""
##
##      p=self.syst.cell.p.get()
##      dist_mat = p*self.dt.get()/self.syst.cell.w.get()
##      eig, eigvals = upper_T.compute_eigp(dist_mat)
##      i_eig = upper_T.compute_ih(eig)
##   
##      exp_mat = numpy.zeros((3,3), float)
##      neg_exp_mat = numpy.zeros((3,3), float)
##      for i in range(3):
##         exp_mat[i,i] = math.exp(eigvals[i])
##         neg_exp_mat[i,i] = math.exp(-eigvals[i])
##      
##      exp_mat = numpy.dot(eig, exp_mat)
##      exp_mat = numpy.dot(exp_mat, i_eig)
##      
##      neg_exp_mat = numpy.dot(eig, neg_exp_mat)
##      neg_exp_mat = numpy.dot(neg_exp_mat, i_eig)
##
##      return exp_mat, neg_exp_mat
#
#   def pstep(self):
#      """Evolves the atom and cell momenta forward in time by a step dt/2"""
#
#      p = self.syst.p.get(); f = self.syst.f.get(); pc=self.syst.cell.p.get()
#      V=self.syst.cell.V.get(); dthalf=self.dt.get()*0.5
#      L = numpy.zeros((3,3), float)
#      for i in range(3):
#         L[i,i] = 3.0-i
#
#      pc[:] = pc[:] + dthalf*(V*( self.syst.stress.get() - self.syst.cell.piext.get()) + 2.0*units.kb*self.temp*L)
#            
#      for i in range(len(self.syst.atoms)):
#         atom_i = self.syst.atoms[i]
#         pc[:] += dthalf**2/(2.0*atom_i.mass.get())*(numpy.outer(atom_i.f.get(), atom_i.p.get()) + numpy.outer(atom_i.p.get(), atom_i.f.get()))
#         pc[:] += dthalf**3/(3.0*atom_i.mass.get())*numpy.outer(atom_i.f.get(), atom_i.f.get())
#      self.syst.cell.p.taint(taintme=False)      
#
#      p[:] += f[:] * dthalf
#      self.syst.p.taint(taintme=False) 
#
#   def rstep(self):
#      """Takes the atom positions, velocities and forces and integrates the 
#         equations of motion forward by a step dt"""
#
#      vel_mat = self.syst.cell.p.get()/self.syst.cell.w.get()
#      exp_mat, neg_exp_mat = upper_T.Crank_Nicolson(vel_mat*self.dt.get())
#      sinh_mat = 0.5*(exp_mat - neg_exp_mat)
#      ip_mat = cell.ut_inverse(vel_mat)
#
#      for i in range(len(self.syst.atoms)):
#         atom_i = self.syst.atoms[i]
#         (atom_i.q.get())[:] = numpy.dot(exp_mat, atom_i.q.get()) + numpy.dot(ip_mat, numpy.dot(sinh_mat, atom_i.p.get()/atom_i.mass.get()))
#         (atom_i.p.get())[:] = numpy.dot(neg_exp_mat, atom_i.p.get())
#
#      self.syst.q.taint(taintme=False);   self.syst.p.taint(taintme=False)          
#      (self.syst.cell.h.get())[:] = numpy.dot(exp_mat, self.syst.cell.h.get())
#      self.syst.cell.h.taint(taintme=False)
#      
#      
#   def step(self):
#      """NST time step, with appropriate thermostatting steps"""
#
#      self.cell_thermo.NST_cell_step()
#      self.thermo.step()
#      self.pstep()
#      self.rstep()
#      self.pstep()
#      self.thermo.step()
#      self.cell_thermo.NST_cell_step()
#   
#   def get_econs(self):
#      """Calculates the conserved energy quantity for the NST ensemble"""
#
#      xv=0.0 # extra term stemming from the Jacobian
#      for i in range(3): xv+=math.log(self.syst.cell.h.get()[i,i])*(3-i) 
#      
#      return nvt_ensemble.get_econs(self)+self.cell_thermo.econs.get()+self.syst.cell.pot.get()+self.syst.cell.kin.get()-2.0*units.kb*self.temp*xv
#     
