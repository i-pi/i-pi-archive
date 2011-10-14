import numpy, math
import engine, io_system, cell

class NST_ens(object):

   @classmethod
   def from_pdbfile(cls, filedesc, thermo, pot_func, temp = 1.0, dt = 0.1, **kwargs):
      cls.dt = dt
      cls.temp = temp
      cls.syst = engine.System.from_pdbfile(filedesc, temp)
      cls.pot_func = pot_func(cls.syst, **kwargs)
      cls.thermo = thermo(temp, dt/2.0)
      return cls()

   @classmethod
   def from_ensemble(cls, ens):
      cls.dt = ens.dt
      cls.temp = ens.temp
      cls.syst = engine.System.from_system(ens.syst)
      cls.pot_func = ens.pot_func
      cls.thermo = ens.thermo
      return cls()

   def exp_p(self):
      dist_mat = self.syst.cell.p*self.dt/self.syst.cell.w
      eig = cell.compute_eigp(dist_mat)
      i_eig = cell.compute_ih(eig)
   
      exp_mat = numpy.zeros((3,3), float)
      neg_exp_mat = numpy.zeros((3,3), float)
      for i in range(3):
         exp_mat[i,i] = math.exp(self.syst.cell.p[i,i]*self.dt/self.syst.cell.w)
         neg_exp_mat[i,i] = math.exp(-self.syst.cell.p[i,i]*self.dt/self.syst.cell.w)
      
      exp_mat = numpy.dot(eig, exp_mat)
      exp_mat = numpy.dot(exp_mat, i_eig)
      
      neg_exp_mat = numpy.dot(eig, neg_exp_mat)
      neg_exp_mat = numpy.dot(neg_exp_mat, i_eig)

      return exp_mat, neg_exp_mat

   def thermo_step(self):
      for i in range(self.syst.natoms):
         self.thermo.step(self.syst.atoms[i])
      self.thermo.cell_step(self.syst.cell)

   def pos_step(self):
      """Takes the atom positions, velocities and forces and integrates the 
         equations of motion forward by a step dt"""

      #equivalent to R-step in paper

      exp_mat, neg_exp_mat = self.exp_p()
      sinh_mat = 0.5*(exp_mat - neg_exp_mat)
      ip_mat = cell.compute_ih(self.syst.cell.p/self.syst.cell.w)

      for i in range(self.syst.natoms):
         atom_i = self.syst.atoms[i]
         q = numpy.dot(exp_mat, atom_i.q) + numpy.dot(ip_mat, numpy.dot(sinh_mat, atom_i.p/atom_i.mass))
         p = numpy.dot(neg_exp_mat, atom_i.p)
         for j in range(3):
            atom_i.q[j] = q[j]
            atom_i.p[j] = p[j]

      self.syst.cell.h = numpy.dot(exp_mat, self.syst.cell.h)

   def vel_step(self):
      #equivalent to P-step in paper

      for i in range(self.syst.natoms*3):
         self.syst.p[i] = self.syst.p[i] + self.syst.f[i]*self.dt/2.0

      L = numpy.zeros((3,3), float)
      for i in range(3):
         L[i,i] = 3.0-i

      p = self.syst.cell.p + self.dt/2.0*(self.syst.cell.V*(self.syst.stress-self.syst.cell.PI_ext) + 2.0*self.thermo.k_Boltz*self.temp*L)
      for i in range(self.syst.natoms):
         atom_i = self.syst.atoms[i]
         p += (self.dt/2.0)**2/(2.0*atom_i.mass)*(numpy.outer(atom_i.f, atom_i.p) + numpy.outer(atom_i.p, atom_i.f))
         p += (self.dt/2.0)**3/(3.0*atom_i.mass)*numpy.outer(atom_i.f, atom_i.f)

      p[1,0]=p[2,0]=p[2,1]=0.0

      self.syst.cell.p = p

   def apply_pbc(self):
      """Takes the system and applies periodic boundary conditions to fold the
         particle positions back into the unit cell"""

      for i in range(self.syst.natoms):
         q = self.syst.cell.apply_pbc(self.syst.atoms[i])
         for j in range(3):
            self.syst.atoms[i].q[j] = q[j]

   def R_update(self):
      print "pre-force update"
      print self.syst.f
      self.pot_func.force_update()
      print "post-force update"
      print self.syst.f
      self.syst.kinetic_update()
      self.syst.stress_update()
      self.syst.cell_update()
      self.syst.tot_E_update()

   def TP_update(self):
      self.syst.kinetic_update()
      self.syst.stress_update()
      self.syst.cell_update()
      self.syst.tot_E_update()

   def simulation(self, maxcount = 5):
      self.R_update()
      print
      print self.syst
      print self.syst.f
      print self.syst.q
      print self.syst.p
      print
      for i in range(maxcount):
         self.thermo_step()
         self.TP_update()
         self.vel_step()
         self.TP_update()
         self.pos_step()
         print "pre r step:"
         print self.syst
         self.R_update()
         print "post r step:"
         print self.syst
         self.vel_step()
         self.TP_update()
         self.thermo_step()
         self.TP_update()
      #   print self.syst
      self.apply_pbc()
      print
      print self.syst
      print self.syst.f
      print self.syst.q
      print self.syst.p
      print
      io_system.print_pdb(self.syst.atoms, self.syst.cell)
