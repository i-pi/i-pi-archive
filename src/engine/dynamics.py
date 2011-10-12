import numpy, math
import engine, io_system, cell

class NST_ens(object):

   @classmethod
   def from_pdbfile(cls, filedesc, thermo, temp = 1.0, dt = 0.1):
      cls.dt = dt
      cls.syst = engine.System.from_pdbfile(filedesc, temp)
      cls.thermo = thermo(temp, dt/2.0)
      return cls()

   @classmethod
   def from_system(cls, system, thermo, dt = 0.1):
      cls.dt = dt
      cls.syst = engine.System.from_system(system)
      cls.thermo = thermo(system.temp, dt/2.0)
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
     # self.thermo.step(syst.cell)

   def pos_step(self):
      """Takes the atom positions, velocities and forces and integrates the 
         equations of motion forward by a step dt"""

      #equivalent to R-step in paper

      exp_mat, neg_exp_mat = self.exp_p()
      sinh_mat = 0.5*(exp_mat - neg_exp_mat)
      ip_mat = cell.compute_ih(self.syst.cell.p/self.syst.cell.w)

      for i in range(self.syst.natoms):
         self.syst.atoms[i].q = numpy.dot(exp_mat, self.syst.atoms[i].q) + numpy.dot(ip_mat, numpy.dot(sinh_mat, self.syst.atoms[i].p/self.syst.atoms[i].mass))
         self.syst.atoms[i].p = numpy.dot(neg_exp_mat, self.syst.atoms[i].p)
      
      self.syst.cell.h = numpy.dot(exp_mat, self.syst.cell.h)

   def vel_step(self):
      #equivalent to P-step in paper
      pass

   def apply_pbc(self):
      """Takes the system and applies periodic boundary conditions to fold the
         particle positions back into the unit cell"""

      for i in range(self.syst.natoms):
         self.syst.atoms[i].q = self.syst.cell.apply_pbc(self.syst.atoms[i])

   def simulation(self, maxcount = 5):
      for i in range(maxcount):
         self.thermo_step()
         self.vel_step()
         self.pos_step()
         self.vel_step()
         self.thermo_step()
      #   print self.syst
      for i in range(self.syst.natoms):
         self.apply_pbc()
      print self.syst
      io_system.print_pdb(self.syst.atoms, self.syst.cell)
