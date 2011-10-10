import numpy, math
import engine, io_system

class NST_ens:

   def __init__(self, filedesc, thermo, temp = 1.0):
      self.dt = thermo.dt*2.0
      self.syst = engine.System(filedesc, temp)
      self.thermo = thermo

   def thermo_step(self):
      for i in range(self.syst.natoms):
         self.thermo.step(self.syst.atoms[i])
     # self.thermo.step(syst.cell)

   def pos_step(self):
      """Takes the atom positions, velocities and forces and integrates the 
         equations of motion forward by a step dt"""

      #equivalent to R-step in paper

      self.syst.q+=self.syst.p*self.dt
      #do something to the cell parameters

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
         print self.syst
      for i in range(self.syst.natoms):
         self.apply_pbc()
      print self.syst
      io_system.print_pdb(self.syst.atoms, self.syst.cell)
