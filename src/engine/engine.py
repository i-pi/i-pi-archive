import numpy
import math
import random
from io_system import *
from atoms import *
from cell import *

class System(object):
   """
   Represents a simulation cell. 
   Includes the cell parameters, the atoms and the like. """

#__qp holds all the positions and momenta for all the atoms in the simulation
#q and p hold the positions and momenta, respectively.
#P_ext will be the external load.
#The initialisation step now takes a pdc-formatted file for the unit cell and atom positions
#step will eventually call the forces from the external program and then do the propagation step. At the moment we simply take free particle trajectories, to test the theory.

   @classmethod
   def from_pdbfile(cls, filedesc, temp = 1.0):
      atoms, cell, natoms = read_pdb(filedesc)
      cls.natoms = natoms
      cls.temp = temp
      cls.k_Boltz = 1.0

      cls.__qpf=numpy.zeros((3*natoms,3),float) 
      for i in range(natoms):
         cls.__qpf[3*i:3*(i+1),0]=atoms[i][1]

      cls.q=cls.__qpf[:,0]
      cls.f = cls.__qpf[:,2]

      cls.atoms = [ Atom(cls.__qpf[3*i:3*(i+1),:], name = atoms[i][0]) for i in range(natoms) ] #Creates a list of atoms from the __qp array

      cls.P_ext = numpy.zeros((3,3),float)
      cls.cell = Cell(cell, cls.P_ext, temp)

      random.seed(12)
      #cls.__qp[:,1]=numpy.arange(0,3*natoms)*0.01
      for i in range(natoms):
         sigma = math.sqrt(cls.atoms[i].mass * cls.k_Boltz * cls.temp)
         cls.__qpf[3*i,1] = random.gauss(0.0, sigma)
         cls.__qpf[3*i+1,1] = random.gauss(0.0, sigma)
         cls.__qpf[3*i+2,1] = random.gauss(0.0, sigma)
      cls.p=cls.__qpf[:,1]

      cls.pot = 0.0
      cls.kinetic = 0.0
      cls.cell_pot = cls.cell.pot() 
      cls.cell_kinetic = cls.cell.kinetic()
      cls.tot_E = 0.0
      cls.virial = numpy.zeros((3,3),float)
      cls.v_stress = numpy.zeros((3,3),float)
      cls.stress = numpy.zeros((3,3),float)
      return cls()

   @classmethod
   def from_system(cls, syst):
      cls.natoms = syst.natoms
      cls.temp = syst.temp
      cls.k_Boltz = syst.k_Boltz

      cls.q = syst.q
      cls.p = syst.p
      cls.f = syst.f
      cls.__qpf = numpy.zeros((3*cls.natoms,3),float) 
      cls.__qpf[:,0]=cls.q
      cls.__qpf[:,1]=cls.p
      cls.__qpf[:,2]=cls.f

      cls.atoms = syst.atoms

      cls.pot = syst.pot
      cls.kinetic = syst.kinetic
      cls.cell_pot = syst.cell_pot
      cls.cell_kinetic = syst.cell_kinetic
      cls.tot_E = syst.tot_E
      
      cls.virial = syst.virial
      cls.stress = syst.stress
      cls.P_ext = syst.P_ext
      cls.cell = syst.cell
      return cls()

   def __str__(self):
      rstr="ATOMS ("+str(self.natoms)+"):\n\n"
      for i in range(self.natoms): 
         rstr=rstr+"Atom %i:" % (i+1) + "\n"
         rstr=rstr+str(self.atoms[i])+"\n"
      rstr = rstr + "Cell:\n" + str(self.cell)
      rstr = rstr + "\n\nTotal energy = " + str(self.tot_E) + ", potential energy = " + str(self.pot) + ", kinetic energy = " + str(self.kinetic)+ ", cell elastic energy = " + str(self.cell_pot) + ", cell kinetic energy = " + str(self.cell_kinetic)
      return rstr
       
#   def pot(self):
#      """Calculates the total potential energy of the system, including
#         cell strain"""
#
#      pe = 0.0
#      for i in range(self.natoms):
#         pe += self.atoms[i].pot()
#      pe += self.cell.pot()
#      return pe

   def kinetic_update(self):
      """Calculates the total kinetic energy of the system, and the kinetic
         contribution to the stress tensor"""

      self.kinetic = 0.0
      self.v_stress = numpy.zeros((3,3),float)
      for i in range(self.natoms):
         p = self.atoms[i].p
         mass = self.atoms[i].mass
         self.kinetic += numpy.inner(p, p)/(2*mass)
         self.v_stress += numpy.outer(p, p)/(mass*self.cell.V)


   def cell_update(self):
      self.cell_kinetic = self.cell.kinetic()
      self.cell_pot = self.cell.pot()

   def stress_update(self):
      self.stress = self.v_stress + self.virial
      self.stress[2,0] = self.stress[1,0] = self.stress[2,1]

   def tot_E_update(self):
      self.tot_E = self.kinetic + self.pot + self.cell_kinetic + self.cell_pot

#   def tot_E(self):
#      """Calculates the total energy of the system"""
#
#      return self.kinetic() + self.pot()

#   def step(self,dt):
#      """Takes the atom positions, velocities and forces and integrates the 
#         equations of motion forward by a step dt"""
#      self.q+=self.p*dt

#   def apply_pbc(self):
#      """Takes the system and applies periodic boundary conditions to fold the
#         particle positions back into the unit cell"""
#
#      for i in range(self.natoms):
#         self.atoms[i].q = self.cell.apply_pbc(self.atoms[i])
