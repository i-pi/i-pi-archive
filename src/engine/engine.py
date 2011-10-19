import numpy
import math
import random
from utils.depend import *
from io_system import *
from atoms import *
from cell import *
from forces import *

class System(object):
   """
   Represents a simulation cell. 
   Includes the cell parameters, the atoms and the like. """

#__qp holds all the positions and momenta for all the atoms in the simulation
#q and p hold the positions and momenta, respectively.
#P_ext will be the external load.
#The initialisation step now takes a pdc-formatted file for the unit cell and atom positions
#step will eventually call the forces from the external program and then do the propagation step. At the moment we simply take free particle trajectories, to test the theory.

   #class properties -- initialized here, then possibly re-defined in the init
   __qpf=numpy.zeros(0)
   q=__qpf;  p=__qpf;  f=__qpf;
   atoms=[]; temp=0.0; dt=1.0; 
   cell=Cell()
   
#   @classmethod
#   def from_system(cls, syst):
#      cls.temp = syst.temp.get()
#      cls.dt = sys.dt.get()

#      cls.__qpf = numpy.copy(syst._system__qpf)

#      natoms = len(syst.atoms)
#      for i in range(natoms):
#         cls.__qpf[3*i:3*(i+1),0]=atoms[i][1]
#      cls.atoms = [ Atom(cls.__qpf[3*i:3*(i+1),:], name = syst.atoms[i].name.get()) for i in range(natoms) ]
#      cls.cell = Cell(h=numpy.copy(syst.cell.h.get()))

#      # must decide what to do with these
##      cls.pot = syst.pot
##      cls.kinetic = syst.kinetic
##      cls.cell_pot = syst.cell_pot
##      cls.cell_kinetic = syst.cell_kinetic
##      cls.tot_E = syst.tot_E
##      
##      cls.virial = syst.virial
##      cls.stress = syst.stress
##      cls.P_ext = syst.P_ext
##      cls.cell = syst.cell
#      return cls()
      
   @classmethod
   def from_pdbfile(cls, filedesc, ffield=forcefield()):   
      atoms, cell, natoms = read_pdb(filedesc)

#      atoms[0][1] = numpy.array([0.5, 0.5, 0.5])
#      atoms[1][1] = numpy.array([0.5, 0.5, 0.5+0.15*2**(1.0/6.0) + 0.1])

      cls.__qpf=numpy.zeros((3*natoms,3),float) 
      for i in range(natoms):
         cls.__qpf[3*i:3*(i+1),0]=atoms[i][1]
      cls.atoms = [ Atom(cls.__qpf[3*i:3*(i+1),:], name = atoms[i][0]) for i in range(natoms) ] #Creates a list of atoms from the __qp array
      cls.cell = Cell.fromSidesAngles(cell)

      return cls(ffield=ffield);
   
   def __init__(self, ffield=forcefield):
      self.q = depend(value=self.__qpf[:,0],name='q')
      self.p = depend(value=self.__qpf[:,1],name='p')
      self.f = depend(value=self.__qpf[:,2],name='f')
      self.kin = depend(func=self.get_kin, name = 'kin')
      self.kstress = depend(name='kstress',func=self.get_kstress,deplist=[self.cell.V])
                  
      depgrps=dict();
      for what in [ 'q', 'p', 'f' ]: depgrps[what]=[];
      
      for atom in self.atoms:
         atom.kin.add_dependant(self.kin)
         atom.kstress.add_dependant(self.kstress)         
         # the global q,p,f and the atomic q,p,f must be synchronized.
         # note that modifying the global vectors in any position will taint ALL the atomic vectors
         for what in [ 'q', 'p', 'f' ]: 
            depgrps[what].append(getattr(atom,what))
            getattr(atom,what).add_dependant(getattr(self,what))
      
      for what in [ 'q', 'p', 'f' ]:  getattr(self,what).add_depgrp(depgrps[what])
      
      self.pot = depend(value=0.0,name='pot')
      self.vir = depend(value=numpy.zeros((3,3),float),name='vir')

      self.ffield=ffield
      self.ffield.bind(cell=self.cell, atoms=self.atoms, pot=self.pot, f=self.f, vir=self.vir)      

      # to be defined somewhere      
      self.stress = depend(name='stress',func=self.get_stress,deplist=[self.vir, self.kstress])
      

   def __str__(self):
      rstr="ATOMS ("+str(len(self.atoms))+"):\n\n"
      for i in range(len(self.atoms)): 
         rstr=rstr+"Atom %i:" % (i+1) + "\n"
         rstr=rstr+str(self.atoms[i])+"\n"
      rstr = rstr + "Cell:\n" + str(self.cell)
      #rstr = rstr + "\n\nTotal energy = " + str(self.pot+self.kinetic+self.thermo.econs) + ", potential energy = " + str(self.pot) + ", kinetic energy = " + str(self.kinetic)+ ", cell elastic energy = " + str(self.cell_pot) + ", cell kinetic energy = " + str(self.cell_kinetic)
      return rstr


   def get_kin(self):
      """Calculates the total kinetic energy of the system,
      by summing the atomic contributions"""
      #print " [ upd. kin ]",
      ke = 0.0
      for at in self.atoms:
         ke += at.kin.get()
      return ke

   def get_stress(self):
      return self.kstress.get()+self.vir.get()
   
   def get_kstress(self):
      ks=numpy.zeros((3,3),float)
      for at in self.atoms:
         ks += at.kstress.get()
      ks/=self.cell.V.get()
      return ks
        
      
      
#   def get_kinetic(self):
#      """Calculates the total kinetic energy of the system, and the kinetic
#         contribution to the stress tensor"""

#      self.kinetic = 0.0
#      self.v_stress = numpy.zeros((3,3),float)
#      for i in range(self.natoms):
#         p = self.atoms[i].p.get()
#         mass = self.atoms[i].mass.get()
#         self.kinetic += numpy.inner(p, p)/(2*mass)
#         self.v_stress += numpy.outer(p, p)/(mass*self.cell.V)

   def cell_update(self):
      self.cell_kinetic = self.cell.kinetic()
      self.cell_pot = self.cell.pot()

   def stress_update(self):
      self.stress = self.v_stress + self.virial
      self.stress[2,0] = self.stress[1,0] = self.stress[2,1] = 0.0

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
