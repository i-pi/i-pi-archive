import numpy
import math
import random
from utils.depend import *
from io_system import *
from atoms import *
from cell import *
from forces import *

class System(object):
   """Represents a simulation cell. Includes the cell parameters, 
      the atoms and the like. 
      Contains: q = atom positions, p = atom momenta, f = atom forces
      __qpf = (q,p,f), atoms = Atom object list, cell = Cell object,
      temp = temperature, ffield = forcefield object, kin = kinetic energy
      kstress = kinetic stress tensor, pot = potential, vir = virial tensor,
      stress = internal stress tensor, press = internal pressure scalar
      Initialised by: syst = System(ffield)
      ffield = forcefield object, default = forcefield()"""

   #class properties -- initialized here, then possibly re-defined in the init
   __qpf=numpy.zeros(0)
   q=__qpf;  p=__qpf;  f=__qpf;
   atoms=[]; temp=0.0; dt=1.0; 
   cell=Cell()
   
   @classmethod
   def from_pdbfile(cls, filedesc, ffield=forcefield()):   
      """A different initialiser, which takes a pdb format file of a system
         and forms the appropriate atom and cell objects.
         Initialised by: syst = System.from_pdbfile(filedesc, ffield)
         ffield = forcefield object, default = forcefield()
         filedesc = file method, eg. filedesc = open(\"./file.pdb\",\"r\")"""

      atoms, cell, natoms = read_pdb(filedesc)

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

      self.stress = depend(name='stress',func=self.get_stress,deplist=[self.vir, self.kstress])
      self.press = depend(name='press',func=self.get_press,deplist=[self.stress])

   def __str__(self):
      rstr="ATOMS ("+str(len(self.atoms))+"):\n\n"
      for i in range(len(self.atoms)): 
         rstr=rstr+"Atom %i:" % (i+1) + "\n"
         rstr=rstr+str(self.atoms[i])+"\n"
      rstr = rstr + "Cell:\n" + str(self.cell)
      return rstr


   def get_kin(self):
      """Calculates the total kinetic energy of the system,
      by summing the atomic contributions"""

      ke = 0.0
      for at in self.atoms:
         ke += at.kin.get()
      return ke

   def get_stress(self):
      """Calculates the internal stress tensor"""
   
      return self.kstress.get()+self.vir.get()
   
   def get_press(self):
      """Calculates the internal pressure scalar"""

      return numpy.trace(self.stress.get())/3.0
      
   def get_kstress(self):
      """Calculates the kinetic stress tensor"""

      ks=numpy.zeros((3,3),float)
      for at in self.atoms:
         ks += at.kstress.get()
      ks/=self.cell.V.get()
      return ks
      
