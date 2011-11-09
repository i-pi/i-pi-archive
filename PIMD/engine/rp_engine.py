import numpy, math
import engine, forces, cell, io_system, atoms, cell
from utils import units
from utils.depend import *

class RP_sys(engine.System):
   """A class of multiple copies of a system for use in PIMD"""

   systems=[]
   __qpf = numpy.zeros(0)

   @classmethod
   def from_pdbfile(cls, filedesc, ffield = forces.forcefield(), nbeads = 8, temp = 1.0):
      """A different initialiser, which takes a pdb formatted file of a system
         and forms the appropriate atom and cell objects.
         Initialised by: 
         sys = RP_sys.from_pdbfile(filedesc, ffield, nbeads, temp)

         ffield = forcefield object, default = forcefield()
         filedesc = file method, eg. filedesc = open(\".file.pdb\",\"r\")
         nbeads = number of copies of the system, default = 8
         temp = temperature, default = 1.0"""
   
      cls.systems = []
      atom_list, cell_list, natoms = io_system.read_pdb(filedesc)
      filedesc.seek(0)
      cls.__qpf = numpy.zeros((nbeads, 3*natoms, 3))

      for i in range(nbeads):
         cls.systems.append(engine.System.from_pdbfile(filedesc, ffield=forces.forcefield(), qpf_slice = cls.__qpf[i,:,:]))
         filedesc.seek(0)

      cls.atoms = [ atoms.Atom(numpy.zeros((3,3)), name = atom_list[i][0]) for i in range(natoms) ]
      cls.cell = cell.Cell.fromSidesAngles(cell_list)

      return cls(ffield = ffield, temp = temp)

   def __init__(self, ffield = forces.forcefield(), temp = 1.0):
      self.temp = depend(name='temp', value=temp)
      self.betan = depend(name='betan', func=self.get_betan, deplist=[self.temp])

      self.q = depend(value=self.__qpf[:,:,0], name='q')
      self.p = depend(value=self.__qpf[:,:,1], name='p')
      self.f = depend(value=self.__qpf[:,:,2], name='f')
      self.kin = depend(func=self.get_kin, name = 'kin')
      self.pot = depend(name='pot')
      self.vir = depend(name='vir')

      self.ffield = ffield
      self.ffield.bind(self)

      depgrps = dict();
      for what in ['q', 'p', 'f']:
         depgrps[what]=[]

      for syst in self.systems:
         syst.kin.add_dependant(self.kin)
         syst.pot.add_dependant(self.pot)
         syst.vir.add_dependant(self.vir)

         for what in ['q', 'p', 'f']:
            depgrps[what].append(getattr(syst,what))
            getattr(syst,what).add_dependant(getattr(self,what))

         depgrps_atom = dict();
         for what in ['q', 'p', 'f']:
            depgrps_atom[what]=[]

         for atom in syst.atoms:
            for what in ['q', 'p', 'f']:
               depgrps_atom[what].append(getattr(atom,what))
               getattr(atom,what).add_dependant(getattr(self,what))

         for what in ['q', 'p', 'f']:
            getattr(self,what).add_depgrp(depgrps_atom[what])

      for what in ['q', 'p', 'f']:
         getattr(self,what).add_depgrp(depgrps[what])

   def init_atom_velocities(self, temp=None):
      for syst in self.systems:
         syst.init_atom_velocities(temp=self.temp.get()*len(self.systems))

   def init_cell_velocities(self, temp=None):
      for syst in self.systems:
         syst.init_cell_velocities(temp=self.temp.get()*len(self.systems))

   def get_betan(self):
      return 1.0/(self.temp.get()*units.kb*len(self.systems))

   def get_kin(self):
      ke = 0.0
      for syst in self.systems:
         ke += syst.kin.get()
#      return ke/float(len(self.systems))
      return ke

   def spring_pot(self):
      pot = 0.0
      natoms = len(self.atoms)
      nbeads = len(self.systems)

      for i in range(natoms):
         spring_const = self.atoms[i].mass.get()/(self.betan.get()*units.hbar)**2

         for j in range(nbeads):
            rij = self.systems[j].cell.minimum_distance(self.systems[j].atoms[i], self.systems[(j+1)%nbeads].atoms[i])
            r_2 = numpy.dot(rij, rij)
            pot += 0.5*spring_const*r_2

      return pot

   def spring_force(self):
      natoms = len(self.atoms)
      nbeads = len(self.systems)
      f = numpy.zeros(nbeads, natoms*3)

      for i in range(natoms):
         spring_const = self.atoms[i].mass.get()/(self.betan.get()*units.hbar)**2
         rij_minus = self.systems[0].cell.minimum_distance(self.systems[0].atoms[i], self.systems[natoms].atoms[i])

         for j in range(nbeads):
            up_index = (j+1)%nbeads

            rij_plus = self.systems[j].cell.minimum_distance(self.systems[j].atoms[i], self.systems[up_index].atoms[i])

            f[j,3*i:3(i+1)] -= spring_const*(rij_plus+rij_minus)

            rij_minus = self.systems[up_index].cell.minimum_distance(self.systems[up_index].atoms[i], self.systems[j].atoms[i])

      return f

