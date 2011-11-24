import numpy, math
import engine, forces, cell, io_system, atoms, cell
from utils import units
from utils.depend import *
from utils.mass_list import *

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

      cls.atoms = [ Necklace(cls.__qpf[:,3*i:3*(i+1),:], name = atom_list[i][0], mass = mlist.masses[atom_list[i][0]]) for i in range(natoms) ]
    
      cls.cell = Cell_necklace(cls.systems)

      return cls(ffield = ffield, temp = temp)

   def __init__(self, ffield = forces.forcefield(), temp = 1.0):
      self.temp = depend(name='temp', value=temp)
      self.betan = depend(name='betan', func=self.get_betan, deplist=[self.temp])
      self.omegan = depend(name='omegan', func=self.get_omegan, deplist=[self.betan])

      self.q = depend(value=self.__qpf[:,:,0], name='q')
      self.p = depend(value=self.__qpf[:,:,1], name='p')
      self.f = depend(value=self.__qpf[:,:,2], name='f')
      self.kin = depend(func=self.get_kin, name = 'kin')
      self.pot = depend(name='pot')
      self.vir = depend(name='vir')

      self.ffield = ffield
      self.ffield.bind(self)

      depgrp_q = []
      for atom in self.atoms:
         atom.q.add_depgrp([self.q])
         depgrp_q.append(atom.q)
      self.q.add_depgrp(depgrp_q)

      depgrps = dict();
      for what in ['q', 'p', 'f']:
         depgrps[what]=[]

      for syst in self.systems:
         syst.kin.add_dependant(self.kin)
         syst.pot.add_dependant(self.pot)
         syst.vir.add_dependant(self.vir)

         for what in ['q', 'p', 'f']:
            depgrps[what].append(getattr(syst,what))
            getattr(syst,what).add_depgrp([getattr(self,what)])

         depgrps_atom = dict();
         for what in ['q', 'p', 'f']:
            depgrps_atom[what]=[]

         for atom in syst.atoms:
            for what in ['q', 'p', 'f']:
               depgrps_atom[what].append(getattr(atom,what))
               getattr(atom,what).add_depgrp([getattr(self,what)])

         for what in ['q', 'p', 'f']:
            getattr(self,what).add_depgrp(depgrps_atom[what])

      for what in ['q', 'p', 'f']:
         getattr(self,what).add_depgrp(depgrps[what])

      self.pot_estimator = depend(func=self.get_pot_estimator, name='pot_estimator', deplist=[self.pot])
      self.kin_estimator = depend(func=self.get_kin_estimator, name='kin_estimator', deplist=[self.betan])
      for syst in self.systems:
         for atom in syst.atoms:
            atom.q.add_dependant(self.kin_estimator)
            atom.f.add_dependant(self.kin_estimator)
      for atom in self.atoms:
         atom.centroid.add_dependant(self.kin_estimator)

      self.trans_mat = depend(value=numpy.zeros((len(self.systems),len(self.systems))), name='trans_mat', func=self.get_trans_mat)
      self.trans_mat.taint(taintme=True)
      self.n_frequencies = depend(value=numpy.zeros(len(self.systems)), name='n_frequencies', func=self.get_n_frequencies, deplist=[self.omegan])
      self.n_frequencies.taint(taintme=True)

   def init_atom_velocities(self, temp=None):
      for syst in self.systems:
         syst.init_atom_velocities(temp=self.temp.get()*len(self.systems))

   def init_cell_velocities(self, temp=None):
      for syst in self.systems:
         syst.init_cell_velocities(temp=self.temp.get()*len(self.systems))

   def get_trans_mat(self):
      nbeads = len(self.systems)
      s_nbeads = math.sqrt(nbeads)
      C = numpy.zeros((nbeads, nbeads))

      C[0,:]=1.0/s_nbeads
      for i in range((nbeads-1)/2):
         for j in range(nbeads):
            C[2*i+1,j] = math.sqrt(2)/s_nbeads*math.sin(2*math.pi*(i+1)*(j+1)/nbeads)
            C[2*i+2,j] = math.sqrt(2)/s_nbeads*math.cos(2*math.pi*(i+1)*(j+1)/nbeads)
         
      if nbeads%2 == 0:
         for j in range(nbeads):
            C[nbeads-1,j] = 1.0/s_nbeads*(-1)**(j+1)

      return C

   def get_n_frequencies(self):
      nbeads = len(self.systems)
      omega = numpy.zeros(nbeads)
      
      omega[0] = 0.0
      for i in range((nbeads-1)/2):
         omega[2*i+1] = omega[2*i+2] = 2*self.omegan.get()*math.sin((i+1)*math.pi/nbeads)
      if nbeads%2 == 0:
         omega[nbeads-1] = 2.0*self.omegan.get()

      return omega

   def get_omegan(self):
      return 1.0/(units.hbar*self.betan.get())

   def get_betan(self):
      return 1.0/(self.temp.get()*units.kb*len(self.systems))

   def get_kin(self):
      ke = 0.0
      for syst in self.systems:
         ke += syst.kin.get()
      return ke

   def get_pot_estimator(self):
      return (self.pot.get()-self.spring_pot())/float(len(self.systems))

   def get_kin_estimator(self):
      f = self.spring_force()
      kin = 0.0
      for j in range(len(self.systems)):
         for i in range(len(self.atoms)):
            dr = self.systems[j].atoms[i].q.get_array()-self.atoms[i].centroid.get_array()
            df = self.systems[j].atoms[i].f.get_array()-f[j,3*i:3*(i+1)]
            kin -= 0.5*numpy.dot(dr,df)

      kin += 3.0*len(self.atoms)/(2.0*self.betan.get())
      kin /= len(self.systems)
      return kin

   def spring_pot(self):
      pot = 0.0
      natoms = len(self.atoms)
      nbeads = len(self.systems)

      for i in range(natoms):
         spring_const = self.atoms[i].mass.get()/(self.betan.get()*units.hbar)**2
         for j in range(nbeads):
            rij = self.systems[j].atoms[i].q.get_array() - self.systems[(j+1)%nbeads].atoms[i].q.get_array()
            r_2 = numpy.dot(rij, rij)
            pot += 0.5*spring_const*r_2

      return pot

   def spring_force(self):
      natoms = len(self.atoms)
      nbeads = len(self.systems)
      f = numpy.zeros((nbeads, natoms*3))

      for i in range(natoms):
         spring_const = self.atoms[i].mass.get()/(self.betan.get()*units.hbar)**2
         rij_minus = self.systems[0].atoms[i].q.get_array() - self.systems[nbeads-1].atoms[i].q.get_array()

         for j in range(nbeads):
            up_index = (j+1)%nbeads

            rij_plus = self.systems[j].atoms[i].q.get_array() - self.systems[up_index].atoms[i].q.get_array()

            f[j,3*i:3*(i+1)] -= spring_const*(rij_plus+rij_minus)

            rij_minus = -rij_plus

      return f

class Necklace:
   def __init__(self, qpf_slice, name="X", mass=1.0):
      self.q = depend(value=qpf_slice[:,:,0], name='q')
      #self.p = self.depend(value=qpf_slice[:,:,1], name='p')
      #self.f = self.depend(value=qpf_slice[:,:,2], name='f')
      self.centroid = depend(value=numpy.zeros(3), name='centroid', func=self.get_centroid)
      self.q.add_dependant(self.centroid)
      self.centroid.taint(taintme=True)

      self.mass = depend(value=mass, name='mass')
      self.name = depend(value=name, name='name')

   def get_centroid(self):
      nbeads = len(self.q.get_array()[:,0])
      q_bar = numpy.zeros(3)

      for i in range(nbeads):
         q_bar += self.q.get_array()[i,:]
      q_bar /= nbeads

      return q_bar

class Cell_necklace:
   def __init__(self, systems):
      self.w = systems[0].cell.w
      nbeads = len(systems)

      p = numpy.zeros((nbeads, 3,3))
      for i in range(nbeads):
         systems[i].cell.p._depend__value = p[i,:,:]

      depgrp_p = []
      self.p = depend(value=p, name='p')
      for i in range(nbeads):
         systems[i].cell.p.add_depgrp([self.p])
         depgrp_p.append(systems[i].cell.p)
      self.p.add_depgrp(depgrp_p)
      #Note that this re-initialises the cell momenta to zero, if they have already been initialised
