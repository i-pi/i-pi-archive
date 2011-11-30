import math, numpy
import io_system
from utils.depend import *

class forcefield(object):
   """Creates an interface between the force, potential and virial calculation
      and the wrapper.
      Contains: _pot = potential, _f = forces, _vir = virial, _ufv = 
      [potential, forces, virial]
      Initialised by ffield = forcefield()"""

   def __init__(self):
      self.ufv = depend(name='ufv', func=self.get_all)

   def bind(self, syst):
      self.syst = syst

#TODO make spring_f an object in syst or RP_syst instead of here.
      self.spring_f = depend(name='spring_f', value=numpy.zeros(3*len(self.syst.atoms)))
      self.ufv.add_dependant(syst.vir)
      self.ufv.add_dependant(syst.f)
      self.ufv.add_dependant(syst.pot)

      syst.q.add_dependant(self.ufv)

      syst.pot._depend__func=self.get_pot
      syst.f._depend__func=self.get_f
      syst.vir._depend__func=self.get_vir

   def get_all(self):
      """Dummy routine where no calculation is done"""

      return [0.0, numpy.zeros(3*len(self.syst.atoms)), numpy.zeros((3,3),float)]

   def get_pot(self):
      """Calls get_all routine of forcefield to update potential"""

      [pot, f, vir] = self.ufv.get()
      return pot

   def get_f(self):
      """Calls get_all routine of forcefield to update potential"""

      [pot, f, vir] = self.ufv.get()
      return f

   def get_vir(self):
      """Calls get_all routine of forcefield to update potential"""

      [pot, f, vir] = self.ufv.get()
      return vir

class pipeforce(forcefield):

   def __init__(self, pars=dict(pipein="pipeforce", pipeout="pipepos")):
      self.ufv = depend(name='ufv', func=self.get_all)
      self.pars=pars

   def get_all(self):
      self.fout=open(self.pars["pipeout"],"w")
      io_system.xml_write(self.syst, self.fout)
      self.fout.close()

      self.fin=open(self.pars["pipein"],"r")
      [pot, f, vir]=io_system.xml_read(self.fin)
      self.fin.close()
      f += self.spring_f.get_array()
      return [pot, f, vir]

class rp_pipeforce(forcefield):

   def __init__(self, pars=dict(pipein="pipeforce", pipeout="pipepos")):
      self.ufv = depend(name='ufv', func=self.get_all)
      self.pars=pars

   def bind(self, rp_syst):
      super(rp_pipeforce,self).bind(rp_syst)
      self.spring_f.set(numpy.zeros_like(rp_syst.f._depend__value))
      rp_syst.temp.add_dependant(self.spring_f)
      for atom in rp_syst.atoms:
         atom.mass.add_dependant(self.spring_f)
      for system in rp_syst.systems:
         for atom in system.atoms:
            atom.q.add_dependant(self.spring_f)
      self.spring_f.add_dependant(self.ufv)
   
      self.ffield_list = []
      depgrp_spring=[]
      for i in range(len(rp_syst.systems)):
         ffield=(pipeforce(pars = self.pars))
         rp_syst.systems[i].ffield = ffield
         ffield.bind(rp_syst.systems[i])
         ffield.spring_f._depend__value=self.spring_f.get_array()[i,:]
         ffield.spring_f.add_depgrp([self.spring_f])
         self.ufv.add_dependant(ffield.ufv)
         depgrp_spring.append(ffield.spring_f)
         self.ffield_list.append(ffield)

      self.spring_f.add_depgrp(depgrp_spring)
      self.spring_f._depend__func=self.get_spring_f
      self.spring_f.taint(taintme=True)

   def get_all(self):
      print "getting..."
      pot = 0.0
      f = numpy.zeros((len(self.syst.systems),3*len(self.syst.systems[0].atoms)))
      vir = numpy.zeros((3,3))
      for i in range(len(self.ffield_list)):
         [sys_pot, sys_f, sys_vir]=self.ffield_list[i].ufv.get()
         pot+=sys_pot;  f[i,:]+=sys_f;  vir+=sys_vir
#Note that spring_f is automatically calculated in sys_f, so we do not need the line f += self.syst.spring_f()
      pot+=self.syst.spring_pot()
      print "got"
      return [pot, f, vir]

   def get_spring_f(self):
      return self.syst.spring_force()

class LJ(forcefield):
   """Creates an interface between the force, potential and virial calculation
      and the wrapper. This uses an internal Lennard-Jones potential
      to do the update.
      Contains: _pot = potential, _f = forces, _vir = virial, _ufv = 
      (potential, forces, virial)
      Initialised by: ffield = pipeforce(dict)
      dict = {\"eps\" : eps_value, \"sigma\" : sigma_value, \"rc\" : cutoff}
      eps_value = energy scale parameter, default = 1.0
      sigma value = length scale parameter, default = 1.0
      cutoff = cutoff radius, default = 2.5"""

   def __init__(self, pars=dict(eps = 1.0, sigma = 1.0, rc = 2.5) ):
      super(LJ,self).__init__() 
      self.eps = pars['eps']; self.sigma = pars['sigma'];  self.rc = pars['rc']
      self.vrc=4*self.eps*((self.sigma/self.rc)**12 - (self.sigma/self.rc)**6)
      self.ufv._depend__func=self.get_all

   def bind(self, syst):
      super(LJ,self).bind(syst)
      self.ufv.func=self.get_all

   def separation(self, atom_i, atom_j):
      """Calculates the vector and scalar separation between two atoms"""

      rij = self.syst.cell.minimum_distance(atom_i, atom_j)
      r = math.sqrt(numpy.dot(rij, rij))
      return r, rij

   def LJ_fdf(self, r):
      """Calculates the force and potential at a given separation"""

      if (r > self.rc):
         return 0.0, 0.0
      else:
         sonr=self.sigma/r; sonr6=sonr**6
         return 4*self.eps*(6/r*sonr6*(2*sonr6-1)), 4*self.eps*(sonr6*(sonr6-1)) - self.vrc
         
   def LJ_fv(self, atom_i, atom_j):
      """Calculates the LJ potential and force between two atoms"""

      r, rij = self.separation(atom_i, atom_j)
      fij, v = self.LJ_fdf(r)
      fij/=r
      fij=rij*fij
      
      return fij, rij, v
      
   def get_all(self):
      """Updates _ufv using an internal LJ potential and force calculator"""

      natoms = len(self.syst.atoms)

      vir = numpy.zeros((3,3),float)
      pot = 0.0
      f = numpy.zeros(3*natoms,float)
   
      for i in range(natoms-1):
         atom_i = self.syst.atoms[i]

         for j in range(i+1, natoms):
            atom_j = self.syst.atoms[j]

            fij, rij, v = self.LJ_fv(atom_i, atom_j)
            f[3*i:3*(i+1)]+=fij
            f[3*j:3*(j+1)]-=fij
            pot += v
            for k in range(3):
               for l in range(k, 3):
                  vir[k,l] += fij[k]*rij[l]
            
      return [pot, f+self.spring_f.get(), vir/self.syst.cell.V.get()]
            
