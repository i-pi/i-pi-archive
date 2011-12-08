import numpy as np
import math
from utils.depend import *

class ForceField(dobject):
   """Creates an interface between the force, potential and virial calculation
      and the wrapper.
      Contains: _pot = potential, _f = forces, _vir = virial, _ufv = 
      [potential, forces, virial]
      Initialised by ffield = forcefield()"""

   def __init__(self):
      dset(self,"ufv", depend_value(name="ufv", deps=depend_func(func=self.get_all)) )
      
   def bind(self, atoms, cell):
      self.atoms = atoms
      self.cell = cell
      depget(self,"ufv").add_dependency(depget(self.atoms,"q"))
      dset( self,"pot",depend_value(name="pot", deps=depend_func(func=self.get_pot)) )
      dset(self,"f", depend_array(name="f",   value=np.zeros(atoms.natoms*3, float),  deps=depend_func(func=self.get_f)) )
      dset(self,"vir", depend_array(name="vir", value=np.zeros((3,3),float),            deps=depend_func(func=self.get_vir)) )


   def get_all(self):
      """Dummy routine where no calculation is done"""
      return [0.0, numpy.zeros(3*self.atoms.natoms), numpy.zeros((3,3),float)]
#      return [1.0, np.array(self.atoms.p*self.atoms.q), np.dot(self.cell.ih, self.cell.h)]

   def get_pot(self):
      """Calls get_all routine of forcefield to update potential"""

      [pot, f, vir] = self.ufv
      return pot

   def get_f(self):
      """Calls get_all routine of forcefield to update force"""

      [pot, f, vir] = self.ufv
      return f

   def get_vir(self):
      """Calls get_all routine of forcefield to update virial"""

      [pot, f, vir] = self.ufv
      return vir


class FFLennardJones(ForceField):
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
      super(FFLennardJones,self).__init__() 
      self.eps = pars['eps']; self.sigma = pars['sigma'];  self.rc = pars['rc']
      self.vrc=4*self.eps*((self.sigma/self.rc)**12 - (self.sigma/self.rc)**6)

   def separation(self, atom_i, atom_j):
      """Calculates the vector and scalar separation between two atoms"""
      rij = self.cell.minimum_distance(atom_i, atom_j)
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

      natoms = self.atoms.natoms

      vir = numpy.zeros((3,3),float)
      pot = 0.0
      f = numpy.zeros(3*natoms,float)
   
      for i in range(natoms-1):
         atom_i = self.atoms[i]

         for j in range(i+1, natoms):
            atom_j = self.atoms[j]

            fij, rij, v = self.LJ_fv(atom_i, atom_j)
            f[3*i:3*(i+1)]+=fij
            f[3*j:3*(j+1)]-=fij
            pot += v
            for k in range(3):
               for l in range(k, 3):
                  vir[k,l] += fij[k]*rij[l]
            
      return [pot, f, vir/self.cell.V]
