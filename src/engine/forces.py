import numpy as np
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

