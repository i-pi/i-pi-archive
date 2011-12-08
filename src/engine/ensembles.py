import numpy as np
from utils.depend import *
from utils import units

class Ensemble(dobject): 
   """General ensemble object, with no particle motion
      Contains: syst = System object, containing the atom and cell coordinates, 
      econs = conserved energy quantity.
      Initialised by: ens = ensemble(system)
      system is a System object, containing the atom and cell coordinates"""

   def __init__(self):
      dset(self,"econs",depend_value(name='econs',deps=depend_func(func=self.get_econs)) )
      
   def bind(self, atoms, cell, force):
      self.atoms=atoms
      self.cell=cell
      self.force=force
      depget(self,"econs").add_dependency(depget(self.atoms, "kin"))
      depget(self,"econs").add_dependency(depget(self.force, "pot"))      
      
   def step(self): 
      """Dummy routine which does nothing"""
      pass

   def get_econs(self):
      """Calculates the conserved energy quantity for constant E ensembles"""
      return self.atoms.kin+self.force.pot
      
      
class NVEEnsemble(Ensemble):
   def __init__(self, dt=1.0):
      super(NVEEnsemble,self).__init__()
      dset(self,"dt",depend_value(name='dt',value=dt))
   
   def step(self):
      """Velocity Verlet time step"""

      self.atoms.p += self.force.f * (self.dt*0.5)
      self.atoms.q += self.atoms.p/self.atoms.m3 *self.dt
      self.atoms.p += self.force.f * (self.dt*0.5)

