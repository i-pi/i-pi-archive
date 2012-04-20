

__all__ = ['RestartBeads']

import numpy as np
from utils.depend import *
from utils.inputvalue import *
from engine.beads import *

class RestartBeads(Input):
   """Beads restart class.

   Handles generating the appropriate beads class from the xml input file,
   and generating the xml checkpoint tags and data from an instance of the 
   object.

   Attributes:
      nbeads: An optional integer giving the number of beads. Defaults to 0.
      natoms: An optional integer giving the number of atoms. Defaults to 0.
      q: An optional array giving the bead positions. Defaults to an empty
         array with no elements.
      p: An optional array giving the bead momenta. Defaults to an empty
         array with no elements.
      m: An optional array giving the bead masses. Defaults to an empty array
         with no elements.
      names: An optional array giving the bead names. Defaults to an empty
         array with no elements.
      init_temp: An optional float giving the kinetic temperature to 
         intialise the bead momenta to.
   """   
   
   fields={ 
      "nbeads" : (InputValue, { "dtype" : int, "default" : 0}), 
      "natoms" : (InputValue, { "dtype" : int, "default" : 0}),  
      "q" : (InputArray, { "dtype" : float, "default" : np.zeros(0)}),
      "p" : (InputArray, { "dtype" : float, "default" : np.zeros(0)}), 
      "m" : (InputArray, { "dtype" : float, "default" : np.zeros(0)}),
      "names" : (InputArray, { "dtype" : str, "default" : np.zeros(0,np.dtype('|S6'))}), 
      "init_temp": (InputValue, { "dtype" : int, "default" : -1.0 }),    }
   
   def __init__(self, beads=None):
      """Initialises RestartBeads.

      Args:
         atoms: An optional Beads object from which to initialise from.
      """

      super(RestartBeads,self).__init__()
      if not beads is None:
         self.store(beads)


   def store(self, beads):
      """Takes a Beads instance and stores a minimal representation of it.

      Args:
         beads: A Beads object from which to initialise from.
      """

      self.natoms.store(beads.natoms)
      self.nbeads.store(beads.nbeads)

      self.q.store(depstrip(beads.q))
      self.p.store(depstrip(beads.p))
      self.m.store(depstrip(beads.m))
      self.names.store(depstrip(beads.names))


   def fetch(self):
      """Creates a beads object.

      Returns:
         A beads object of the appropriate type and with the appropriate
         properties given the attributes of the RestartBeads object.
      """

      beads = Beads(self.natoms.fetch(),self.nbeads.fetch())
      beads.q = self.q.fetch()
      beads.p = self.p.fetch()      
      beads.m = self.m.fetch()   
      beads.names = self.names.fetch()
      return beads

