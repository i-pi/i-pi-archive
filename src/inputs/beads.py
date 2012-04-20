"""Deals with creating the beads class.

Classes:
   RestartBeads: Deals with creating the Beads object from a file, and 
      writing the checkpoints.
"""

import numpy as np
from engine.beads import *
from utils.inputvalue import *
import utils.io.io_pdb
from utils.depend import *
from utils.units import *

__all__ = ['RestartBeads']
      
class RestartBeads(Input):
   """Beads input class.

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

   fields={ "natoms"    : (InputValue, {"dtype"     : int,
                                        "default"   : 0,
                                        "help"      : "The number of atoms"}), 
            "nbeads"    : (InputValue, {"dtype"     : int,
                                        "default"   : 0,
                                        "help"      : "The number of beads"}), 
            "q"         : (InputArray, {"dtype"     : float,
                                        "default"   : np.zeros(0),
                                        "help"      : "The positions of the atoms, in the format [x1, y1, z1, x2, ... ]",
                                        "dimension" : "length"}),
            "p"         : (InputArray, {"dtype"     : float,
                                        "default"   : np.zeros(0),
                                        "help"      : "The momenta of the atoms, in the format [px1, py1, pz1, px2, ... ]",
                                        "dimension" : "momentum"}),
            "m"         : (InputArray, {"dtype"     : float, 
                                        "default"   : np.zeros(0),
                                        "help"      : "The masses of the atoms, in the format [m1, m2, ... ]",
                                        "dimension" : "mass"}),
            "names"     : (InputArray, {"dtype"     : str,
                                        "default"   : np.zeros(0, np.dtype('|S6')),
                                        "help"      : "The names of the atoms, in the format [name1, name2, ... ]"}),
            "init_temp" : (InputValue, {"dtype"     : float, 
                                        "default"   : -1.0,
                                        "help"      : "The temperature at which the initial velocity distribution is taken, if applicable."})  }
   
   def __init__(self, beads=None):
      """Initialises RestartBeads.

      Args:
         atoms: An optional Beads object from which to initialise from.
      """

      super(RestartBeads,self).__init__()
      self._optional = True
      if not beads is None:
         self.store(beads)
                       
   def store(self, beads):
      """Takes a Beads instance and stores a minimal representation of it.

      Args:
         beads: A Beads object from which to initialise from.
      """

      super(RestartBeads,self).store(beads)
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

      super(RestartBeads,self).fetch()
      beads = Beads(self.natoms.fetch(),self.nbeads.fetch())
      beads.q = self.q.fetch()
      beads.p = self.p.fetch()      
      beads.m = self.m.fetch()   
      beads.names = self.names.fetch()
      return beads
