"""Deals with creating the beads class.

Classes:
   InputBeads: Deals with creating the Beads object from a file, and 
      writing the checkpoints.
"""

import numpy as np
from engine.beads import *
from engine.atoms import Atoms
from utils.inputvalue import *
import utils.io.io_pdb
from utils.depend import *
from utils.units import *
from inputs.atoms import *

__all__ = ['InputBeads']
      
class InputBeads(Input):
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
                                        "help"      : "The number of beads"}), 
            "start_centroid"     : (InputAtoms, {"help"    : "An atoms object from which the centroid coordinates can be initialized", 
                                                 "default" : Atoms(0) }),
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
                                        "help"      : "The temperature at which the initial velocity distribution is taken, if applicable.",
                                        "dimension" : "temperature"})  }

   def write(self,  name="", indent=""):
      """Overloads Input write() function so that nothing is written if
      no beads are present. This will happen if only the classical configuration
      has been specified.

      Returns:
         A string giving the appropriate xml tags for the checkpoint file.
      """

      if self.nbeads._explicit and self.nbeads.fetch() > 0:
         return super(InputBeads,self).write(name=name,indent=indent)
      else:
         return ""
                       
   def store(self, beads):
      """Takes a Beads instance and stores a minimal representation of it.

      Args:
         beads: A Beads object from which to initialise from.
      """

      super(InputBeads,self).store()
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
         properties given the attributes of the InputBeads object.
      """

      super(InputBeads,self).fetch()
      beads = Beads(self.natoms.fetch(),self.nbeads.fetch())
      beads.q = self.q.fetch()
      beads.p = self.p.fetch()  
      beads.m = self.m.fetch()   
      beads.names = self.names.fetch()
      return beads
      
      
   def check(self):
      """Function that deals with optional arguments.

      Deals with deciding which values to initialize from the centroid 
      configurations, and which values to initialize from an explicit array. 
      """

      super(InputBeads,self).check()
      if not (self.start_centroid._explicit or self.q._explicit):
         raise ValueError("Must provide explicit positions or give start_config.")
      if self.start_centroid._explicit:
         # reads the start_config atom tag and created dummy tags for the full bead object
         atoms = self.start_centroid.fetch()
         self.natoms.store(atoms.natoms)
         q = np.zeros((self.nbeads.fetch(),self.natoms.fetch()*3))
         p = np.zeros((self.nbeads.fetch(),self.natoms.fetch()*3))
         m = atoms.m
         names = atoms.names
         for b in range(self.nbeads.fetch()):
            q[b] = atoms.q[:]
            p[b] = atoms.p[:]

         # We can overwrite any of the properties in start_centroid
         # by specifying them in beads.
         if not self.q._explicit:
            self.q.store(q)
         if not self.p._explicit:
            self.p.store(p)
         if not self.m._explicit:
            self.m.store(m)
         if not self.names._explicit:
            self.names.store(names)

      if not 3*self.natoms.fetch()*self.nbeads.fetch() == self.q.fetch().size == self.p.fetch().size == 3*self.nbeads.fetch()*len(self.m.fetch()) == 3*self.nbeads.fetch()*len(self.names.fetch()):
            raise ValueError("Incompatible dimensions of the beads' data arrays.")
      for mass in self.m.fetch():
         if mass <= 0:
            raise ValueError("Unphysical atom mass")
