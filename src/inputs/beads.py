"""Deals with creating the beads class.

Classes:
   InputStartBeads: Deals with starting a simulation from a different
      number of beads.
   InputBeads: Deals with creating the Beads object from a file, and
      writing the checkpoints.
"""

import numpy as np
import math
from engine.beads import *
from engine.atoms import Atoms
from utils.inputvalue import *
import utils.io.io_pdb
from utils.depend import *
from utils.units import *
from inputs.atoms import *

__all__ = ['InputBeads', 'InputStartBeads']

class InputStartBeads(Input):
   """Beads restart input class.

   Generates a stripped down input beads class, that can be used to
   initialize another simulation from a different number of beads.

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
   """

   fields={ "natoms"    : (InputValue, {"dtype"     : int,
                                        "default"   : 0,
                                        "help"      : "The number of atoms."}),
            "nbeads"    : (InputValue, {"dtype"     : int,
                                        "help"      : "The number of beads."}),
            "q"         : (InputArray, {"dtype"     : float,
                                        "default"   : np.zeros(0),
                                        "help"      : "The positions of the beads. In an array of size [nbeads, 3*natoms].",
                                        "dimension" : "length"}),
            "p"         : (InputArray, {"dtype"     : float,
                                        "default"   : np.zeros(0),
                                        "help"      : "The momenta of the beads. In an array of size [nbeads, 3*natoms].",
                                        "dimension" : "momentum"}),
            "m"         : (InputArray, {"dtype"     : float,
                                        "default"   : np.zeros(0),
                                        "help"      : "The masses of the atoms, in the format [m1, m2, ... ].",
                                        "dimension" : "mass"}),
            "names"     : (InputArray, {"dtype"     : str,
                                        "default"   : np.zeros(0, np.dtype('|S6')),
                                        "help"      : "The names of the atoms, in the format [name1, name2, ... ]."})  }

   default_help = "Bead configurations from which to restart a simulation from. Used if the number of beads should be changed after a restart."
   default_label = "RESTART BEADS"

   def write(self,  name="", indent=""):
      """Overloads Input write() function so that we don't restart using this         configuration again.

      Returns:
         An empty string.
      """

      return ""

   def store(self, beads):
      """Takes a Beads instance and stores a minimal representation of it.

      Args:
         beads: A Beads object from which to initialise from.
      """

      super(InputStartBeads,self).store()
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

      super(InputStartBeads,self).fetch()
      beads = Beads(self.natoms.fetch(),self.nbeads.fetch())
      beads.q = self.q.fetch()
      beads.p = self.p.fetch()
      beads.m = self.m.fetch()
      beads.names = self.names.fetch()
      return beads

   def check(self):
      """Function that deals with checking for incorrect input."""

      super(InputStartBeads,self).check()
      if not (self.q._explicit):
         raise ValueError("You must provide a way of generating the starting configuration.")
      if not (self.nbeads.fetch(),3*self.natoms.fetch()) == self.q.fetch().shape:
         raise ValueError("q array is the wrong shape in restart beads object.")
      if not (self.nbeads.fetch(),3*self.natoms.fetch()) == self.p.fetch().shape:
         raise ValueError("p array is the wrong shape in restart beads object.")
      if not (self.natoms.fetch(),) == self.m.fetch().shape:
         raise ValueError("m array is the wrong shape in restart beads object.")
      if not (self.natoms.fetch(),) == self.names.fetch().shape:
         raise ValueError("names array is the wrong shape in restart beads object.")

      for mass in self.m.fetch():
         if mass <= 0:
            raise ValueError("Unphysical atom mass")

class InputBeads(InputStartBeads):
   """Beads input class.

   Handles generating the appropriate beads class from the xml input file,
   and generating the xml checkpoint tags and data from an instance of the
   object.

   Attributes:
      nbeads: An optional integer giving the number of beads. Defaults to 0.
      natoms: An optional integer giving the number of atoms. Defaults to 0.
      start_centroid: An atoms object to initialize the centroid postions from.
      start_beads: A beads object to initialize the normal mode
         coordinates from.
      q: An optional array giving the bead positions. Defaults to an empty
         array with no elements.
      p: An optional array giving the bead momenta. Defaults to an empty
         array with no elements.
      m: An optional array giving the bead masses. Defaults to an empty array
         with no elements.
      names: An optional array giving the bead names. Defaults to an empty
         array with no elements.
   """

   fields={ "natoms"    : (InputValue, {"dtype"     : int,
                                        "default"   : 0,
                                        "help"      : "The number of atoms."}),
            "nbeads"    : (InputValue, {"dtype"     : int,
                                        "help"      : "The number of beads."}),
            "start_centroid"     : (InputAtoms, {"help"    : "An atoms object from which the centroid coordinates can be initialized. Any parameters given here can be overwritten by specifying them explicitly.",
                                                 "default" : Atoms(0) }),
            "start_beads"     : (InputStartBeads, {"help"    : "An beads object from which the lower frequency normal mode coordinates can be initialized. Any parameters given here can be overwritten by specifying them explicitly.",
                                                 "default" : Beads(0, 1) }),
            "q"         : (InputArray, {"dtype"     : float,
                                        "default"   : np.zeros(0),
                                        "help"      : "The positions of the beads. In an array of size [nbeads, 3*natoms].",
                                        "dimension" : "length"}),
            "p"         : (InputArray, {"dtype"     : float,
                                        "default"   : np.zeros(0),
                                        "help"      : "The momenta of the beads. In an array of size [nbeads, 3*natoms].",
                                        "dimension" : "momentum"}),
            "m"         : (InputArray, {"dtype"     : float,
                                        "default"   : np.zeros(0),
                                        "help"      : "The masses of the atoms, in the format [m1, m2, ... ].",
                                        "dimension" : "mass"}),
            "names"     : (InputArray, {"dtype"     : str,
                                        "default"   : np.zeros(0, np.dtype('|S6')),
                                        "help"      : "The names of the atoms, in the format [name1, name2, ... ]."})  }

   default_help = "Deals with the configurations of path integral simulations."
   default_label = "BEADS"

   def write(self,  name="", indent=""):
      """Overloads Input write() function so that nothing is written if
      no beads are present. This will happen if only the classical configuration
      has been specified.

      Returns:
         A string giving the appropriate xml tags for the checkpoint file.
      """

      if self.nbeads._explicit and self.nbeads.fetch() > 0:
         return super(InputStartBeads,self).write(name=name,indent=indent)
      else:
         return ""

   def check(self):
      """Function that deals with optional arguments.

      Deals with deciding which values to initialize from the centroid
      configurations, and which values to initialize from an explicit array.
      """

      if not (self.start_centroid._explicit or self.q._explicit or self.start_beads._explicit):
         raise ValueError("You must provide a way of generating the starting configuration.")

      #!TODO refurbish this way of setting up the necklace
      #~ if self.start_beads._explicit:
         #~ beads = self.start_beads.fetch()
         #~ self.natoms.store(beads.natoms)
         #~ dbeads = Beads(beads.natoms, self.nbeads.fetch())
         #~ for b in range(beads.nbeads):
            #~ dbeads.qnm[b] = depstrip(beads.qnm[b])*math.sqrt(dbeads.nbeads/float(beads.nbeads))
            #~ dbeads.pnm[b] = depstrip(beads.pnm[b])*math.sqrt(dbeads.nbeads/float(beads.nbeads))
         #~ dbeads.qnm[beads.nbeads:,:] = 0.0
         #~ dbeads.pnm[beads.nbeads:,:] = 0.0

         # We can overwrite any of the properties in start_beads
         # by specifying them in beads.
         if not self.q._explicit:
            self.q.store(depstrip(dbeads.q))
         if not self.p._explicit:
            self.p.store(depstrip(dbeads.p))
         if not self.m._explicit:
            self.m.store(depstrip(beads.m))
         if not self.names._explicit:
            self.names.store(depstrip(beads.names))

      elif self.start_centroid._explicit:
         atoms = self.start_centroid.fetch()
         self.natoms.store(atoms.natoms)
         dbeads = Beads(atoms.natoms, self.nbeads.fetch())

         for b in range(self.nbeads.fetch()):
            dbeads.q[b,:]=depstrip(atoms.q)
            dbeads.p[b,:]=depstrip(atoms.p)

         # We can overwrite any of the properties in start_centroid
         # by specifying them in beads.
         if not self.q._explicit:
            self.q.store(depstrip(dbeads.q))
         if not self.p._explicit:
            self.p.store(depstrip(dbeads.p))
         if not self.m._explicit:
            self.m.store(depstrip(atoms.m))
         if not self.names._explicit:
            self.names.store(depstrip(atoms.names))

      super(InputBeads,self).check()
