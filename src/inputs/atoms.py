"""Deals with creating the atoms class.

Generates an atoms class either from a set of positions and momenta, or from 
a configuration file. This class is only used if no beads tag is present in
the xml file.

Classes:
   RestartAtoms: Deals with creating the Atoms object from a file, and 
      writing the checkpoints.
"""

import numpy as np
from engine.atoms import *
from utils.inputvalue import *
import utils.io.io_pdb
from utils.depend import *
from utils.units import *

__all__ = ['RestartAtoms']
      
class RestartAtoms(Input):
   """Atoms input class.

   Handles generating the appropriate atoms class from the xml input file,
   and generating the xml checkpoint tags and data from an instance of the 
   object.

   Attributes:
      natoms: An optional integer giving the number of atoms. Defaults to 0.
      q: An optional array giving the atom positions. Defaults to an empty
         array with no elements.
      p: An optional array giving the atom momenta. Defaults to an empty
         array with no elements.
      m: An optional array giving the atom masses. Defaults to an empty
         array with no elements.
      names: An optional array giving the atom names. Defaults to an empty
         array with no elements.
      from_file: An optional string giving a pdb format file with the atom
         positions. Defaults to ''.
      init_temp: An optional float giving the kinetic temperature to 
         initialise the atom momenta to.
   """

   fields={ "natoms"    : (InputValue, {"dtype"     : int,
                                        "default"   : 0,
                                        "help"      : "The number of atoms" }), 
            "q"         : (InputArray, {"dtype"     : float,
                                        "default"   : np.zeros(0),
                                        "help"      : "The positions of the atoms, in the format [x1, y1, z1, x2, ... ]",
                                        "dimension" : "length" }),
            "p"         : (InputArray, {"dtype"     : float,
                                        "default"   : np.zeros(0),
                                        "help"      : "The momenta of the atoms, in the format [px1, py1, pz1, px2, ... ]",
                                        "dimension" : "momentum" }),
            "m"         : (InputArray, {"dtype"     : float, 
                                        "default"   : np.zeros(0),
                                        "help"      : "The masses of the atoms, in the format [m1, m2, ... ]",
                                        "dimension" : "mass" }),
            "names"     : (InputArray, {"dtype"     : str,
                                        "default"   : np.zeros(0, np.dtype('|S6')),
                                        "help"      : "The names of the atoms, in the format [name1, name2, ... ]" }),
            "from_file" : (InputValue, {"dtype"     : str, 
                                        "default"   : "", 
                                        "help"      : "Gives the name of the file from which the configurations are taken, if present." }),
            "file_units": (InputValue, {"dtype"     : str,
                                        "default"   : "",
                                        "help"      : "The units in which the lengths in the configuration file are given.",
                                        "options"   : ['', 'atomic_unit', 'angstrom', 'nanometer'] }),
            "init_temp" : (InputValue, {"dtype"     : float, 
                                        "default"   : -1.0,
                                        "help"      : "The temperature at which the initial velocity distribution is taken, if applicable.",
                                        "dimension" : "temperature"})  }
       
   def __init__(self, atoms=None, filename=""):
      """Initialises RestartAtoms.

      Args:
         atoms: An optional Atoms object from which to initialise from.
         filename: An optional string giving a filename to take the atom 
            positions from. Defaults to ''.
      """

      super(RestartAtoms,self).__init__()
      self._optional = True
      if not atoms is None:
         self.store(atoms, filename="")
                       
   def store(self, atoms, filename=""):
      """Takes an Atoms instance and stores a minimal representation of it.

      Args:
         atoms: An Atoms object from which to initialise from.
         filename: An optional string giving a filename to take the atom 
            positions from. Defaults to ''.
      """

      super(RestartAtoms,self).store(atoms)
      self.natoms.store(atoms.natoms)
      self.q.store(depstrip(atoms.q))
      self.p.store(depstrip(atoms.p))
      self.m.store(depstrip(atoms.m))
      self.names.store(depstrip(atoms.names))
      self.from_file.store(filename)
      
   def fetch(self):
      """Creates an atoms object.

      Returns:
         An atoms object of the appropriate type and with the appropriate
         properties given the attributes of the RestartAtoms object.
      """

      super(RestartAtoms,self).fetch()
      atoms = Atoms(self.natoms.fetch())
      atoms.q = self.q.fetch()
      atoms.p = self.p.fetch() 
      atoms.m = self.m.fetch()   
      atoms.names = self.names.fetch()
      return atoms
   
   def write(self,  name="", indent=""):
      """Overloads Restart write() function so that nothing is written if
      no atoms are present.

      Returns:
         A string giving the appropriate xml tags for the checkpoint file.
      """

      if self.natoms.fetch() > 0:
         return super(RestartAtoms,self).write(name=name,indent=indent)
      else:
         return ""
      
   
   def check(self): 
      """Function that deals with optional arguments.

      Deals with the init_temp and from_file arguments, and uses them to 
      intialise some of the atoms parameters depending on which ones have
      been specified explicitly.
      """

      super(RestartAtoms,self).check()
      if self.from_file.fetch() != "":
         myatoms, mycell = utils.io.io_pdb.read_pdb(open(self.from_file.fetch(),"r"))
         myatoms.q *= UnitMap["length"][self.file_units.fetch()]
         if len(self.q.fetch()) == 0:
            self.q.store(depstrip(myatoms.q))
         if len(self.p.fetch()) == 0:
            self.p.store(depstrip(myatoms.p))
         if len(self.m.fetch()) == 0:
            self.m.store(depstrip(myatoms.m))
         if len(self.names.fetch()) == 0:
            self.names.store(depstrip(myatoms.names))
         if self.natoms.fetch() == 0:
            self.natoms.store(myatoms.natoms)
