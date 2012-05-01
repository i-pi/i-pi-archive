"""Deals with creating the atoms class.

Generates an atoms class either from a set of positions and momenta, or from 
a configuration file. This class is only used if no beads tag is present in
the xml file.

Classes:
   InputAtoms: Deals with creating the Atoms object from a file, and 
      writing the checkpoints.
"""

import numpy as np
from engine.atoms import *
from utils.inputvalue import *
import utils.io.io_pdb, utils.io.io_xyz
from utils.depend import *
from utils.units import *

__all__ = ['InputAtoms']
      
class InputAtoms(Input):
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
      file_units: An optional string giving the length units that the file is
         specified by. Defaults to ''.
   """

   fields={ "natoms"    : (InputValue, {"dtype"     : int,
                                        "default"   : 0,
                                        "help"      : "The number of atoms." }), 
            "q"         : (InputArray, {"dtype"     : float,
                                        "default"   : np.zeros(0),
                                        "help"      : "The positions of the atoms, in the format [x1, y1, z1, x2, ... ].",
                                        "dimension" : "length" }),
            "p"         : (InputArray, {"dtype"     : float,
                                        "default"   : np.zeros(0),
                                        "help"      : "The momenta of the atoms, in the format [px1, py1, pz1, px2, ... ].",
                                        "dimension" : "momentum" }),
            "m"         : (InputArray, {"dtype"     : float, 
                                        "default"   : np.zeros(0),
                                        "help"      : "The masses of the atoms, in the format [m1, m2, ... ].",
                                        "dimension" : "mass" }),
            "names"     : (InputArray, {"dtype"     : str,
                                        "default"   : np.zeros(0, np.dtype('|S6')),
                                        "help"      : "The names of the atoms, in the format [name1, name2, ... ]." }),
            "from_file" : (InputValue, {"dtype"     : str, 
                                        "default"   : "", 
                                        "help"      : "Gives the name of the file from which the configurations are taken, if present. Any value given in this file can be overwritten by specifying it explicitly." }),
            "file_units": (InputValue, {"dtype"     : str,
                                        "default"   : "",
                                        "help"      : "The units in which the lengths in the configuration file are given.",
                                        "options"   : [unit for unit in UnitMap["length"]]})  }

   default_help = "Deals with single replicas of the system or classical simulations."
   default_label = "ATOMS"
       
   def store(self, atoms, filename=""):
      """Takes an Atoms instance and stores a minimal representation of it.

      Args:
         atoms: An Atoms object from which to initialise from.
         filename: An optional string giving a filename to take the atom 
            positions from. Defaults to ''.
      """

      super(InputAtoms,self).store()
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
         properties given the attributes of the InputAtoms object.
      """

      super(InputAtoms,self).fetch()
      atoms = Atoms(self.natoms.fetch())
      atoms.q = self.q.fetch()
      atoms.p = self.p.fetch() 
      atoms.m = self.m.fetch()   
      atoms.names = self.names.fetch()
      return atoms
   
   def write(self,  name="", indent=""):
      """Overloads Input write() function so that nothing is written if
      no atoms are present. This occurs if the beads object has been specified,
      so that the classical atoms object is not initialized.

      Returns:
         A string giving the appropriate xml tags for the checkpoint file.
      """

      if self.natoms.fetch() > 0:
         return super(InputAtoms,self).write(name=name,indent=indent)
      elif self.from_file.fetch() != "":
         rstr = InputValue(dtype=str)
         rstr.store(self.from_file.fetch())
         return indent + "<" + name + ">" + rstr.write("from_file","") + "</" + name + ">\n"
      else:
         return ""
      
   
   def check(self): 
      """Function that deals with optional arguments.

      Deals with the from_file argument, and uses it to 
      intialise some of the atoms parameters depending on which ones have
      been specified explicitly.
      """

      super(InputAtoms,self).check()
      if not (self.from_file._explicit or self.q._explicit):
         raise ValueError("Must provide explicit positions or give from_file.")
      if self.from_file._explicit:
      
         filename=self.from_file.fetch()
         ext=filename[len(filename)-3:]
         if (ext == "pdb"):
            myatoms, mycell = utils.io.io_pdb.read_pdb(open(self.from_file.fetch(),"r"))
         elif (ext == "xyz"):
            myatoms = utils.io.io_pdb.read_xyz(open(self.from_file.fetch(),"r"))
         else:
            raise ValueError("Unrecognized extension for atomic configuration file")
            
         myatoms.q *= UnitMap["length"][self.file_units.fetch()]

         # We can overwrite any of the properties in from_file
         # by specifying them in atoms.
         if not self.q._explicit:
            self.q.store(depstrip(myatoms.q))
         if not self.p._explicit:
            self.p.store(depstrip(myatoms.p))
         if not self.m._explicit:
            self.m.store(depstrip(myatoms.m))
         if not self.names._explicit:
            self.names.store(depstrip(myatoms.names))
         self.natoms.store(myatoms.natoms)

      if not (3*self.natoms.fetch(),) == self.q.fetch().shape:
         raise ValueError("q array is the wrong shape in atoms object.")
      if not (3*self.natoms.fetch(),) == self.p.fetch().shape:
         raise ValueError("p array is the wrong shape in atoms object.")
      if not (self.natoms.fetch(),) == self.m.fetch().shape:
         raise ValueError("m array is the wrong shape in atoms object.")
      if not (self.natoms.fetch(),) == self.names.fetch().shape:
         raise ValueError("names array is the wrong shape in atoms object.")

      for mass in self.m.fetch():
         if mass <= 0:
            raise ValueError("Unphysical atom mass")
