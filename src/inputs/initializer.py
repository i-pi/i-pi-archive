"""Deals with creating the initiliazer class.

Classes:
   InputInitializer: Initializes the classes that initialize the simulation
      data.
   InputInitFile: Initializes the classes that initialize the simulation data
      from a file.
"""
import numpy as np
from utils.inputvalue import *
from copy import copy
from inputs.beads import InputBeads
from inputs.cell import InputCell
from utils.io import io_xml
import engine.initializer as ei

__all__ = ['InputInitializer', 'InputInitFile']

class InputInitFile(InputValue):
   """Class to handle input from file.

   Attributes:
      format: The format of the file to read data from.
   """

   attribs = copy(InputValue.attribs)
   attribs["cell_units"] = (InputAttribute,{ "dtype" : str, "default": "",
                                       "help": "The units for the cell dimensions." } )
   attribs["format"] =     (InputAttribute,{ "dtype" : str, "default": "xyz", "help": "The input file format. 'xyz' and 'pdb' stand for xyz and pdb input files respectively. 'chk' and 'checkpoint' are both aliases for input from a restart file.", "options": ['xyz', 'pdb', 'chk', 'checkpoint']} )

   default_label = "INITFILE"
   default_help = "This is a simple utility class to deal with initialization from file. Holds all the data needed to open the file and read its contents. The data held between its tags corresponds the the name of the file."

   def __init__(self, help=None, default=None, dtype=None, dimension=None):
      """Initializes InputInitFile.

      Just calls the parent initialize function with appropriate arguments.
      """

      super(InputInitFile,self).__init__(dtype=str, dimension=dimension, default=default, help=help)

   def store(self, iif):
      """Takes a InitFile instance and stores a minimal representation of it.

      Args:
         iif: An input file object.
      """

      super(InputInitFile,self).store(iif.filename, units=iif.units)
      self.format.store(iif.format)
      self.cell_units.store(iif.cell_units)

   def fetch(self):
      """Creates an input file object.

      Returns:
         An input file object.
      """

      return ei.InitFile(filename=super(InputInitFile,self).fetch(), format=self.format.fetch(),
                         units=self.units.fetch(), cell_units=self.cell_units.fetch() )

class InputInitBase(InputValue):
   """Base class to handle input.

   Attributes:
      format: The format of the file to read data from.
   """

   attribs = copy(InputValue.attribs)
   attribs["mode"] =     (InputAttribute,{ "dtype" : str, "default": "array", "help": "The input file format. 'xyz' and 'pdb' stand for xyz and pdb input files respectively. 'chk' and 'checkpoint' are both aliases for input from a restart file.", "options": ['array', 'list', 'value', 'string']} )

   default_label = "INITBASE"
   default_help = "This is the base class for initialization. Initializers for different aspects of the simulation can be inherit for it for the base methods."

   def __init__(self, help=None, default=None, dtype=None, dimension=None):
      """Initializes InputInitFile.

      Just calls the parent initialize function with appropriate arguments.
      """

      super(InputInitBase,self).__init__(dtype=str, dimension=dimension, default=default, help=help)

   def store(self, ibase, mode=None):
      """Takes a InitBase instance and stores a minimal representation of it.

      Args:
         ibase: An input base object.
      """

      if mode is None: mode = ibase.mode
      if mode == "array" or mode == "list":
         value = io_xml.write_list(ibase.value)
      elif mode == "value":
         value = str(ibase.value)
      elif mode == "string":
         value = ibase.value

      super(InputInitBase,self).store(value, units=ibase.units)
      self.mode.store(ibase.mode)

   def fetch(self, mode=None):
      """Creates an input base object.

      Returns:
         An input base object.
      """

      if mode is None : mode = self.mode.fetch()
      value = super(InputInitBase,self).fetch()
      if mode == "list":
         value = io_xml.read_list(value)
      if mode == "array":
         value = io_xml.read_array(np.float, value)
      elif mode == "value":
         value = float(value)
      elif mode == "string":
         value = str(value)  # typically this will be a no-op

      return ei.InitBase(value=value, mode=mode, units=self.units.fetch())

class InputInitVector(InputInitBase):

   attribs = copy(InputInitBase.attribs)
   attribs["index"] =     (InputAttribute,{ "dtype" : int, "default": -1, "help": "The index of the atom of which we are to set the coordinate." } )
   attribs["bead"]  =     (InputAttribute,{ "dtype" : int, "default": -1, "help": "The index of the bead of which we are to set the coordinate." } )
   attribs["mode"] =     (InputAttribute,{ "dtype" : str, "default": "manual", "help": "The input file format. 'xyz' and 'pdb' stand for xyz and pdb input files respectively. 'chk' and 'checkpoint' are both aliases for input from a restart file.", "options": ['manual', 'atom', 'xyz', 'pdb', 'chk']} )


   default_label = "INITVECTOR"
   default_help = "This is a helper class to initialize a vector."

   initclass = None

   def store(self, ipos):

      if ipos.mode == "manual" or ipos.mode == "atom":
         super(InputInitVector,self).store(ipos, "manual")
      else:
         super(InputInitVector,self).store(ipos, "string")

      self.index.store(ipos.index)
      self.bead.store(ipos.bead)

   def check(self):

      if self.mode.fetch() == "atom" and self.index.fetch()<0:
         raise ValueError("Atom initialization requires index to be specified")
      elif self.mode.fetch() != "atom" and self.index.fetch()>=0:
         raise ValueError("Index makes sense only with atom initialization")
      super(InputInitVector,self).check()

   def fetch(self):

      mode = self.mode.fetch()
      if mode == "manual" or mode == "atom":
         ibase=super(InputInitVector,self).fetch("manual")
      else:
         ibase=super(InputInitVector,self).fetch("string")

      return self.initclass(value=ibase.value, mode=mode, units=self.units.fetch(), index=self.index.fetch(), bead=self.bead.fetch())

class InputInitPositions(InputInitVector):

   attribs = copy(InputInitVector.attribs)

   default_label = "INITPOSITIONS"
   default_help = "This is the class to initialize positions."

   initclass = ei.InitPositions

class InputInitVelocities(InputInitVector):

   attribs = copy(InputInitVector.attribs)
   attribs["mode"][1]["options"].append( "thermal" )

   default_label = "INITVELOCITIES"
   default_help = "This is the class to initialize velocities."

   initclass = ei.InitVelocities

class InputInitMomenta(InputInitVector):

   attribs = copy(InputInitVector.attribs)
   attribs["mode"][1]["options"].append( "thermal" )

   default_label = "INITMOMENTA"
   default_help = "This is the class to initialize momenta."

   initclass = ei.InitMomenta

class InputInitCell(InputInitBase):

   attribs = copy(InputInitBase.attribs)
   attribs["mode"][1]["options"]= ['manual', 'pdb', 'chk']

   default_label = "INITCELL"
   default_help = "This is the class to initialize cell."

   initclass = ei.InitMomenta

class InputInitializer(Input):
   """Input class to handle initialization.

   Attributes:
      nbeads: The number of beads to be used in the simulation.
      extra: A list of all the initialize objects read in dynamically from
         the xml input file.
   """

   attribs = { "nbeads"    : (InputAttribute, {"dtype"     : int,
                                        "help"      : "The number of beads. Will override any provision from inside the initializer. A ring polymer contraction scheme is used to scale down the number of beads if required. If instead the number of beads is scaled up, higher normal modes will be initialized to zero."})
            }

   dynamic = {
           "positions"  : (InputInitPositions, { "help" : "Initializes atomic positions" }),
           "velocities" : (InputInitVelocities, { "help" : "Initializes atomic velocities" }),
           "momenta"    : (InputInitMomenta, { "help" : "Initializes atomic momenta" }),

           "beads" : (InputBeads, { "help" : "Initializes the configuration of the path from a Beads object" }),
           "cell" : (InputCell, { "help" : "Initializes the configuration of the cell from a Cell object" }),

           "file" : (InputInitFile, {"help" : "Initializes bead(s) and cell configuration from an external file" }),
           "file_v" : (InputInitFile, {"help" : "Initializes bead(s) velocities from an external file" }),
           "file_p" : (InputInitFile, {"help" : "Initializes bead(s) momenta from an external file" }),
           "resample_v" : (InputValue, {"dimension": "temperature",
                                        "dtype" : float,
                           "help" : "Re-sample the beads' velocities from a Maxwell distribution at the given temperature - or the ensemble temperature if a negative temperature is specified.", "default": -1.0})
           }

   default_help = "Specifies the number of beads, and how the system should be initialized."
   default_label = "INITIALIZER"

   def write(self,  name="", indent=""):
      """Overloads Input write() function so that we never write out
      InputInitializer to restart files.

      Returns:
         An empty string.
      """

      return ""

   def store(self, ii):
      """Takes a Initializer instance and stores a minimal representation of it.

      Args:
         iif: An initializer object.
      """

      self.extra = []

      for (k, el) in ii.queue:
         if k == "file" :
            ip = InputInitFile()
            ip.store(el)
            self.extra.append(("file", ip))
         elif k == "file_v" :
            ip = InputInitFile()
            ip.store(el)
            self.extra.append(("file_v", ip))
         elif k == "file_p" :
            ip = InputInitFile()
            ip.store(el)
            self.extra.append(("file_p", ip))
         elif k == "beads" :
            ip = InputBeads()
            ip.store(el)
            self.extra.append(("beads", ip))
         elif k == "cell" :
            ip = InputCell()
            ip.store(el)
            self.extra.append(("cell", ip))
         elif k == "resample_v" :
            ip = InputValue(dtype=float)
            ip.store(el)
            self.extra.append(("resample_v", ip))

      self.nbeads.store(ii.nbeads)

   def fetch(self):
      """Creates an initializer object.

      Returns:
         An initializer object.
      """

      super(InputInitializer,self).fetch()

      return ei.Initializer(self.nbeads.fetch(), [ (k,v.fetch())  for (k,v) in self.extra ] )

