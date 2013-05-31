"""Deals with creating the initiliazer class.

Classes:
   InputInitializer: Initializes the classes that initialize the simulation
      data.
   InputInitFile: Initializes the classes that initialize the simulation data
      from a file.
"""
import numpy as np
from utils.inputvalue import *
from copy import copy, deepcopy
from inputs.beads import InputBeads
from inputs.cell import InputCell
from utils.io import io_xml
import utils.mathtools as mt
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

   attribs = deepcopy(InputValue.attribs)
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

   attribs = deepcopy(InputInitBase.attribs)
   attribs["index"] =     (InputAttribute,{ "dtype" : int, "default": -1, "help": "The index of the atom of which we are to set the coordinate." } )
   attribs["bead"]  =     (InputAttribute,{ "dtype" : int, "default": -1, "help": "The index of the bead of which we are to set the coordinate." } )
   attribs["mode"] =     (InputAttribute,{ "dtype" : str, "default": "manual", "help": "The input file format. 'xyz' and 'pdb' stand for xyz and pdb input files respectively. 'chk' and 'checkpoint' are both aliases for input from a restart file.", "options": ['manual', 'xyz', 'pdb', 'chk']} )


   default_label = "INITVECTOR"
   default_help = "This is a helper class to initialize a vector."

   _initclass = None
   _storageclass = float

   def store(self, ipos):

      if ipos.mode == "manual":
         if self._storageclass is float:
            super(InputInitVector,self).store(ipos, "array")
         else:
            super(InputInitVector,self).store(ipos, "list")
      else:
         super(InputInitVector,self).store(ipos, "string")

      self.index.store(ipos.index)
      self.bead.store(ipos.bead)

   def fetch(self, initclass=None):

      mode = self.mode.fetch()
      if mode == "manual":
         if self._storageclass is float:
            ibase=super(InputInitVector,self).fetch("array")
         else:
            ibase=super(InputInitVector,self).fetch("list")
      else:
         ibase=super(InputInitVector,self).fetch("string")

      if initclass is None:
         initclass = self._initclass
      return initclass(value=ibase.value, mode=mode, units=self.units.fetch(), index=self.index.fetch(), bead=self.bead.fetch())

class InputInitPositions(InputInitVector):

   attribs = deepcopy(InputInitVector.attribs)

   default_label = "INITPOSITIONS"
   default_help = "This is the class to initialize positions."

   _initclass = ei.InitPositions

class InputInitMomenta(InputInitVector):

   attribs = deepcopy(InputInitVector.attribs)
   attribs["mode"][1]["options"].append( "thermal" )

   default_label = "INITMOMENTA"
   default_help = "This is the class to initialize momenta."

   _initclass = ei.InitMomenta

   def fetch(self):
      if self.mode.fetch() == "thermal":
         return self._initclass(value=float(InputValue.fetch(self)),  mode=self.mode.fetch(), units=self.units.fetch(), index=self.index.fetch(), bead=self.bead.fetch())
      else:
         return super(InputInitMomenta,self).fetch()


class InputInitVelocities(InputInitMomenta):

   attribs = deepcopy(InputInitVector.attribs)
   attribs["mode"][1]["options"].append( "thermal" )

   default_label = "INITVELOCITIES"
   default_help = "This is the class to initialize velocities."

   _initclass = ei.InitVelocities

class InputInitCell(InputInitBase):

   attribs = deepcopy(InputInitBase.attribs)
   attribs["mode"] = (InputAttribute, { "dtype"  : str,
                                        "default": "manual",
                                        "options": ["manual", "pdb", "chk", "abc", "abcABC"],
                                        "help"   : "This decides whether the system box is created from a cell parameter matrix, or from the side lengths and angles between them. If 'mode' is 'manual', then 'cell' takes a 9-elements vector containing the cell matrix (row-major). If 'mode' is 'abcABC', then 'cell' takes an array of 6 floats, the first three being the length of the sides of the system parallelopiped, and the last three being the angles (in degrees) between those sides. Angle A corresponds to the angle between sides b and c, and so on for B and C. If mode is 'abc', then this is the same as ffor 'abcABC', but the cell is assumed to be orthorhombic. 'pdb' and 'chk' read the cell from a PDB or a checkpoint file, respectively."} )

   default_label = "INITCELL"
   default_help = "This is the class to initialize cell."

   def store(self, ipos):
      if ipos.mode == "manual":
         super(InputInitCell,self).store(ipos, "array")
      else:
         super(InputInitCell,self).store(ipos, "string")

   def fetch(self):

      mode = self.mode.fetch()
      if mode == "manual" or mode == "abc" or mode == "abcABC":
         ibase=super(InputInitCell,self).fetch("array")
         h = ibase.value

         if mode == "abc":
            if h.size != 3:
               raise ValueError("If you are initializing cell from cell side lengths you must pass the 'cell' tag an array of 3 floats.")
            else:
               h = mt.abc2h(h[0], h[1], h[2], np.pi/2, np.pi/2, np.pi/2)
         elif mode == "abcABC":
            if h.size != 6:
               raise ValueError("If you are initializing cell from cell side lengths and angles you must pass the 'cell' tag an array of 6 floats.")
            else:
               h = mt.abc2h(h[0], h[1], h[2], h[3]*np.pi/180.0, h[4]*np.pi/180.0, h[5]*np.pi/180.0)

         if h.size != 9:
               raise ValueError("Cell objects must contain a 3x3 matrix describing the cell vectors.")

         h.shape = (9,)
         if not (h[3] == 0.0 and h[6] == 0.0 and h[7] == 0.0):
            warning("Cell vector matrix must be upper triangular, all elements below the diagonal being set to zero.", verbosity.low)
            h[3] = h[6] = h[7] = 0
         ibase.value = h

      else:
         ibase=super(InputInitCell,self).fetch("string")

      return ei.InitCell(value=ibase.value, mode=mode, units=self.units.fetch())

class InputInitMasses(InputInitVector):

   attribs = deepcopy(InputInitVector.attribs)
   attribs["mode"][1]["options"]= ['manual', 'xyz', 'pdb', 'chk']

   default_label = "INITMASSES"
   default_help = "This is the class to initialize atomic masses."

   _initclass = ei.InitMasses

class InputInitLabels(InputInitVector):

   attribs = deepcopy(InputInitVector.attribs)
   attribs["mode"][1]["options"]= ['manual', 'xyz', 'pdb', 'chk']

   default_label = "INITLABELS"
   default_help = "This is the class to initialize atomic labels."

   _storageclass = str
   _initclass = ei.InitMomenta


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
           "masses"    : (InputInitMasses, { "help" : "Initializes atomic masses" }),
           "labels"    : (InputInitLabels, { "help" : "Initializes atomic labels" }),
           "cell" :     (InputInitCell, { "help" : "Initializes the configuration of the cell" }),
           "all" :      (InputInitVector, { "help" : "Initializes everything possible for the given mode" }),

           "beads" : (InputBeads, { "help" : "Initializes the configuration of the path from a Beads object" }),


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

      # TODO: this is broken!
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

      initlist = []
      for (k,v) in self.extra:
         if k == "all":
            mode = v.mode.fetch()
            if mode == "xyz" or mode == "manual" or mode == "pdb" or mode == "chk":
               initlist.append( ( "positions", v.fetch(initclass=ei.InitPositions) ) )
            if mode == "xyz" or mode == "pdb" or mode == "chk":
               rm=v.fetch(initclass=ei.InitMasses); rm.units = ""
               initlist.append( ( "masses",   rm ) )
               initlist.append( ( "labels",   v.fetch(initclass=ei.InitLabels) ) )
            if mode == "pdb" or mode == "chk":
               initlist.append( ( "cell", InputInitCell(base=v).fetch() ) )
            if mode == "chk":
               initlist.append( ( "momenta", InputInitMomenta(base=v).fetch() ) )
         else:
            initlist.append( (k, v.fetch()) )
      print initlist
      return ei.Initializer(self.nbeads.fetch(), initlist )

