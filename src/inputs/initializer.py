"""Deals with creating the initiliazer class.

Classes:
   InputInitializer: Initializes the classes that initialize the simulation
      data.
   InputInitFile: Initializes the classes that initialize the simulation data
      from a file.
"""

from utils.inputvalue import *
from copy import copy
from inputs.beads import InputBeads
from inputs.cell import InputCell
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

