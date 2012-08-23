from utils.inputvalue import *
from copy import copy
from inputs.beads import InputBeads
from inputs.cell import InputCell
import engine.initializer as ei

__all__ = ['InputInitializer', 'InputInitFile']


class InputInitFile(InputValue):
   """Class to handle input from file. """

   attribs=copy(InputValue.attribs)
   attribs["format"]=(InputValue,{ "dtype" : str, "default": "xyz"} )

   def __init__(self, help=None, dimension=None, units=None, default=None, dtype=None):
      """ Initializes an InputInitFile object by just calling the parent
          with appropriate arguments. """

      super(InputInitFile,self).__init__(dtype=str, dimension=dimension, default=default, units=units, help=help)

   def store(self,iif):
      super(InputInitFile,self).store(iif.filename)
      self.format.store(iif.format)

   def fetch(self):
      return ei.InitFile(filename=super(InputInitFile,self).fetch(), format=self.format.fetch() )

class InputInitializer(Input):
   """Input class to handle initialization.
   """

   attribs = { "nbeads"    : (InputValue, {"dtype"     : int,
                                        "help"      : "The number of beads. Will override any provision from inside the initializer."})
            }

   dynamic = {
           "beads" : (InputBeads, { "help" : "Inputs the configuration of the path as a Beads object" }),
           "cell" : (InputCell, { "help" : "Inputs the configuration of the cell as a Cell object" }),
           "file" : (InputInitFile, {"help" : "Reads bead(s) and cell configuration from an external file" })
           }

   default_help = "Specifies the number of beads, and how the system should be initialized."
   default_label = "INITIALIZER"

   def write(self,  name="", indent=""):
      """Overloads Input write() function so that we never write out InputInitializer to restart files.

      Returns:
         An empty string.
      """

      return ""

   def store(self, ii):

      print "Storing input thing", ii
      self.extra = []

      for (k, el) in ii.queue:
         if k == "file" :
            ip = InputInitFile()
            ip.store(el)
            self.extra.append(("file", ip))
         elif k == "beads" :
            ip = InputBeads()
            ip.store(el)
            self.extra.append(("beads", ip))
         elif k == "cell" :
            ip = InputCell()
            ip.store(el)
            self.extra.append(("cell", ip))

      self.nbeads.store(ii.nbeads)



   def fetch(self):

      super(InputInitializer,self).fetch()

      return ei.Initializer(self.nbeads.fetch(), [ (k,v.fetch())  for (k,v) in self.extra ] )

