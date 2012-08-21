import numpy as np
from engine.normalmodes import *
from utils.inputvalue import *
from utils.units import *
from copy import copy

__all__ = ['InputNormalModes']

class InputNormalModes(InputArray):
   """ Storage class for NormalModes engine.

   Describes how normal-modes transformation and integration should be
   performed.
   """

   attribs=copy(InputArray.attribs)
   attribs["mode"]= (InputValue, {"dtype"   : str,
                                    "default" : "rpmd",
                                    "help"    : "How the dynamical masses for the .",
                                    "options" : ['pa-cmd', 'wmax-cmd', 'manual', 'rpmd']})

   #~ fields = {
               #~ "frequencies" : (InputArray, {"dtype"        : float,
                                    #~ "default"        : np.identity(0),
                                    #~ "help"           : "Manual frequencies for the ring polymer normal modes. Just one number if doing CMD.",
                                    #~ "dimension"      : "frequency"})
             #~ }

   default_help = "Describes how the dynamic of the path must be performed, by manipulating the frequencies of normal modes."
   default_label = "NORMALMODES"

   def __init__(self, help=None, dimension=None, units=None, default=None, dtype=None):
      """ Initializes an InputProperties object by just calling the parent
          with appropriate arguments. """

      super(InputNormalModes,self).__init__(dtype=float, dimension="frequency", default=default, units=units, help=help)

   def store(self, nm):

      super(InputNormalModes,self).store(nm.nm_freqs)
      self.mode.store(nm.mode)

   def fetch(self):

      super(InputNormalModes,self).check()
      return NormalModes(self.mode.fetch(), super(InputNormalModes,self).fetch() )



