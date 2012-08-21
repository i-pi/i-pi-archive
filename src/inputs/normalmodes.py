import numpy as np
from engine.normalmodes import *
from utils.inputvalue import *
from utils.units import *

__all__ = ['InputNormalModes']

class InputNormalModes(Input):
   """ Storage class for NormalModes engine.

   Describes how normal-modes transformation and integration should be
   performed.
   """

   attribs = {
               "mode"  : (InputValue, {"dtype"   : str,
                                    "default" : "rpmd",
                                    "help"    : "How the dynamical masses for the .",
                                    "options" : ['pa-cmd', 'wmax-cmd', 'manual', 'rpmd']})
             }

   fields = {
               "frequencies" : (InputArray, {"dtype"        : float,
                                    "default"        : ClassDefault(type=np.identity, args=(0,)),
                                    "help"           : "Manual frequencies for the ring polymer normal modes. Just one number if doing CMD.",
                                    "dimension"      : "frequency"})
             }

   default_help = "Describes how the dynamic of the path must be performed, by manipulating the frequencies of normal modes."
   default_label = "NORMALMODES"

   def store(self, nm):

      self.mode.store(nm.mode)
      self.frequencies.store(nm.nm_freqs)

   def fetch(self):

      super(InputNormalModes,self).check()
      return NormalModes(self.mode.fetch(), self.frequencies.fetch())



