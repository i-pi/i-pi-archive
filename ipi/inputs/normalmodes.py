"""Deals with creating the normal mode representation arrays.

Copyright (C) 2013, Joshua More and Michele Ceriotti

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


Classes:
   InputNormalModes: Deals with creating the normal mode objects.
"""

import numpy as np
from copy import copy
from ipi.engine.normalmodes import *
from ipi.utils.inputvalue import *
from ipi.utils.units import *

__all__ = ['InputNormalModes']

class InputNormalModes(InputArray):
   """ Storage class for NormalModes engine.

   Describes how normal-modes transformation and integration should be
   performed.

   Attributes:
      mode: Specifies the method by which the dynamical masses are created.
      transform: Specifies whether the normal mode calculation will be
         done using a FFT transform or a matrix multiplication.
   """

   attribs = copy(InputArray.attribs)
   attribs["mode"] = (InputAttribute, {"dtype"   : str,
                                       "default" : "rpmd",
                                       "help"    : "Specifies the technique to be used to calculate the dynamical masses. 'rpmd' simply assigns the bead masses the physical mass. 'manual' sets all the normal mode frequencies except the centroid normal mode manually. 'pa-cmd' takes an argument giving the frequency to set all the non-centroid normal modes to. 'wmax-cmd' is similar to 'pa-cmd', except instead of taking one argument it takes two ([wmax,wtarget]). The lowest-lying normal mode will be set to wtarget for a free particle, and all the normal modes will coincide at frequency wmax. ",
                                       "options" : ['pa-cmd', 'wmax-cmd', 'manual', 'rpmd']})
   attribs["transform"] = (InputValue,{"dtype"   : str,
                                       "default" : "fft",
                                       "help"    : "Specifies whether to calculate the normal mode transform using a fast Fourier transform or a matrix multiplication. For small numbers of beads the matrix multiplication may be faster.",
                                       "options" : ['fft', 'matrix']})

   default_help = "Deals with the normal mode transformations, including the adjustment of bead masses to give the desired ring polymer normal mode frequencies if appropriate. Takes as arguments frequencies, of which different numbers must be specified and which are used to scale the normal mode frequencies in different ways depending on which 'mode' is specified."
   default_label = "NORMALMODES"

   def __init__(self, help=None, dimension=None, default=None, dtype=None):
      """ Initializes InputNormalModes.

      Just calls the parent initialization function with appropriate arguments.
      """

      super(InputNormalModes,self).__init__(help=help, default=default, dtype=float, dimension="frequency")

   def store(self, nm):
      """Takes a normal modes instance and stores a minimal representation
      of it.

      Args:
         nm: A normal modes object.
      """

      super(InputNormalModes,self).store(nm.nm_freqs)
      self.mode.store(nm.mode)
      self.transform.store(nm.transform_method)

   def fetch(self):
      """Creates a normal modes object.

      Returns:
         A normal modes object.
      """

      super(InputNormalModes,self).check()
      return NormalModes(self.mode.fetch(), self.transform.fetch(), super(InputNormalModes,self).fetch() )
