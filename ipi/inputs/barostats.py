"""Deals with creating the barostat class.

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
   InputBaro: Deals with creating the Barostat object from a file, and
      writing the checkpoints.
"""

import numpy as np
import ipi.engine.thermostats
from ipi.engine.barostats import *
from ipi.utils.inputvalue import *
from ipi.inputs.thermostats import *

__all__ = ['InputBaro']

class InputBaro(Input):
   """Barostat input class.

   Handles generating the appropriate barostat class from the xml input file,
   and generating the xml checkpoint tags and data from an
   instance of the object.

   Attributes:
      mode: An optional string giving the type of barostat used. Defaults to
         'rigid'.

   Fields:
      thermostat: A thermostat object giving the cell thermostat.
      tau: The time constant associated with the dynamics of the piston.
      p: The conjugate momentum to the volume degree of freedom.
   """

   attribs={ "mode": (InputAttribute, {"dtype"    : str,
                                   "default" : "dummy",
                                   "help"     : """The type of barostat.  Currently, only a 'isotropic' barostat is implemented, that combines
                                   ideas from the Bussi-Zykova-Parrinello barostat for classical MD with ideas from the
                                   Martyna-Hughes-Tuckerman centroid barostat for PIMD; see Ceriotti, More, Manolopoulos, Comp. Phys. Comm. 2013 for
                                   implementation details.""",
                                   "options"  : ["dummy", "isotropic"]}) }
   fields={ "thermostat": (InputThermo, {"default" : input_default(factory=ipi.engine.thermostats.Thermostat),
                                         "help"    : "The thermostat for the cell. Keeps the cell velocity distribution at the correct temperature. Note that the 'pile_l', 'pile_g', 'nm_gle' and 'nm_gle_g' options will not work for this thermostat."}),
            "tau": (InputValue, {"default" : 1.0,
                                  "dtype" : float,
                                  "dimension" : "time",
                                  "help"    : "The time constant associated with the dynamics of the piston."}),
            "p": (InputArray, {  "dtype"     : float,
                                 "default"   : input_default(factory=np.zeros, args = (0,)),
                                 "help"      : "Momentum (or momenta) of the piston.",
                                 "dimension" : "momentum" })
           }

   default_help = "Simulates an external pressure bath."
   default_label = "BAROSTAT"

   def store(self, baro):
      """Takes a barostat instance and stores a minimal representation of it.

      Args:
         baro: A barostat object.
      """

      super(InputBaro,self).store(baro)
      self.thermostat.store(baro.thermostat)
      self.tau.store(baro.tau)
      if type(baro) is BaroBZP:
         self.mode.store("isotropic")
         self.p.store(baro.p)
      elif type(baro) is Barostat:
         self.mode.store("dummy")
      else:
         raise TypeError("The type " + type(baro).__name__ + " is not a valid barostat type")


   def fetch(self):
      """Creates a barostat object.

      Returns:
         A barostat object of the appropriate type and with the appropriate
         thermostat given the attributes of the InputBaro object.
      """

      super(InputBaro,self).fetch()
      if self.mode.fetch() == "isotropic":
         baro = BaroBZP(thermostat=self.thermostat.fetch(), tau=self.tau.fetch())
         if self.p._explicit: baro.p = self.p.fetch()
      elif self.mode.fetch() == "dummy":
         baro = Barostat(thermostat=self.thermostat.fetch(), tau=self.tau.fetch())
      else:
         raise ValueError(self.mode.fetch() + " is not a valid mode of barostat")

      return baro
