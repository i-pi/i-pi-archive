"""Deals with creating the barostat class.

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
                                   "help"     : "The type of barostat. 'bzp' gives a Bussi-Zykova-Parrinello isotropic barostat. 'mht' gives a Martyna-Hughes-Tuckerman isotropic barostat.",
                                   "options"  : ["dummy", "bzp", "mht"]}) }
   fields={ "thermostat": (InputThermo, {"default" : input_default(factory=ipi.engine.thermostats.Thermostat),
                                         "help"    : "The thermostat for the cell. Keeps the cell velocity distribution at the correct temperature."}),
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
         self.mode.store("bzp")
         self.p.store(baro.p)
      elif type(baro) is BaroMHT:
         self.mode.store("mht")
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
      if self.mode.fetch() == "bzp":
         baro = BaroBZP(thermostat=self.thermostat.fetch(), tau=self.tau.fetch())
         if self.p._explicit: baro.p = self.p.fetch()
      elif self.mode.fetch() == "mht":
         baro = BaroMHT(thermostat=self.thermostat.fetch(), tau=self.tau.fetch())
         if self.p._explicit: baro.p = self.p.fetch()
      elif self.mode.fetch() == "dummy":
         baro = Barostat(thermostat=self.thermostat.fetch(), tau=self.tau.fetch())
      else:
         raise ValueError(self.mode.fetch() + " is not a valid mode of barostat")

      return baro
