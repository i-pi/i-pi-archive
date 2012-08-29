"""Deals with creating the barostat class.

Classes:
   InputBaro: Deals with creating the Barostat object from a file, and
      writing the checkpoints.
"""

from engine.barostats import *
import engine.thermostats
from utils.inputvalue import *
from inputs.thermostats import *

__all__ = ['InputBaro']

class InputBaro(Input):
   """Barostat input class.

   Handles generating the appropriate barostat class from the xml input file,
   and generating the xml checkpoint tags and data from an
   instance of the object.

   Attributes:
      kind: An optional string giving the type of barostat used. Defaults to
         'rigid'.
      thermostat: A thermostat object giving the cell thermostat.
   """

   attribs={ "kind": (InputAttribute, {"dtype"    : str,
                                   "default"  : "rigid",
                                   "help"     : "The type of barostat. 'Rigid' gives a barostat that keeps the internal pressure constant by allowing cell volume changes, whereas flexible allows the shape of the cell to fluctuate too.",
                                   "options"  : ["rigid"]}) }
   fields={ "thermostat": (InputThermo, {"default" : input_default(factory=engine.thermostats.Thermostat),
                                         "help"    : "The thermostat for the cell. Keeps the cell velocity distribution at the correct temperature."}) }

   default_help = "Simulates an external pressure bath to keep the pressure or stress at the external values."
   default_label = "BAROSTAT"

   def store(self, baro):
      """Takes a barostat instance and stores a minimal representation of it.

      Args:
         baro: A barostat object.
      """

      super(InputBaro,self).store(baro)
      if type(baro) is BaroRigid or type(baro) is Barostat:
         self.kind.store("rigid")
      else:
         raise TypeError("The type " + type(baro).__name__ + " is not a valid barostat type")
      self.thermostat.store(baro.thermostat)

   def fetch(self):
      """Creates a barostat object.

      Returns:
         A barostat object of the appropriate type and with the appropriate
         thermostat given the attributes of the InputBaro object.
      """

      super(InputBaro,self).fetch()
      if self.kind.fetch() == "rigid":
         baro=BaroRigid(thermostat=self.thermostat.fetch())
      else:
         raise ValueError("Kind " + self.kind.fetch() + " is not a valid kind of barostat")

      return baro
