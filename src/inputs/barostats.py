"""Deals with creating the barostat class.

Classes:
   RestartBaro: Deals with creating the Barostat object from a file, and 
      writing the checkpoints.
"""

from engine.barostats import *
import engine.thermostats
from utils.inputvalue import *
from inputs.thermostats import *

__all__ = ['RestartBaro']
      
class RestartBaro(Input):
   """Barostat input class.

   Handles generating the appropriate barostat class from the xml input file, 
   and generating the xml checkpoint tags and data from an 
   instance of the object.

   Attributes:
      kind: An optional string giving the type of barostat used. Defaults to
         'rigid'.
      thermostat: A thermostat object giving the cell thermostat.
   """

   attribs={ "kind": (InputValue, {"dtype"    : str, 
                                   "default"  : "rigid",
                                   "help"     : "The type of barostat.",
                                   "options"  : ["rigid", "flexible"]}) }
   fields={ "thermostat": (RestartThermo, {"default" : engine.thermostats.Thermostat(),
                                           "help"    : "The thermostat for the cell. Keeps the cell velocity distribution at the correct temperature."}) }

   def store(self, baro):
      """Takes a barostat instance and stores a minimal representation of it.

      Args:
         baro: A barostat object.
      """

      super(RestartBaro,self).store(baro)
      if type(baro) is BaroRigid:
         self.kind.store("rigid")
      if type(baro) is BaroFlexi:
         self.kind.store("flexible")
      else:
         self.kind.store("unknown")      
      self.thermostat.store(baro.thermostat)
      
   def fetch(self):
      """Creates a barostat object.

      Returns:
         A barostat object of the appropriate type and with the appropriate 
         thermostat given the attributes of the RestartBaro object.
      """

      super(RestartBaro,self).fetch()
      if self.kind.fetch().upper() == "RIGID":
         baro=BaroRigid(thermostat=self.thermostat.fetch())
      elif self.kind.fetch().upper() == "FLEXIBLE":
         baro=BaroFlexi(thermostat=self.thermostat.fetch())
      else:
         baro=Barostat(thermostat=self.thermostat.fetch())

      return baro
