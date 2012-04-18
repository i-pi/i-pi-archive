"""Contains the classes that connect the driver to the python code.

Communicates with the driver code, obtaining the force, virial and potential.
Deals with creating the jobs that will be sent to the driver, and 
returning the results to the python code.

Classes:
   RestartForce: Deals with creating the ForceField object from a file, and
      writing the checkpoints.
   ForceField: Base forcefield class with the generic methods and attributes.
   FFSocket: Deals with a single replica of the system
   ForceBeads: Deals with the parallelization of the force calculation for
      a PI simulation.
"""

__all__ = ['RestartForce']

import numpy as np
import math, time
from utils.depend import *
from engine.forces import *
from inputs.interface import RestartInterface
from utils.restart import *
from utils.inputvalue import *

class RestartForce(Input):
   """Forcefield restart class.

   Handles generating the appropriate forcefield class from the xml
   input file, and generating the xml checkpoint tags and data from an 
   instance of the object.

   Attributes:
      type: A string indicating the type being used. 'socket' is currently
         the only available option.
      interface: A restart interface instance.
      parameters: A dictionary of the parameters used by the driver. Of the
         form {"name": value}.
   """

   attribs = { "type" : ( InputValue, { "dtype"   :  str, 
                                       "default" : "socket",
                                       "options" : [ "socket" ],
                                       "help"    : "Well, we'll write a help string some day"  }  )}
   fields =  { "interface"  : ( RestartInterface, { } ),                      # THESE WILL HAVE TO BE RE-WRITTEN WHEN THE CORRESPONDING INTERFACES ARE UPDATED
               "parameters" : ( RestartValue, { "dtype" : dict, 
                                                "default" : {} }) }
   
   def store(self, force):
      """Takes a ForceField instance and stores a minimal representation of it.

      Args:
         force: A forcefield object.
      """

      if (type(force) is FFSocket):  
         self.type.store("socket")
         self.interface.store(force.socket)
         self.parameters.store(force.pars)
      else: 
         self.type.store("unknown")
         

   def fetch(self):
      """Creates a forcefield object.

      Returns:
         A forcefield object of the appropriate type and with the appropriate
         interface given the attributes of the RestartForce object.
      """

      if self.type.fetch().upper() == "SOCKET": 
         force = FFSocket(pars=self.parameters.fetch(), interface=self.interface.fetch())
      else: 
         force = ForceField()
      return force
