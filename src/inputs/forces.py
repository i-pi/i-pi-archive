"""Deals with creating the forcefield class.

Classes:
   InputForce: Deals with creating the ForceField object from a file, and
      writing the checkpoints.
"""

__all__ = ['InputForce']

from engine.forces import *
from inputs.interface import InputInterface
from utils.inputvalue import *

class InputForce(Input):
   """Forcefield input class.

   Handles generating the appropriate forcefield class from the xml
   input file, and generating the xml checkpoint tags and data from an 
   instance of the object.

   Attributes:
      type: A string indicating the type being used. 'socket' is currently
         the only available option.
      interface: A restart interface instance.
      parameters: A dictionary of the parameters used by the driver. Of the
         form {"name": value}.
      nreduced: An integer giving the number of beads to use if a ring polymer
         contraction scheme is being used.
   """

   attribs = { "type"     : ( InputValue, { "dtype"   : str, 
                                            "default" : "socket",
                                            "options" : ["socket"],
                                            "help"    : "Specifies which kind of force object is created."  }  ), 
               "nreduced" : ( InputValue, { "dtype"   : int,
                                            "default" : 0,
                                            "help"    : "If the forcefield is to be evaluated on a contracted ring polymer, this gives the number of beads that are used. If not specified, the forcefield will be evaluated on the full ring polymer." } ) }
   fields =  { "interface"  : ( InputInterface, {"help" : InputInterface.default_help } ),
               "parameters" : ( InputValue, { "dtype"   : dict, 
                                              "default" : input_default(factory=dict),
                                              "help"    : "Deprecated dictionary of initialization parameters. May be removed in the future." }) }

   default_help = "Deals with the assigning of jobs to different driver codes, and collecting the data."
   default_label = "FORCES"
   
   def store(self, force):
      """Takes a ForceField instance and stores a minimal representation of it.

      Args:
         force: A forcefield object.
      """

      super(InputForce,self).store(force)
      if (type(force) is FFSocket):  
         self.type.store("socket")
         self.interface.store(force.socket)
         self.parameters.store(force.pars)
         self.nreduced.store(force.nbeads)
      else: 
         raise TypeError("The type " + type(force).__name__ + " is not a valid forcefield type")

   def fetch(self):
      """Creates a forcefield object.

      Returns:
         A forcefield object of the appropriate type and with the appropriate
         interface given the attributes of the InputForce object.
      """

      super(InputForce,self).fetch()
      if self.type.fetch() == "socket": 
         force = FFSocket(pars=self.parameters.fetch(), 
            interface=self.interface.fetch())
         if self.nreduced._explicit:
            force.nbeads = self.nreduced.fetch()
      else: 
         raise ValueError("Kind " + self.kind.fetch() + " is not a valid kind of forcefield")

      return force
