"""Deals with creating the thermostats class.

Chooses between the different possible thermostat options and creates the
appropriate thermostat object, with suitable parameters.

Classes:
   InputThermo: Deals with creating the thermostat object from a file, and
      writing the checkpoints.
"""

__all__ = ['InputThermo']

import numpy as np
import math
from utils.depend   import *
from utils.inputvalue  import *
from engine.thermostats import *
            
class InputThermo(Input):
   """Thermostat input class.

   Handles generating the appropriate thermostat class from the xml input file,
   and generating the xml checkpoiunt tags and data from an instance of the
   object.

   Attributes:
      kind: An optional string giving the type of the thermostat used. Defaults
         to 'langevin'.
      ethermo: An optional float giving the amount of heat energy transferred
         to the bath. Defaults to 0.0.
      tau: An optional float giving the damping time scale. Defaults to 1.0.
      A: An optional array of floats giving the drift matrix. Defaults to 0.0.
      C: An optional array of floats giving the static covariance matrix. 
         Defaults to 0.0.
   """

   attribs = { "kind": (InputValue, { "dtype"   : str, 
                                      "default" : "langevin",
                                      "options" : [ "langevin", "svr", "pile_l", "pile_g", "gle", "nm_gle" ],
                                      "help"    : "The style of thermostatting."
                                         }) }
   fields = { "ethermo" : (InputValue, {  "dtype"     : float, 
                                          "default"   : 0.0,
                                          "help"      : "The initial value of the thermostat energy. Only useful in restarts to guarantee continuity of the conserved quantity. ",
                                          "dimension" : "energy" }), 
            "tau" : (InputValue, {  "dtype"     : float, 
                                    "default"   : 0.0,
                                    "help"      : "The target temperature for the thermostat.",
                                    "dimension" : "temperature" }) ,
            "A" : (InputArray, {    "dtype"     : float, 
                                    "default"   : np.zeros(0),
                                    "help"      : "The friction matrix for GLE thermostats.",
                                    "dimension" : "frequency" }),
            "C" : (InputArray, {    "dtype"     : float, 
                                    "default"   : np.zeros(0),
                                    "help"      : "The covariance matrix for GLE thermostats.",
                                    "dimension" : "temperature" }),
            "s" : (InputArray, {    "dtype"     : float, 
                                    "default"   : np.zeros(0),
                                    "help"      : "Input values for the additional momenta in GLE.",
                                    "dimension" : "ms-momentum" })
             }
   
   def store(self, thermo):
      """Takes a thermostat instance and stores a minimal representation of it.

      Args:
         thermo: A thermostat object.

      Raises:
         TypeError: Raised if the thermostat is not a recognized type.
      """

      super(InputThermo,self).store(thermo)
      if type(thermo) is ThermoLangevin: 
         self.kind.store("langevin")
         self.tau.store(thermo.tau)
      elif type(thermo) is ThermoSVR: 
         self.kind.store("svr")
         self.tau.store(thermo.tau)         
      elif type(thermo) is ThermoPILE_L: 
         self.kind.store("pile_l")
         self.tau.store(thermo.tau)
      elif type(thermo) is ThermoPILE_G: 
         self.kind.store("pile_g")
         self.tau.store(thermo.tau)     
      elif type(thermo) is ThermoGLE: 
         self.kind.store("gle")
         self.A.store(thermo.A)
         if dget(thermo,"C")._func is None:
            self.C.store(thermo.C)
         self.s.store(thermo.s)
      elif type(thermo) is ThermoNMGLE: 
         self.kind.store("nm_gle")
         self.A.store(thermo.A)
         if dget(thermo,"C")._func is None:
            self.C.store(thermo.C)
         self.s.store(thermo.s)
      else:
         raise TypeError("Unknown thermostat kind " + type(thermo).__name__)
      self.ethermo.store(thermo.ethermo)
      
   def fetch(self):
      """Creates a thermostat object.

      Returns:
         A thermostat object of the appropriate type and with the appropriate
         parameters given the attributes of the InputThermo object.

      Raises:
         TypeError: Raised if the thermostat type is not a recognized option.
      """

      super(InputThermo,self).fetch()
      if self.kind.fetch() == "langevin":
         thermo = ThermoLangevin(tau=self.tau.fetch())
      elif self.kind.fetch() == "svr":
         thermo = ThermoSVR(tau=self.tau.fetch())
      elif self.kind.fetch() == "pile_l":
         thermo = ThermoPILE_L(tau=self.tau.fetch())
      elif self.kind.fetch() == "pile_g":
         thermo = ThermoPILE_G(tau=self.tau.fetch())
      elif self.kind.fetch() == "gle":
         rC = self.C.fetch()
         if len(rC) == 0:
            rC = None
         thermo = ThermoGLE(A=self.A.fetch(),C=rC)
         thermo.s = self.s.fetch()
      elif self.kind.fetch() == "nm_gle":
         rC = self.C.fetch()
         if len(rC) == 0:
            rC = None
         thermo = ThermoNMGLE(A=self.A.fetch(),C=rC)
         thermo.s = self.s.fetch()
      else:
         raise TypeError("Invalid thermostat kind " + self.kind.fetch())
         
      thermo.ethermo = self.ethermo.fetch()

      return thermo
