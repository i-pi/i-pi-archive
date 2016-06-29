"""Deals with creating the ParaTemp class."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

from ipi.utils.inputvalue import *
from ipi.engine.paratemp import ParaTemp


__all__=[ "InputParaTemp" ]


class InputParaTemp(Input):
   """Input class for the ParaTemp object.

   Contains all the options for a parallel tempering simulation.

   Fields:
      temp_list: The list of temperatures for the replicas. Must match
          the size of the system list.
      stride: How often -- on average -- an exchange should be made.
   """


   fields={"temp_list" : (InputArray, {"dtype": float,
                                       "default"   : input_default(factory=np.zeros, args = (0,)),
                                         "help"      : "List of temperatures for a parallel tempering simulation",
                                         "dimension" : "temperature" }),
           "temp_index" : (InputArray, {"dtype": int,
                                       "default"   : input_default(factory=np.zeros, args = (0,int)),
                                         "help"      : "Maps the temperatures to the list of systems."
                                         }),
           "stride" : (InputValue, {"dtype"        : float,
                                      "default"      : 0.0,
                                      "help"         : "Every how often to try exchanges (on average)."
                                      }),
         }

   default_help = "Contains all the options for a parallel tempering simulation."
   default_label = "PARATEMP"

   def __init__(self, help=None,  default=None):
      """Initializes InputParaTemp.

      Just calls the parent initialization function with appropriate arguments.
      """

      super(InputParaTemp,self).__init__(help=help, default=default)

   def store(self, pt):
      """Stores a ParaTemp object."""

      self.temp_list.store(pt.temp_list)
      self.temp_index.store(pt.temp_index)
      self.stride.store(pt.stride)


   def fetch(self):
      """Creates a ParaTemp object based on the input parameters."""

      return ParaTemp(self.temp_list.fetch(), self.temp_index.fetch(), self.stride.fetch())
