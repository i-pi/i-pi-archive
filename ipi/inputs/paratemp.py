"""Deals with creating the ParaTemp class.

Copyright (C) 2013, Michele Ceriotti

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http.//www.gnu.org/licenses/>.

Classes:
   InputParaTemp: Defines all the input parameters for a ParaTemp object.
"""

__all__=[ "InputParaTemp" ]

import numpy as np
from ipi.utils.inputvalue import *
from ipi.engine.paratemp import ParaTemp

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


