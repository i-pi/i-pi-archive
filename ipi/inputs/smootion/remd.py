"""Deals with creating the ensembles class.

Copyright (C) 2013, Joshua More and Michele Ceriotti

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
   InputEnsemble: Deals with creating the Ensemble object from a file, and
      writing the checkpoints.
"""

import numpy as np
import ipi.engine.initializer
from ipi.engine.motion import *
from ipi.utils.inputvalue import *
from ipi.inputs.thermostats import *
from ipi.inputs.initializer import *
from ipi.utils.units import *

__all__ = ['InputREMD']

class InputREMD(Input):
    """Geometry optimization options.

    Contains options related with replica exchange, such as method,
    steps on which REMD should be performed, etc.

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

    default_help = "Replica Exchange"
    default_label = "REMD"

    def store(self, remd):
        if remd == {}: return
        self.temp_list.store(remd.temp_list)
        self.temp_index.store(remd.temp_index)
        self.stride.store(remd.stride)

    def fetch(self):
        rv = super(InputREMD,self).fetch()
        rv["mode"] = self.mode.fetch()
        return rv
