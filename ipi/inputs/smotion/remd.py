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
from ipi.engine.motion import *
from ipi.utils.inputvalue import *
from ipi.utils.units import *

__all__ = ['InputReplicaExchange']

class InputReplicaExchange(InputDictionary):
    """Geometry optimization options.

    Contains options related with replica exchange, such as method,
    steps on which REMD should be performed, etc.

    """

    fields={
           "stride" : (InputValue, {"dtype"        : float,
                                      "default"      : 1.0,
                                      "help"         : "Every how often to try exchanges (on average)."
                                      }),
         }

    default_help = "Replica Exchange"
    default_label = "REMD"

    def store(self, remd):
        if remd == {}: return
        self.stride.store(remd.stride)

    def fetch(self):
        rv = super(InputReplicaExchange,self).fetch()
        rv["stride"] = self.stride.fetch()
        return rv
