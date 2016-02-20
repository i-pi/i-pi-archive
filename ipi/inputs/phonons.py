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

Inputs created by Michele Ceriotti and Benjamin Helfrecht, 2015

Classes:
   InputGeop: Deals with creating the Geop object from a file, and
      writing the checkpoints.
"""

import numpy as np
import ipi.engine.initializer
from ipi.engine.motion import *
from ipi.utils.inputvalue import *
from ipi.inputs.thermostats import *
from ipi.inputs.initializer import *
from ipi.utils.units import *

__all__ = ['InputDynMatrix']

class InputDynMatrix(InputDictionary):
    """Dynamic matrix calculation options.
    
       Contains options related with finite difference computation of force constats. 

    """

    attribs={"mode"  : (InputAttribute, {"dtype"   : str, "default": "std",
                                    "help"    : "The algorithm to be used",
                                    "options" : ['std', 'nrg', 'ref']}) }
    fields = { 
                "pos_shift"  : (InputValue, {"dtype"   : float, "default": 0.01, 
                                    "help"    : "The finite deviation in position used to compute derivative of force."
                                    }), 
                "energy_shift"  : (InputValue, {"dtype"   : float, "default": 0.000, 
                                    "help"    : "The finite deviation in energy used to compute deribvative of force."
                                    }), 
                "output_shift"  : (InputValue, {"dtype"   : float, "default": 0.000, 
                                    "help"    : "Shift by the dynamical matrix diagonally before outputting."
                                    }),
                "prefix"  : (InputValue, {"dtype"   : str, "default": "PHONONS", 
                                    "help"    : "Shift by this much the dynamical matrix in the output."
                                    }),  
                "asr"  : (InputValue, {"dtype"   : str, "default": "none", "options" : ["none", "simple", "balanced", "crystal" ],
                                    "help"    : "Shift by this much the dynamical matrix in the output."
                                    }),   
                "dynmat" : ( InputArray, {"dtype" : float, 
                              "default" :  np.zeros(0, float),
                              "help"    : "Portion of the dynamical matrix known up to now."}),
                "dynmat_r" : ( InputArray, {"dtype" : float, 
                              "default" :  np.zeros(0, float),
                              "help"    : "Portion of the dynamical matrix known up to now (refining)."})              
             }
                   
    dynamic = {  }

    default_help = "Fill in."
    default_label = "PHONONS"

    def store(self, phonons):
        if phonons  == {}: return
        self.mode.store(phonons.mode)
        self.pos_shift.store(phonons.deltax)
        self.energy_shift.store(phonons.deltae)        
        self.output_shift.store(phonons.deltaw)     
        self.prefix.store(phonons.prefix)
        self.asr.store(phonons.asr)
        self.dynmat.store(phonons.dynmatrix)
        self.dynmat_r.store(phonons.dynmatrix_r)
        
    def fetch(self):		
        rv = super(InputDynMatrix,self).fetch()
        rv["mode"] = self.mode.fetch()        
        return rv
