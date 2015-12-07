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
from ipi.engine.mover import *
from ipi.utils.inputvalue import *
from ipi.inputs.thermostats import *
from ipi.inputs.initializer import *
from ipi.utils.units import *

__all__ = ['InputForceConst']

class InputForceConst(InputDictionary):
    """Dynamic matrix calculation options.
    
       Contains options related with finite difference computation of force constats. 

    """

    #TODO: RENAME BFGS/L-BFGS ONLY OPTIONS TO INDICATE SPECIFIC TO BFGS/L-BFGS
    attribs={"epsilon"  : (InputAttribute, {"dtype"   : float, "default": 0.01, 
                                    "help"    : "The finite deviation used to compute deribvative of force."
                                    }) }
   
    fields = { 
                "oldk" : ( InputDictionary, {"dtype" : int, 
                              "default" : 0.01,
                              "help"    : "Number of rows of the Hessian calculated until previous step."}),
                "oldhessian" : ( InputArray, {"dtype" : float, 
                              "default" :  np.zeros(0, float),
                              "help"    : "Hessian known until previous step."})
             }
                   
    dynamic = {  }

    default_help = "Fill in."
    default_label = "PHONONS"

    def store(self, phonons):
        print "inside phonons"
        if phonons  == {}: return
        self.epsilon.store(phonons.epsilon)
        self.oldk.store(phonons.oldk)
        self.oldhessian.store(phonons.oldhessian)
        
    def fetch(self):		
        rv = super(phonons,self).fetch()
        rv.epsilon = self.epsilon.fetch()        
        rv.oldk = self.oldk.fetch()        
        rv.oldhessian = self.oldhessian.fetch()        
        return rv
