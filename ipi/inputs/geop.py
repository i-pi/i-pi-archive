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
from ipi.engine.mover import *
from ipi.utils.inputvalue import *
from ipi.inputs.thermostats import *
from ipi.inputs.initializer import *
from ipi.utils.units import *

__all__ = ['InputGeop']

class InputGeop(InputDictionary):
    """Geometry optimization options.
    
    Contains options related with geometry optimization, such as method, 
    thresholds, linear search strategy, etc. 

    """

    attribs={"mode"  : (InputAttribute, {"dtype"   : str, "default": "cg", 
                                    "help"    : "The geometry optimization algorithm to be used",
                                    "options" : ['sd', 'cg', 'bfgs']}) }
   
    fields = { "ls_options" : ( InputDictionary, {"dtype" : float, 
                              "help" : """Options for line search methods. Includes: 
                              tolerance: stopping tolerance (for position) for the search,
                              grad_tolerance: stopping tolerance on gradient for 
                              BFGS line search,
                              iter: the maximum number of iterations,
                              step: initial step for bracketing,
                              adaptive: whether to update line_step.
                              """, 
                              "options" : ["tolerance", "gradtolerance",  "iter", "step", "adaptive"],
                              "default" : [1e-6, 1e-6, 100, 1e-3, True] }),       
                "tolerances" : ( InputDictionary, {"dtype" : float, 
                              "help"    : """Requirements for converged relaxation""",
                              "options" : [ "energy", "force", "position" ],
                              "default" : [ 1e-8, 1e-8, 1e-8 ] }),
                "cg_old_force": (InputArray, {"dtype" : float,
                              "default"   : input_default(factory=np.zeros, args = (0,)),
                              "help"      : "The previous force in a CG optimization.",
                              "dimension" : "force"}),
                "cg_old_direction": (InputArray, {"dtype" : float,
                              "default" : input_default(factory=np.zeros, args = (0,)),
                              "help"    : "The previous direction in a CG optimization."}),
                "maximum_step": (InputValue, {"dtype" : float,
                              "default" : 100.0,
                              "help"    : "The maximum step size for BFGS line minimizations."}),
                "invhessian" : (InputArray, {"dtype" : float, 
                              "default" : input_default(factory=np.eye, args = (0,)),
                              "help"    : "Approximate inverse Hessian for BFGS, if known."})
                     }
                   
    dynamic = {  }

    default_help = "TODO EXPLAIN WHAT THIS IS"
    default_label = "GEOP"   

    def store(self, geop):
        if geop == {}: return
        self.ls_options.store(geop.ls_options)
        self.tolerances.store(geop.tolerances)
        self.mode.store(geop.mode)
        self.cg_old_force.store(geop.cg_old_f)
        self.cg_old_direction.store(geop.cg_old_d)
        self.maximum_step.store(geop.max_step)
        
    def fetch(self):		
        rv = super(InputGeop,self).fetch()
        rv["mode"] = self.mode.fetch()        
        return rv
