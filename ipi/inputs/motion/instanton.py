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
   InputInst: Deals with creating the Inst object from a file, and
      writing the checkpoints.
"""

import numpy as np
from ipi.utils.inputvalue import *


__all__ = ['InputInst']

class InputInst(InputDictionary):
    """Instanton optimization options.

    Contains options related with instanton, such as method,
    thresholds, hessian update strategy, etc.

    """
    attribs={"mode"  : (InputAttribute, {"dtype"   : str, "default": "full",
                                    "help"    : "Use the full or half of the ring polymer during the optimization",
                                    "options" : ['full','half']}) }

    fields = { "tolerances" : ( InputDictionary, {"dtype" : float,
                              "options"  : [ "energy", "force", "position" ],
                              "default"  : [ 1e-6, 1e-6, 1e-6 ],
                              "help"     : "Convergence criteria for optimization. Default values are extremely conservative. Set them to appropriate values for production runs.",
                              "dimension": [ "energy", "force", "length" ] }),
               "biggest_step": (InputValue, {"dtype" : float,
                              "default"  : 100.0,
                              "help"     : "The maximum step size for (L)-BFGS line minimizations."}),
               "old_pos":    (InputArray, {"dtype" : float,
                              "default"  : input_default(factory=np.zeros, args = (0,)),
                              "help"     : "The previous positions in an optimization step.",
                              "dimension": "length"}),
               "old_pot":    (InputArray, {"dtype" : float,
                              "default"  : input_default(factory=np.zeros, args = (0,)),
                              "help"     : "The previous potential energy in an optimization step.",
                              "dimension": "energy"}),
               "old_force":  (InputArray, {"dtype" : float,
                              "default"  : input_default(factory=np.zeros, args = (0,)),
                              "help"     : "The previous force in an optimization step.",
                              "dimension": "force"}),
               "hessian":    (InputArray, {"dtype" : float,
                              "default"  : input_default(factory=np.eye, args = (0,)),
                              "help"     : "Approximate Hessian for trm, if known."}),
               "delta":       (InputValue, {"dtype": float,
                               "default": 0.1,
                               "help": "Initial stretch amplitude."}),
               "hessian_init": (InputValue, {"dtype": str,
                                "default": "expand",
                                "options": ["true","expand"],
                                "help": "How to initialize the hessian if it is not fully provided. If true computes the Hessian from scratch and if expand, reads from a Hessian that was provided and interpolates if necessary."}),
               "hessian_update": (InputValue, {"dtype": str,
                                "default": "powell",
                                "options": ["powell", "recompute"],
                                "help": "How to update the hessian in each step. If powell uses the Powell method and if recompute, recomputes the Hessian from scratch."}),
               "hessian_asr": (InputValue, {"dtype": str,
                                "default": "none",
                                "options": ["none","poly","crystal"],
                                "help": "Removes the zero frequency vibrational modes depending on the symmerty of the system."}),
               "final_rates": (InputValue, {"dtype": str,
                                               "default": "false",
                                               "options": ["false", "true"],
                                               "help": "How to update the hessian in each step."})

               } 

    dynamic = {  }

    default_help = "TODO EXPLAIN WHAT THIS IS"
    default_label = "Instanton"

    def store(self, geop):
        if geop == {}: 
            return


        self.mode.store(geop.mode)
        self.tolerances.store(geop.tolerances)
        self.hessian.store(geop.hessian)
        self.biggest_step.store(geop.big_step)
        self.delta.store(geop.delta)
        self.hessian_init.store(geop.hessian_init)
        self.hessian_update.store(geop.hessian_update)
        self.hessian_asr.store(geop.hessian_asr)
        self.final_rates.store(geop.final_rates)



    def fetch(self):
        rv = super(InputInst,self).fetch()
        rv["mode"] = self.mode.fetch()
        return rv
