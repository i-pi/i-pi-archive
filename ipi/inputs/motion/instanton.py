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
                              "default"  : [ 1e-5, 1e-5, 5e-3 ],
                              "help"     : "Convergence criteria for optimization.",
                              "dimension": [ "energy", "force", "length" ] }),
               "biggest_step": (InputValue, {"dtype" : float,
                              "default"  : 0.3,
                              "help"     : "The maximum step size during the optimization."}),
               "old_pos":    (InputArray, {"dtype" : float,
                              "default"  : input_default(factory=np.zeros, args = (0,)),
                              "help"     : "The previous step positions during the optimization. ",
                              "dimension": "length"}),
               "old_pot":    (InputArray, {"dtype" : float,
                              "default"  : input_default(factory=np.zeros, args = (0,)),
                              "help"     : "The previous step potential energy during the optimization",
                              "dimension": "energy"}),
               "old_force":  (InputArray, {"dtype" : float,
                              "default"  : input_default(factory=np.zeros, args = (0,)),
                              "help"     : "The previous step force during the optimization",
                              "dimension": "force"}),
               "hessian":    (InputArray, {"dtype" : float,
                              "default"  : input_default(factory=np.eye, args = (0,)),
                              "help"     : "(Approximate) Hessian."}),
               "delta":       (InputValue, {"dtype": float,
                               "default": 0.1,
                               "help": "Initial stretch amplitude."}),
               "hessian_init": (InputValue, {"dtype": str,
                                "default": 'None',
                                "options": ["true",'None'],
                                "help": "How to initialize the hessian if it is not fully provided. If true computes the Hessian from scratch and if expand, reads from a Hessian that was provided and interpolates if necessary."}),
               "hessian_update": (InputValue, {"dtype": str,
                                "default": "powell",
                                "options": ["powell", "recompute"],
                                "help": "How to update the hessian in each step. 'powell' uses the Powell method. 'recompute' recomputes the Hessian from scratch."}),
               "hessian_asr": (InputValue, {"dtype": str,
                                "default": "none",
                                "options": ["none","poly","crystal"],
                                "help": "Removes the zero frequency vibrational modes depending on the symmerty of the system."}),
               "final_rates": (InputValue, {"dtype": str,
                                "default": "false",
                                "options": ["false", "true"],
                                "help": "Decide if we are going to compute the Qvib and therefore the big-hessian by finite difference."}),
               "action": (InputArray, {"dtype": float,
                                        "default": input_default(factory=np.zeros, args=(0,)),
                                        "help": "Vector containing the 2  components ('spring' and 'physical') of the actions. Unitless "}),
               "prefix": (InputValue, {"dtype": str,
                                       "default": "INSTANTON",
                                       "help": "Prefix of the output files."
                                       })
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
        self.action.store(geop.action)
        self.prefix.store(geop.prefix)

    def fetch(self):
        rv = super(InputInst,self).fetch()
        rv["mode"] = self.mode.fetch()
        return rv
