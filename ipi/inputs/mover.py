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

__all__ = ['InputMover', 'InputGEOP']

class InputGEOP(Input):
    """Geometry optimization options.
    
    Contains options related with geometry optimization, such as method, 
    thresholds, linear search strategy, etc. 

    """

    attribs={"mode"  : (InputAttribute, {"dtype"   : str,
                                    "help"    : "The geometry optimization algorithm to be used",
                                    "options" : ['steepest']}) }
   
    fields = { "line_tolerance": (InputValue, {"dtype"         : float,
                                     "default"       : 1.0e-5,
                                     "help"          : "The tolerance for line search procedures."})  }
              
    dynamic = {  }

    default_help = "TODO EXPLAIN WHAT THIS IS"
    default_label = "GEOP"   

    def store(self, geop):
        self.line_tolerance.store(geop.tol) # todo set the proper parameter
        self.mode.store(geop.mode)
		
    def fetch(self):		
        ngeo = Ensemble() # this should be a geop object we don't have here
        ngeo.tol = self.line_tolerance.fetch()


		
class InputMover(Input):
   """Mover calculation input class.

   A class to encompass the different "mover" (non-MD) calculations. 

   Attributes:
      mode: An optional string giving the kind of mover calculation to be performed.

   Fields:      
      fixcom: An optional boolean which decides whether the centre of mass
         motion will be constrained or not. Defaults to False.
      fixatoms: A list of the indices of atoms that should not be moved.
      
   """

   attribs={"mode"  : (InputAttribute, {"dtype"   : str,
                                    "help"    : "The ensemble hat will be sampled during the simulation. 'replay' means that a simulation is restarted from a previous simulation.",
                                    "options" : ['minimize', 'replay', 'dummy']}) }
   fields={"fixcom": (InputValue, {"dtype"           : bool,
                                   "default"         : True,
                                   "help"            : "This describes whether the centre of mass of the particles is fixed."}),
           "fixatoms" : (InputArray, {"dtype"        : int,
                                    "default"      : np.zeros(0,int),
                                    "help"         : "Indices of the atmoms that should be held fixed."}),
           "replay_file": (InputInitFile, {"default" : input_default(factory=ipi.engine.initializer.InitBase),
                           "help"            : "This describes the location to read a trajectory file from."})
         }
         
   dynamic = {  }

   default_help = "Holds all the information that is calculation specific, such as geometry optimization parameters, etc."
   default_label = "STATIC"

   def store(self, sc):
      """Takes a mover calculation instance and stores a minimal representation of it.

      Args:
         sc: A mover calculation class.
      """

      super(InputMover,self).store(sc)    
      tsc = -1
      if type(sc) is Mover:
          self.mode.store("dummy")
      elif type(sc) is ReplayMover:
         self.mode.store("replay")
         tsc = 0
      elif type(sc) is GeopMover:
         self.mode.store("minimize")
         tsc = 1
      else: 
         raise ValueError("Cannot store Mover calculator of type "+str(type(sc)))
      
      if tsc == 0:
         self.replay_file.store(sc.intraj)
      if tsc > 0:
         self.fixcom.store(sc.fixcom)
         self.fixatoms.store(sc.fixatoms)

   def fetch(self):
      """Creates an mover calculator object.

      Returns:
         An ensemble object of the appropriate mode and with the appropriate
         objects given the attributes of the InputEnsemble object.
      """

      super(InputMover,self).fetch()

      if self.mode.fetch() == "replay" :
         sc = ReplayMover(fixcom=False, fixatoms=None, intraj=self.replay_file.fetch() )
      elif self.mode.fetch() == "minimize":
         sc = GeopMover(fixcom=False, fixatoms=None )
      else:
         raise ValueError("'" + self.mode.fetch() + "' is not a supported mover calculation mode.")

      return sc
