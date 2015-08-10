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
from ipi.engine.geop import GeopMover
from ipi.engine.neb import NEBMover
from ipi.engine.mover import *
from ipi.utils.inputvalue import *
from ipi.inputs.thermostats import *
from ipi.inputs.initializer import *
from ipi.inputs.geop import InputGeop
from ipi.inputs.neb import InputNEB
from ipi.utils.units import *

__all__ = ['InputMover']
		
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
                                    "options" : ['minimize', 'replay', 'neb', 'dummy']}) }
   fields={"fixcom": (InputValue, {"dtype"           : bool,
                                   "default"         : True,
                                   "help"            : "This describes whether the centre of mass of the particles is fixed."}),
           "fixatoms" : (InputArray, {"dtype"        : int,
                                    "default"      : np.zeros(0,int),
                                    "help"         : "Indices of the atmoms that should be held fixed."}),
           "optimizer" : ( InputGeop, { "default" : {}, 
                                     "help":  "Option for geometry optimization" } ),
           "neb_optimizer" : ( InputGeop, { "default" : {}, 
                                     "help":  "Option for geometry optimization" } ),
           "dynamics" : ( InputGeop, { "default" : {}, 
                                     "help":  "Option for (path integral) molecular dynamics" } ),                          
           "file": (InputInitFile, {"default" : input_default(factory=ipi.engine.initializer.InitBase,kwargs={"mode":"xyz"}),
                           "help"            : "This describes the location to read a trajectory file from."})
         }
         
   dynamic = {  }

   default_help = "Holds all the information that is calculation specific, such as geometry optimization parameters, etc."
   default_label = "MOVER"

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
         self.optimizer.store(sc)
         tsc = 1
      elif type(sc) is NEBMover:
         self.mode.store("neb")
         self.neb_optimizer.store(sc)
         tsc = 1
      elif type(sc) is DynMover:
         self.mode.store("dynamics")
         self.dynamics.store(sc)
         tsc = 1   
      else: 
         raise ValueError("Cannot store Mover calculator of type "+str(type(sc)))
      
      if tsc == 0:
         self.file.store(sc.intraj)
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
         sc = ReplayMover(fixcom=False, fixatoms=None, intraj=self.file.fetch() )
      elif self.mode.fetch() == "minimize":
         sc = GeopMover(fixcom=False, fixatoms=None, **self.optimizer.fetch() )
      elif self.mode.fetch() == "neb":
         sc = NEBMover(fixcom=False, fixatoms=None, **self.neb_optimizer.fetch() )
      elif self.mode.fetch() == "dynamics":
         sc = DynMover(fixcom=False, fixatoms=None, **self.dynamics.fetch() )
      else:
         sc = Mover()
         #raise ValueError("'" + self.mode.fetch() + "' is not a supported mover calculation mode.")

      return sc
