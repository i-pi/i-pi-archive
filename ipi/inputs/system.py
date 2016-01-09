"""Deals with creating a representation of a system."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import os.path
import sys

import numpy as np

import ipi.engine.system
from ipi.utils.depend import *
from ipi.utils.inputvalue import *
from ipi.utils.units  import *
from ipi.utils.prng   import *
from ipi.utils.io     import *
from ipi.utils.io.inputs.io_xml import *
from ipi.utils.messages import verbosity
from ipi.inputs.forces import InputForces
from ipi.inputs.beads import InputBeads
from ipi.inputs.cell import InputCell
from ipi.inputs.ensembles import InputEnsemble
from ipi.inputs.motion import InputMotion
from ipi.inputs.normalmodes import InputNormalModes
from ipi.engine.normalmodes import NormalModes
from ipi.engine.atoms import Atoms
from ipi.engine.beads import Beads
from ipi.engine.cell import Cell
from ipi.engine.forces import Forces
from ipi.engine.ensembles import Ensemble
from ipi.engine.motion import Motion
from ipi.inputs.initializer import InputInitializer
from ipi.engine.initializer import Initializer


__all__ = ['InputSystem']


class InputSystem(Input):
   """Physical system input class.

   Handles generating the appropriate forcefield class from the xml input file,
   and generating the xml checkpoint tags and data from an instance of the
   object.

   Attributes:
      copies: Decides how many of each system to create.
      prefix: A string to prepend to the output file names for this system.

   Fields:
      forces: A restart force instance. Used as a model for all the replicas.
      ensemble: A restart ensemble instance.
      beads: A restart beads instance.
      normal_modes: Setup of normal mode integrator.
      cell: A restart cell instance.
      initialize: An array of strings giving all the quantities that should
         be output.
   """

   fields = {
             "initialize" : (InputInitializer, { "help" : InputInitializer.default_help,
                                                "default" : input_default(factory=Initializer) } ),
             "forces" :   (InputForces,    { "help"  : InputForces.default_help }),
             "bias" :   (InputForces,    { "help"  : InputForces.default_help,
                                           "default" : [] }),
             "ensemble": (InputEnsemble, { "help"  : InputEnsemble.default_help ,
                             "default" : input_default(factory=Ensemble, kwargs={'temp':1.0})} ),
             "motion": (InputMotion, { "help"  : InputMotion.default_help, "default" : input_default(factory=Motion) } ),
             "beads" :   (InputBeads, { "help"     : InputBeads.default_help,
                                        "default"  : input_default(factory=Beads, kwargs={'natoms': 0, 'nbeads': 0}) } ),
             "normal_modes" :   (InputNormalModes, { "help"     : InputNormalModes.default_help,
                                        "default"  : input_default(factory=NormalModes, kwargs={'mode': "rpmd"}) } ),
             "cell" :    (InputCell,   { "help"    : InputCell.default_help,
                                        "default"  : input_default(factory=Cell) })
             }

   attribs = {
    "copies": (InputAttribute, {"help" : "Create multiple copies of the system. This is handy for initialising simulations with multiple systems.", "default": 1, "dtype": int}) ,
    "prefix": (InputAttribute, {"help" : "Prepend this string to output files generated for this system. If 'copies' is greater than 1, a trailing number will be appended.", "default": "", "dtype": str})
   }

   default_help = "This is the class which holds all the data which represents a single state of the system."
   default_label = "SYSTEM"

   def store(self, psys):
      """Takes a System instance and stores a minimal representation of it.

      Args:
         psys: A physical system object.
      """

      super(InputSystem,self).store()

      self.prefix.store(psys.prefix)
      self.forces.store(psys.fcomp)
      self.bias.store(psys.bcomp)
      self.ensemble.store(psys.ensemble)
      self.motion.store(psys.motion)
      self.beads.store(psys.beads)
      self.normal_modes.store(psys.nm)
      self.cell.store(psys.cell)

   def fetch(self):
      """Creates a physical system object.

      Returns:
         A System object of the appropriate type and with the appropriate
         properties and other objects given the attributes of the
         InputSystem object.

      Raises:
         TypeError: Raised if one of the file types in the stride keyword
            is incorrect.
      """

      super(InputSystem,self).fetch()

      # this creates a simulation object which gathers all the little bits
      #TODO use named arguments since this list is a bit too long...
      rsys = ipi.engine.system.System( init = self.initialize.fetch(),
                                       beads = self.beads.fetch(),
                                       nm = self.normal_modes.fetch(),
                                       cell = self.cell.fetch(),
                                       fcomponents = self.forces.fetch(),
                                       bcomponents = self.bias.fetch(),
                                       ensemble = self.ensemble.fetch(),
                                       motion = self.motion.fetch(),
                                       prefix = self.prefix.fetch()
                                       )

      return rsys
