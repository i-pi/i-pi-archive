"""Contains the class that deals with storing the state of a physical system.

Contains code used to hold the information which represents the state of
a system, including the particle positions and momenta, and the
forcefields which govern the interaction potential.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import os.path
import sys
import time

import numpy as np

from ipi.utils.depend import *
from ipi.utils.units  import *
from ipi.utils.prng   import *
from ipi.utils.io     import *
from ipi.utils.io.inputs.io_xml import *
from ipi.utils.messages import verbosity, info
from ipi.utils.softexit import softexit
from ipi.engine.atoms import *
from ipi.engine.cell import *
from ipi.engine.forces import Forces
from ipi.engine.properties import Properties, Trajectories


__all__ = ['System']


class System(dobject):
   """Physical system object.

   Contains all the phsyical information. Also handles stepping and output.

   Attributes:
      beads: A beads object giving the atom positions.
      cell: A cell object giving the system box.
      flist: A list of forcefield objects giving different ways to partially
         calculate the forces.
      fproto: A Forces object for calculating the total force for all the
         replicas.
      bproto: A Forces object for calculating the total bias for all the
         replicas.
      ensemble: An ensemble object giving the objects necessary for producing
         the correct ensemble.
      outputs: A list of output objects that should be printed during the run
      nm:  A helper object dealing with normal modes transformation
      properties: A property object for dealing with property output.
      trajs: A trajectory object for dealing with trajectory output.
      init: A class to deal with initializing the system.
      simul: The parent simulation object.
   """

   def __init__(self, init, beads, cell, force_proto, ensemble, nm, prefix="", bias_proto=[]):
      """Initialises System class.

      Args:
         init: A class to deal with initializing the system.
         beads: A beads object giving the atom positions.
         cell: A cell object giving the system box.
         force_proto: A forcefield object giving the force calculator for each
            replica of the system.
         bias_proto: A forcefield object giving the bias calculator for each
            replica of the system.
         ensemble: An ensemble object giving the objects necessary for
            producing the correct ensemble.
         nm: A class dealing with path NM operations.
         prefix: A string used to differentiate the output files of different
            systems.
      """

      info(" # Initializing system object ", verbosity.low )
      self.prefix = prefix
      self.init = init
      self.ensemble = ensemble
      self.beads = beads
      self.cell = cell
      self.nm = nm

      self.fproto = force_proto
      self.forces = Forces()

      self.bproto = bias_proto
      self.bias = Forces()


      self.properties = Properties()
      self.trajs = Trajectories()

   def bind(self, simul):
      """Calls the bind routines for all the objects in the system."""

      self.simul = simul # keeps a handle to the parent simulation object

      # binds important computation engines
      self.nm.bind(self.beads, self.ensemble)
      self.forces.bind(self.beads, self.cell, self.fproto, self.simul.fflist)

      self.bias.bind(self.beads, self.cell, self.bproto, self.simul.fflist)

      self.ensemble.bind(self.beads, self.nm, self.cell, self.forces, self.bias, self.prng)

      self.init.init_stage2(self)

      # binds output management objects
      self.properties.bind(self)
      self.trajs.bind(self)
