"""Contains the class that deals with storing the state of a physical system.

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


The root class for the whole simulation. Contains references to all the top
level objects used in the simulation, and controls all the steps that are
not inherently system dependent, like the running of each time step,
choosing which properties to initialise, and which properties to output.

Classes:
   System: Deals with storing and outputting information on a physical system.
"""

__all__ = ['System']

import numpy as np
import os.path, sys, time
from ipi.utils.depend import *
from ipi.utils.units  import *
from ipi.utils.prng   import *
from ipi.utils.io     import *
from ipi.utils.io.io_xml import *
from ipi.utils.messages import verbosity, info
from ipi.utils.softexit import softexit
from ipi.engine.atoms import *
from ipi.engine.cell import *
from ipi.engine.forces import Forces
from ipi.engine.beads import Beads
from ipi.engine.normalmodes import NormalModes
from ipi.engine.properties import Properties, Trajectories

class System(dobject):
   """Physical system object.

   Contains all the phsyical information. Also handles stepping and output.

   Attributes:
      beads: A beads object giving the atom positions.
      cell: A cell object giving the system box.
      prng: A random number generator object.
      flist: A list of forcefield objects giving different ways to partially
         calculate the forces.
      forces: A Forces object for calculating the total force for all the
         replicas.
      ensemble: An ensemble object giving the objects necessary for producing
         the correct ensemble.
      outputs: A list of output objects that should be printed during the run
      nm:  A helper object dealing with normal modes transformation
      properties: A property object for dealing with property output.
      trajs: A trajectory object for dealing with trajectory output.

   Depend objects:
      step: The current simulation step.
   """

   def __init__(self, beads, cell, forces, ensemble, nm):
      """Initialises System class.

      Args:
         beads: A beads object giving the atom positions.
         cell: A cell object giving the system box.
         forces: A forcefield object giving the force calculator for each
            replica of the system.
         ensemble: An ensemble object giving the objects necessary for
            producing the correct ensemble.
         nm: A class dealing with path NM operations.
      """

      info(" # Initializing system object ", verbosity.low )
      self.ensemble = ensemble
      self.beads = beads
      self.cell = cell
      self.nm = nm

      self.flist = forces
      self.forces = Forces()

      self.properties = Properties()
      self.trajs = Trajectories()

   def bind(self, simul):
      """Calls the bind routines for all the objects in the system."""

      self.simul = simul # keeps a handle to the parent simulation object

      # binds important computation engines
      self.nm.bind(self.beads, self.ensemble)
      self.forces.bind(self.beads, self.cell, self.flist)
      self.ensemble.bind(self.beads, self.nm, self.cell, self.forces, self.prng)

      # binds output management objects
      self.properties.bind(self)
      self.trajs.bind(self)
