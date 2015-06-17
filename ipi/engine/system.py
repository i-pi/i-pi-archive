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


Contains code used to hold the information which represents the state of
a system, including the particle positions and momenta, and the
forcefields which govern the interaction potential.

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
from ipi.engine.properties import Properties, Trajectories

class System(dobject):
   """Physical system object.

   Contains all the phsyical information. Also handles stepping and output.

   Attributes:
      beads: A beads object giving the atom positions.
      cell: A cell object giving the system box.
      fcomp: A list of force components that must act on each replica
      bcomp: A list of bias components that must act on each replica
      forces: A Forces object that actually compute energy and forces
      bias: A Forces object that compute the bias components
      ensemble: An ensemble object giving the objects necessary for producing
         the correct ensemble.
      outputs: A list of output objects that should be printed during the run
      nm:  A helper object dealing with normal modes transformation
      properties: A property object for dealing with property output.
      trajs: A trajectory object for dealing with trajectory output.
      init: A class to deal with initializing the system.
      simul: The parent simulation object.
   """

   def __init__(self, init, beads, nm, cell, fcomponents, bcomponents=[], ensemble=None, mover=None, prefix=""):
      """Initialises System class.

      Args:
         init: A class to deal with initializing the system.
         beads: A beads object giving the atom positions.
         cell: A cell object giving the system box.
         fcomponents: A list of force components that are active for each
            replica of the system.
         bcomponents: A list of force components that are considered as bias, and act on each 
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
      self.mover = mover
      self.beads = beads
      self.cell = cell
      self.nm = nm

      self.fcomp = fcomponents
      self.forces = Forces()

      self.bcomp = bcomponents
      self.bias = Forces()


      self.properties = Properties()
      self.trajs = Trajectories()

   def bind(self, simul):
      """Calls the bind routines for all the objects in the system."""

      self.simul = simul # keeps a handle to the parent simulation object

      # binds important computation engines
      self.nm.bind(self.beads, self.ensemble)
      self.forces.bind(self.beads, self.cell, self.fcomp, self.simul.fflist)
      
      self.bias.bind(self.beads, self.cell, self.bcomp, self.simul.fflist)
      
      self.ensemble.bind(self.beads, self.nm, self.cell, self.forces, self.bias, self.prng)

      self.init.init_stage2(self)

      # binds output management objects
      self.properties.bind(self)
      self.trajs.bind(self)
