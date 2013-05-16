"""Contains the class that deals with the running of the simulation and
outputting the results.

The root class for the whole simulation. Contains references to all the top
level objects used in the simulation, and controls all the steps that are
not inherently system dependent, like the running of each time step,
choosing which properties to initialise, and which properties to output.

Classes:
   Simulation: Deals with running the simulation and outputting the results.
"""

__all__ = ['Simulation']

import numpy as np
import math
import os.path, sys
from utils.depend import *
from utils.units  import *
from utils.prng   import *
from utils.io     import *
from utils.io.io_xml import *
from utils.messages import verbosity, info
from atoms import *
import time
from cell import *
from forces import Forces
from beads import Beads
from normalmodes import NormalModes
from properties import Properties, Trajectories
from outputs import CheckpointOutput
import inputs.simulation

class Simulation(dobject):
   """Main simulation object.

   Contains all the references and the main dynamics loop. Also handles the
   initialisation and output.

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
      tsteps: The total number of steps.
      ttime: The wall clock time (in seconds).
      format: A string specifying both the format and the extension of traj
         output.
      outputs: A list of output objects that should be printed during the run
      nm:  A helper object dealing with normal modes transformation
      properties: A property object for dealing with property output.
      trajs: A trajectory object for dealing with trajectory output.
      chk: A checkpoint object for dealing with checkpoint output.

   Depend objects:
      step: The current simulation step.
   """

   def __init__(self, beads, cell, forces, ensemble, prng, outputs, nm, init, step=0, tsteps=1000, ttime=0):
      """Initialises Simulation class.

      Args:
         beads: A beads object giving the atom positions.
         cell: A cell object giving the system box.
         forces: A forcefield object giving the force calculator for each
            replica of the system.
         ensemble: An ensemble object giving the objects necessary for
            producing the correct ensemble.
         prng: A random number object.
         outputs: A list of output objects.
         nm: A class dealing with path NM operations.
         init: A class to deal with initializing the simulation object.
         step: An optional integer giving the current simulation time step.
            Defaults to 0.
         tsteps: An optional integer giving the total number of steps. Defaults
            to 1000.
      """

      info(" # Initializing simulation object ", verbosity.low )

      self.prng = prng
      self.ensemble = ensemble
      self.beads = beads
      self.cell = cell
      self.nm = nm

      # initialize the configuration of the system
      init.init(self)

      self.flist = forces
      self.forces = Forces()
      self.outputs = outputs

      dset(self, "step", depend_value(name="step", value=step))
      self.tsteps = tsteps
      self.ttime = ttime

      self.properties = Properties()
      self.trajs = Trajectories()
      self.chk = None


   def bind(self):
      """Calls the bind routines for all the objects in the simulation."""

      # binds important computation engines
      self.nm.bind(self.beads, self.ensemble)
      self.forces.bind(self.beads, self.cell, self.flist, softexit=self.soft_exit)
      self.ensemble.bind(self.beads, self.nm, self.cell, self.forces, self.prng)

      # binds output management objects
      self.properties.bind(self)
      self.trajs.bind(self)
      for o in self.outputs:
         o.bind(self)

      self.chk = CheckpointOutput("RESTART", 1, True, 0)
      self.chk.bind(self)

   def soft_exit(self, rollback=True):
      """Deals with a soft exit request.

      Tries to ensure that a consistent restart checkpoint is
      written out.
      """

      if self.step < self.tsteps:
         self.step += 1
      if not rollback:
         self.chk.store()
      self.chk.write(store=False)

      self.forces.stop()

      sys.exit()

   def run(self):
      """Runs the simulation.

      Does all the simulation steps, and outputs data to the appropriate files
      when necessary. Also deals with starting and cleaning up the threads used
      in the communication between the driver and the PIMD code.
      """

      self.forces.run()

      # prints inital configuration -- only if we are not restarting
      if (self.step == 0):
         self.step = -1
         for o in self.outputs:
            o.write()
         self.step = 0

      steptime = 0.0
      simtime =  time.time()
      # main MD loop
      for self.step in range(self.step,self.tsteps):
         # stores the state before doing a step.
         # this is a bit time-consuming but makes sure that we can honor soft
         # exit requests without screwing the trajectory

         steptime = -time.time()
         self.chk.store()

         self.ensemble.step()

         for o in self.outputs:
            o.write()

         if os.path.exists("EXIT"): # soft-exit
            self.soft_exit(rollback=False)

         steptime += time.time()
         print " # MD step % 7d complete. Timings -->  %10.5e [p: %10.5e  q: %10.5e  t: %10.5e]" % (self.step, steptime, self.ensemble.ptime, self.ensemble.qtime, self.ensemble.ttime )
         print " # MD diagnostics: V: %10.5e    Kcv: %10.5e   Ecns: %10.5e" % (self.properties["potential"], self.properties["kinetic_cv"], self.properties["conserved"] )

         if (self.ttime > 0 and time.time() - simtime > self.ttime):
            print " # Wall clock time expired! Bye bye"
            break

      self.soft_exit(rollback=False)
