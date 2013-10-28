"""Contains the class that deals with the running of the simulation and
outputting the results.

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
   Simulation: Deals with running the simulation and outputting the results.
"""

__all__ = ['Simulation']

import numpy as np
import os.path, sys, time, threading
from copy import deepcopy
from ipi.utils.depend import *
from ipi.utils.units  import *
from ipi.utils.prng   import *
from ipi.utils.io     import *
from ipi.utils.io.io_xml import *
from ipi.utils.messages import verbosity, info, warning
from ipi.utils.softexit import softexit
from ipi.engine.atoms import *
from ipi.engine.cell import *
from ipi.engine.forces import Forces
from ipi.engine.beads import Beads
from ipi.engine.normalmodes import NormalModes
from ipi.engine.properties import Properties, Trajectories
from ipi.engine.outputs import CheckpointOutput

import objgraph

class Simulation(dobject):
   """Main simulation object.

   Contains all the references and the main dynamics loop. Also handles the
   initialisation and output.

   Attributes:
      prng: A random number generator object.
      fflist: A list of forcefield objects that can be called to compute energy and forces
      syslist: A list of physical systems
      tsteps: The total number of steps.
      ttime: The wall clock time (in seconds).
      outputs: A list of output objects that should be printed during the run
      paratemp: A helper object for parallel tempering simulations
      chk: A checkpoint object which is kept up-to-date in case of emergency exit
      rollback: If set to true, the state of the simulation at the start
         of the step will be output to a restart file rather than
         the current state of the simulation. This is because we cannot
         restart from half way through a step, only from the beginning of a
         step, so this is necessary for the trajectory to be continuous.

   Depend objects:
      step: The current simulation step.
   """

   def __init__(self, mode, syslist, fflist, outputs, prng, paratemp, step=0, tsteps=1000, ttime=0):
      """Initialises Simulation class.

      Args:
         mode: What kind of simulation is this
         syslist: A list of system objects
         fflist: A list of forcefield objects
         prng: A random number object.
         paratemp: A parallel tempering helper class.
         outputs: A list of output objects.
         init: A class to deal with initializing the simulation object.
         step: An optional integer giving the current simulation time step.
            Defaults to 0.
         tsteps: An optional integer giving the total number of steps. Defaults
            to 1000.
         ttime: The simulation running time. Used on restart, to keep a
            cumulative total.
      """

      info(" # Initializing simulation object ", verbosity.low )
      self.prng = prng
      self.mode = mode

      self.syslist = syslist
      for s in syslist:
         s.prng = self.prng # binds the system's prng to self prng
         s.init.init_stage1(s)

      if self.mode == "md" and len(syslist)>1:
         warning("Multiple systems will evolve independently in a '"+self.mode+"' simulation.")

      self.fflist = {}
      for f in fflist:
         self.fflist[f.name] = f

      self.outtemplate = outputs

      dset(self, "step", depend_value(name="step", value=step))
      self.tsteps = tsteps
      self.ttime = ttime
      self.paratemp = paratemp

      self.chk = None
      self.rollback = True

   def bind(self):
      """Calls the bind routines for all the objects in the simulation."""

      for s in self.syslist:
         # binds important computation engines
         s.bind(self)

      self.outputs = []
      for o in self.outtemplate:
         if type(o) is CheckpointOutput:    # checkpoints are output per simulation
            o.bind(self)
            self.outputs.append(o)
         else:   # properties and trajectories are output per system
            isys=0
            for s in self.syslist:   # create multiple copies
               no = deepcopy(o)
               if s.prefix != "":
                  no.filename = s.prefix+"_"+no.filename
               no.bind(s)
               self.outputs.append(no)
               isys+=1

      self.chk = CheckpointOutput("RESTART", 1, True, 0)
      self.chk.bind(self)

      if self.mode == "paratemp":
          self.paratemp.bind(self.syslist, self.prng)
          softexit.register_function(self.paratemp.softexit)

   def softexit(self):
      """Deals with a soft exit request.

      Tries to ensure that a consistent restart checkpoint is
      written out.
      """

      if self.step < self.tsteps:
         self.step += 1
      if not self.rollback:
         self.chk.store()

      print "WRITING CHECKPOINT", self.chk.status.extra
      self.chk.write(store=False)

   def run(self):
      """Runs the simulation.

      Does all the simulation steps, and outputs data to the appropriate files
      when necessary. Also deals with starting and cleaning up the threads used
      in the communication between the driver and the PIMD code.
      """

      # registers the softexit routine
      softexit.register_function(self.softexit)
      softexit.start(self.ttime)

      for (k,f) in self.fflist.iteritems():
         f.run()

      # prints inital configuration -- only if we are not restarting
      if (self.step == 0):
         self.step = -1
         for o in self.outputs:
            o.write()
         if self.mode == "paratemp":
            self.paratemp.parafile.write("%10d" % (self.step+1))
            for i in self.paratemp.temp_index:
               self.paratemp.parafile.write(" %5d" %i)
            self.paratemp.parafile.write("\n")
            if self.paratemp.wtefile != None:
               self.paratemp.wtefile.write("%10d" % (self.step+1))
               for v in self.paratemp.system_v:
                  self.paratemp.wtefile.write(" %12.7e" % v)
               self.paratemp.wtefile.write("\n")


         self.step = 0

      steptime = 0.0
      simtime =  time.time()

      cstep = 0
      tptime = 0.0
      tqtime = 0.0
      tttime = 0.0
      ttot = 0.0
      # main MD loop
      for self.step in range(self.step,self.tsteps):
         # stores the state before doing a step.
         # this is a bit time-consuming but makes sure that we can honor soft
         # exit requests without screwing the trajectory

         steptime = -time.time()
         if softexit.triggered: break

         self.chk.store()

         stepthreads = []
         # steps through all the systems
         for s in self.syslist:
            # creates separate threads for the different systems
            st = threading.Thread(target=s.ensemble.step, name=s.prefix)
            st.daemon = True
            st.start()
            stepthreads.append(st)

         for st in stepthreads:
            st.join()

         if softexit.triggered: break # don't continue if we are about to exit!

         for o in self.outputs:
            o.write()

         if self.mode == "paratemp": # does parallel tempering and/or WTE

            self.paratemp.wtestep(self.step)

            # because of where this is in the loop, we must write out BEFORE doing the swaps.
            self.paratemp.parafile.write("%10d" % (self.step+1))
            for i in self.paratemp.temp_index:
               self.paratemp.parafile.write(" %5d" %i)
            self.paratemp.parafile.write("\n")

            # applies the WTE forces, if they are defined.
            if self.paratemp.wtefile != None:
               self.paratemp.wtefile.write("%10d" % (self.step+1))
               for v in self.paratemp.system_v:
                  self.paratemp.wtefile.write(" %12.7e" % v)
               self.paratemp.wtefile.write("\n")

            self.paratemp.swap(self.step)
            self.paratemp.wtestep(self.step)


         if softexit.triggered: break # don't write if we are about to exit!

         steptime += time.time()
         ttot += steptime
         cstep += 1

         if verbosity.high or (verbosity.medium and self.step%100 == 0) or (verbosity.low and self.step%1000 == 0):
            info(" # Average timings at MD step % 7d. t/step: %10.5e" %
               ( self.step, ttot/cstep ) )
            cstep = 0
            ttot = 0.0
         #   info(" # MD diagnostics: V: %10.5e    Kcv: %10.5e   Ecns: %10.5e" %
         #      (self.properties["potential"], self.properties["kinetic_cv"], self.properties["conserved"] ) )

         if os.path.exists("EXIT"): # soft-exit
            info(" # EXIT file detected! Bye bye!", verbosity.low )
            break

         if (self.ttime > 0 and time.time() - simtime > self.ttime):
            info(" # Wall clock time expired! Bye bye!", verbosity.low )
            break

         objgraph.show_growth(limit=3)

      self.rollback = False
      softexit.trigger(" @ SIMULATION: Exiting cleanly.")
