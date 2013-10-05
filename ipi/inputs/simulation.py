"""Deals with creating the simulation class.

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
   InputSimulation: Deals with creating the Simulation object from a file, and
      writing the checkpoints.
"""

__all__ = ['InputSimulation']

import numpy as np
import os.path, sys
import ipi.engine.simulation
from ipi.utils.depend import *
from ipi.utils.inputvalue import *
from ipi.utils.units  import *
from ipi.utils.prng   import *
from ipi.utils.io     import *
from ipi.utils.io.io_xml import *
from ipi.utils.messages import verbosity
from ipi.inputs.prng import InputRandom
from ipi.inputs.system import InputSystem
from ipi.inputs.forcefields import InputFFSocket
from ipi.inputs.outputs import InputOutputs

class InputSimulation(Input):
   """Simulation input class.

   Handles generating the appropriate forcefield class from the xml input file,
   and generating the xml checkpoint tags and data from an instance of the
   object.

   Attributes:
      verbosity: A string saying how much should be output to standard output.

   Fields:
      force: A restart force instance. Used as a model for all the replicas.
      ensemble: A restart ensemble instance.
      beads: A restart beads instance.
      normal_modes: Setup of normal mode integrator.
      cell: A restart cell instance.
      output: A list of the required outputs.
      prng: A random number generator object.
      step: An integer giving the current simulation step. Defaults to 0.
      total_steps: The total number of steps. Defaults to 1000
      total_time:  The wall clock time limit. Defaults to 0 (no limit).
      initialize: An array of strings giving all the quantities that should
         be output.
   """

   fields = {
             "prng" :    (InputRandom,   { "help"  : InputRandom.default_help,
                                         "default" : input_default(factory=Random)} ),
             "output" :  (InputOutputs, { "help"   : InputOutputs.default_help,
                                          "default": input_default(factory=InputOutputs.make_default)  }),
             "step" :       ( InputValue, { "dtype"    : int,
                                            "default"  : 0,
                                            "help"     : "The current simulation time step." }),
             "total_steps": ( InputValue, { "dtype"    : int,
                                            "default"  : 1000,
                                            "help"     : "The total number of steps that will be done. If 'step' is equal to or greater than 'total_steps', then the simulation will finish." }),
             "total_time" :       ( InputValue, { "dtype"    : float,
                                            "default"  : 0,
                                            "help"     : "The maximum wall clock time (in seconds)." }),
                                             }

   attribs = { "verbosity" : (InputAttribute, { "dtype"   : str,
                                      "default" : "low",
                                      "options" : [ "quiet", "low", "medium", "high", "debug" ],
                                      "help"    : "The level of output on stdout."
                                         }),
               "mode"  : (InputAttribute, {"dtype"   : str,
                                    "default" : "md",
                                    "help"    : "What kind of simulation should be run.",
                                    "options" : ['md', 'paratemp']})
             }

   dynamic = {
             "system" :   (InputSystem,    { "help"  : InputSystem.default_help }), 
             "ffsocket": (InputFFSocket, { "help": InputFFSocket.default_help} )
             }

   default_help = "This is the top level class that deals with the running of the simulation, including holding the simulation specific properties such as the time step and outputting the data."
   default_label = "SIMULATION"

   def store(self, simul):
      """Takes a simulation instance and stores a minimal representation of it.

      Args:
         simul: A simulation object.
      """

      super(InputSimulation,self).store()


      self.output.store(simul.outtemplate)
      self.prng.store(simul.prng)
      self.step.store(simul.step)
      self.total_steps.store(simul.tsteps)
      self.total_time.store(simul.ttime)

      # this we pick from the messages class. kind of a "global" but it seems to
      # be the best way to pass around the (global) information on the level of output.
      if verbosity.debug:
         self.verbosity.store("debug")
      elif verbosity.high:
         self.verbosity.store("high")
      elif verbosity.medium:
         self.verbosity.store("medium")
      elif verbosity.low:
         self.verbosity.store("low")
      elif verbosity.quiet:
         self.verbosity.store("quiet")
      else:
         raise ValueError("Invalid verbosity level")

      self.extra = []

      for f in simul.fflist:
         iff = InputFFSocket()         
         iff.store(simul.fflist[f])
         self.extra.append(("ffsocket",iff))
      for s in simul.syslist:
         isys = InputSystem()
         isys.store(s)
         self.extra.append(("system",isys))


   def fetch(self):
      """Creates a simulation object.

      Returns:
         A simulation object of the appropriate type and with the appropriate
         properties and other objects given the attributes of the
         InputSimulation object.

      Raises:
         TypeError: Raised if one of the file types in the stride keyword
            is incorrect.
      """

      super(InputSimulation,self).fetch()

      # small hack: initialize here the verbosity level -- we really assume to have
      # just one simulation object
      verbosity.level=self.verbosity.fetch()

      syslist=[]
      fflist=[]
      for (k,v) in self.extra:
         if k=="system" : 
            for isys in range(v.copies.fetch()): # creates multiple copies of system if desired               
               syslist.append(v.fetch())
               if (v.copies.fetch() > 1):
                  syslist[-1].prefix = syslist[-1].prefix+( ("%0" + str(int(1 + np.floor(np.log(v.copies.fetch())/np.log(10)))) + "d") % (isys) )
         elif k=="ffsocket": fflist.append(v.fetch())


      # this creates a simulation object which gathers all the little bits
      #TODO use named arguments since this list is a bit too long...
      rsim = ipi.engine.simulation.Simulation(
                  syslist, fflist,
                  self.output.fetch(),
                  self.prng.fetch(),
                      self.step.fetch(),
                        tsteps=self.total_steps.fetch(),
                           ttime=self.total_time.fetch())

      # this does all of the piping between the components of the simulation
      rsim.bind()

      return rsim

   def check(self):
      """Function that deals with optional arguments.

      Deals with the difference between classical and PI dynamics. If there is
      no beads argument, the bead positions are generated from the atoms, with
      the necklace being fixed at the atom position. Similarly, if no nbeads
      argument is specified a classical simulation is done.

      Raises:
         TypeError: Raised if no beads or atoms attribute is defined.
      """

      super(InputSimulation,self).check()
      if self.total_steps.fetch() <= self.step.fetch():
         raise ValueError("Current step greater than total steps, no dynamics will be done.")
