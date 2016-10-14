"""Creates objects that hold the whole simulation."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import os.path
import sys
import time

import numpy as np

from ipi.utils.depend import *
from ipi.utils.inputvalue import *
from ipi.utils.units  import *
from ipi.utils.prng   import *
from ipi.utils.io     import *
from ipi.utils.io.inputs.io_xml import *
from ipi.utils.messages import verbosity
from ipi.engine.paratemp import ParaTemp
from ipi.inputs.prng import InputRandom
from ipi.inputs.system import InputSystem
import ipi.inputs.forcefields as iforcefields
import ipi.engine.forcefields as eforcefields
import ipi.inputs.outputs as ioutputs
from ipi.inputs.paratemp import InputParaTemp


__all__ = ['InputSimulation']


class InputSimulation(Input):
   """Simulation input class.

   Handles generating the appropriate forcefield class from the xml input file,
   and generating the xml checkpoint tags and data from an instance of the
   object.

   Attributes:
      verbosity: A string saying how much should be output to standard output.
      mode: A string which determines what type of simulation will be run.

   Fields:
      output: A list of the required outputs.
      prng: A random number generator object.
      step: An integer giving the current simulation step. Defaults to 0.
      total_steps: The total number of steps. Defaults to 1000
      total_time:  The wall clock time limit. Defaults to 0 (no limit).
      paratemp: A helper object for parallel tempering simulations

   Dynamic fields:
      system: Holds the data needed to specify the state of a single system.
      ffsocket: Gives a forcefield which will use a socket interface to
         communicate with the driver code.
      fflj: Gives a forcefield which uses the internal Python Lennard-Jones
         script to calculate the potential and forces.
   """

   fields = {
             "prng" :    (InputRandom,   { "help"  : InputRandom.default_help,
                                         "default" : input_default(factory=Random)} ),
             "output" :  (ioutputs.InputOutputs, { "help"   : ioutputs.InputOutputs.default_help,
                                          "default": input_default(factory=ioutputs.InputOutputs.make_default)  }),
             "step" :       ( InputValue, { "dtype"    : int,
                                            "default"  : 0,
                                            "help"     : "The current simulation time step." }),
             "total_steps": ( InputValue, { "dtype"    : int,
                                            "default"  : 1000,
                                            "help"     : "The total number of steps that will be done. If 'step' is equal to or greater than 'total_steps', then the simulation will finish." }),
             "total_time" :       ( InputValue, { "dtype"    : float,
                                            "default"  : 0,
                                            "help"     : "The maximum wall clock time (in seconds)." }),
             "paratemp" : (InputParaTemp, {"default"   : input_default(factory=ParaTemp),
                                         "help"      : "Options for a parallel tempering simulation"})
            }

   attribs = { "verbosity" : (InputAttribute, { "dtype"   : str,
                                      "default" : "low",
                                      "options" : [ "quiet", "low", "medium", "high", "debug" ],
                                      "help"    : "The level of output on stdout."
                                         }),
               "mode"  : (InputAttribute, {"dtype"   : str,
                                    "default" : "md",
                                    "help"    : "What kind of simulation should be run.",
                                    "options" : ['md', 'paratemp', 'static']})
             }

   dynamic = {
             "system" :   (InputSystem,    { "help"  : InputSystem.default_help }),
             "ffsocket": (iforcefields.InputFFSocket, { "help": iforcefields.InputFFSocket.default_help} ),
             "fflj": (iforcefields.InputFFLennardJones, { "help": iforcefields.InputFFLennardJones.default_help} )
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
      self.paratemp.store(simul.paratemp)

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

      self.mode.store(simul.mode)

      self.extra = []

      for fname in simul.fflist:
         ff=simul.fflist[fname]
         if type(ff) is eforcefields.FFSocket:
            iff = iforcefields.InputFFSocket()
            iff.store(ff)
            self.extra.append(("ffsocket",iff))
         elif type(ff) is eforcefields.FFLennardJones:
            iff = iforcefields.InputFFLennardJones()
            iff.store(ff)
            self.extra.append(("fflj",iff))


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
         if k == "system":
            for isys in range(v.copies.fetch()): # creates multiple copies of system if desired
               syslist.append(v.fetch())
               if (v.copies.fetch() > 1):
                  syslist[-1].prefix = syslist[-1].prefix + ( ("%0" + str(int(1 + np.floor(np.log(v.copies.fetch())/np.log(10)))) + "d") % (isys) )
         elif k == "ffsocket":
            fflist.append(v.fetch())
         elif k == "fflj":
            fflist.append(v.fetch())


      # this creates a simulation object which gathers all the little bits
      import ipi.engine.simulation as esimulation   # import here as otherwise this is the mother of all circular imports...
      rsim = esimulation.Simulation(
                  mode = self.mode.fetch(),
                  syslist = syslist,
                  fflist = fflist,
                  outputs = self.output.fetch(),
                  prng = self.prng.fetch(),
                  paratemp = self.paratemp.fetch(),
                  step = self.step.fetch(),
                  tsteps = self.total_steps.fetch(),
                  ttime = self.total_time.fetch())

      return rsim
