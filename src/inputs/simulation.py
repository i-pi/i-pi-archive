"""Deals with creating the simulation class.

Classes:
   InputSimulation: Deals with creating the Simulation object from a file, and
      writing the checkpoints.
"""

__all__ = ['InputSimulation']

import numpy as np
import math, random
import os.path, sys
from utils.depend import *
from utils.inputvalue import *
from utils.units  import *
from utils.prng   import *
from utils.io     import *
from utils.io.io_xml import *
from atoms import *
from cell import *
from inputs.forces import InputForce
from inputs.prng import InputRandom
from inputs.initializer import InputInitializer
from inputs.beads import InputBeads
from inputs.cell import InputCell
from inputs.ensembles import InputEnsemble
from inputs.outputs import InputOutputs
from inputs.normalmodes import InputNormalModes
from engine.normalmodes import NormalModes
from engine.atoms import Atoms
from engine.beads import Beads
from engine.initializer import Initializer

import engine.outputs
import engine.simulation

class InputSimulation(Input):
   """Simulation input class.

   Handles generating the appropriate forcefield class from the xml input file,
   and generating the xml checkpoint tags and data from an instance of the
   object.

   Attributes:
      force: A restart force instance. Used as a model for all the replicas.
      ensemble: A restart ensemble instance.
      atoms: A restart atoms instance.
      beads: A restart beads instance.
      normal_modes: Setup of normal mode integrator.
      cell: A restart cell instance.
      prng: A random number generator object.
      step: An integer giving the current simulation step. Defaults to 0.
      total_steps: The total number of steps. Defaults to 0.
      stride: A dictionary giving the number of steps between printing out
         data for the different types of data. Defaults to _DEFAULT_STRIDES.
      traj_format: A string giving the format of the trajectory output files.
         Defaults to 'pdb'.
      trajectories: An array of strings giving all the trajectory data that
         should be output space separated. Defaults to _DEFAULT_TRAJ.
      initialize: An array of strings giving all the quantities that should
         be output.
      fd_delta: A float giving the size of the finite difference
         parameter used in the Yamamoto kinetic energy estimator. Defaults
         to 0.
   """

   #TODO all of these defaults set to objects are bad practice because the defaults will be instanciated at parse time.
   #should
   fields= { "force" :   (InputForce,    { "help"  : InputForce.default_help }),
             "ensemble": (InputEnsemble, { "help"  : InputEnsemble.default_help } ),
             "prng" :    (InputRandom,   { "help"  : InputRandom.default_help + " Optional.",
                                         "default" : input_default(factory=Random)} ),
             "initialize" : (InputInitializer, { "help" : InputInitializer.default_help,
                                                "default" : input_default(factory=Initializer) } ),
             "beads" :   (InputBeads, { "help"     : InputBeads.default_help,
                                        "default"  : input_default(factory=Beads, kwargs={'natoms': 0, 'nbeads': 0}) } ),
             "normal_modes" :   (InputNormalModes, { "help"     : InputNormalModes.default_help,
                                        "default"  : input_default(factory=NormalModes, kwargs={'mode': "rpmd"}) } ),
             "cell" :    (InputCell,   { "help"    : InputCell.default_help }),
             "output" :  (InputOutputs, { "help"   : InputOutputs.default_help,
                                          "default": input_default(factory=InputOutputs.make_default)  }),
             "step" :       ( InputValue, { "dtype"    : int,
                                            "default"  : 0,
                                            "help"     : "How many time steps have been done." }),
             "total_steps": ( InputValue, { "dtype"    : int,
                                            "default"  : 1000,
                                            "help"     : "The total number of steps that will be done." })
                                             }

   default_help = "This is the top level class that deals with the running of the simulation, including holding the simulation specific properties such as the time step and outputting the data."
   default_label = "SIMULATION"

   def store(self, simul):
      """Takes a simulation instance and stores a minimal representation of it.

      Args:
         simul: A simulation object.
      """

      super(InputSimulation,self).store()
      self.force.store(simul._forcemodel)
      self.ensemble.store(simul.ensemble)

      self.beads.store(simul.beads)

      self.normal_modes.store(simul.nm)
      self.cell.store(simul.cell)
      self.prng.store(simul.prng)
      self.step.store(simul.step)
      self.total_steps.store(simul.tsteps)
      self.output.store(simul.outputs)

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

      nbeads = self.beads.fetch()
      ncell = self.cell.fetch()
      nprng = self.prng.fetch()


      # this creates a simulation object which gathers all the little bits
      rsim = engine.simulation.Simulation(nbeads, ncell, self.force.fetch(),
                     self.ensemble.fetch(), nprng, self.output.fetch(),
                     self.normal_modes.fetch(), self.initialize.fetch(), self.step.fetch(),
                     tsteps=self.total_steps.fetch())

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

      #~ if self.beads._explicit:
         #~ # nothing to be done here! user/restart provides a beads object
         #~ pass
      #~ elif self.atoms._explicit:
         #~ # user is providing atoms: assume a classical simulation
         #~ atoms = self.atoms.fetch()
         #~ nbeads = 1
         #~ rbeads = Beads(atoms.natoms, nbeads)
         #~ rbeads[0] = atoms.copy()
         #~ # we create a dummy beads storage so that fetch can proceed as if a
         #~ # beads object had been specified
         #~ self.beads.store(rbeads)
      #~ else:
         #~ raise TypeError("Either a <beads> or a <atoms> block must be provided")

      if self.total_steps.fetch() <= self.step.fetch():
         raise ValueError("Current step greater than total steps, no dynamics will be done.")

      #~ for init in self.initialize.fetch():
         #~ if not init in ["velocities", "normal_modes", "cell_velocities"]:
            #~ raise ValueError("Initialization parameter " + init + " is not a valid keyword for initialize.")

