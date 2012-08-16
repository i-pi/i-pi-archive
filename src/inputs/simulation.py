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
from inputs.atoms import InputAtoms
from inputs.beads import InputBeads
from inputs.cell import InputCell
from inputs.ensembles import InputEnsemble
from inputs.outputs import InputOutputs
from engine.atoms import Atoms
from engine.beads import Beads
from engine.outputs import *
import engine.simulation


_DEFAULT_STRIDES = {"checkpoint": 1000, "properties": 10, "progress": 100, "centroid": 20,  "trajectory": 100}
_DEFAULT_OUTPUT = [ "time", "conserved", "kinetic_cv", "potential" ]
_DEFAULT_TRAJ = [ "positions" ]

_DefaultOutput=None #InputOutputs()
#~ _DefaultOutput.store([
                     #~ engine.simulation.PropertyOutput("md", 10, [ "time", "step", "conserved", "temperature", "potential", "kinetic_cv" ]),
                     #~ engine.simulation.TrajectoryOutput("pos", 100, "positions", "xyz"),
                     #~ engine.simulation.CheckpointOutput("checkpoint",1000,overwrite=True)
                    #~ ])
#~ _DefaultOutput.store("wrap_pi")

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

   fields= { "force" :   (InputForce,    { "help"  : InputForce.default_help }),
             "ensemble": (InputEnsemble, { "help"  : InputEnsemble.default_help } ),
             "prng" :    (InputRandom,   { "help"  : InputRandom.default_help + " It is not necessary to specify this tag.",
                                         "default" : Random() } ),
             "atoms" :   (InputAtoms, { "help"     : "Deals with classical simulations. Only needs to be specified if a classical simulation is required, and should be left blank otherwise.",
                                        "default"  : Atoms(0) } ),
             "beads" :   (InputBeads, { "help"     : InputBeads.default_help + " Only needs to be specified if the atoms tag is not, but overwrites it otherwise.",
                                        "default"  : Beads(0,1) } ),
             "cell" :    (InputCell,   { "help"    : InputCell.default_help }),
             "output" :  (InputOutputs, { "help" : InputOutputs.default_help ,
                                          "default" : _DefaultOutput  }),
             "step" :       ( InputValue, { "dtype"    : int,
                                            "default"  : 0,
                                            "help"     : "How many time steps have been done." }),
             "total_steps": ( InputValue, { "dtype"    : int,
                                            "default"  : 1000,
                                            "help"     : "The total number of steps that will be done." }),
             "initialize":  ( InputValue, { "dtype"    : dict,
                                            "default"  : {},
                                            "help"     : "A dictionary giving the properties of the system that need to be initialized, and their initial values. The allowed keywords are ['velocities', 'cell_velocities', 'normal_modes']. The initial value of 'velocities' corresponds to the temperature to initialise the velocity distribution from. If 0, then the system temperature is used. 'cell_velocities' is the same but for the cell velocity. The initial value of 'normal_modes' corresponds to the temperature from which to initialize the higher normal mode frequencies from, if we start a simulation from a configuration with a smaller number of beads. If 0, then the system temperature is used." }),
             "fd_delta":    ( InputValue, { "dtype"    : float,
                                            "default"  : 0.0,
                                            "help"     : "The parameter used in the finite difference differentiation in the calculation of the scaled path velocity estimator. Defaults to 1e-5." })
#             "traj_format": ( InputValue, { "dtype"    : str,
#                                            "default"  : "pdb",
#                                            "help"     : "The file format for the output file. Allowed keywords are ['pdb', 'xyz'].",
#                                            "options"  : ["pdb", "xyz"] }),
#
#             "trajectories": ( InputArray, { "dtype"   : str,
#                                             "default" : np.zeros(0, np.dtype('|S12')),
#                                             "help"    :
#                                             "A list of the properties to print out the per-atom or per-bead trajectories of. Allowed values are ['positions', 'velocities', 'forces', 'kinetic_cv', 'kodterms_cv', 'momentum_centroid', 'centroid', 'gyration', 'spring']. 'kinetic_cv' gives the quantum kinetic energy estimator for each degree of freedom, whereas 'kodterms_cv' gives the off-diagonal elements of the kinetic stress tensor estimator for each degree of freedom. 'gyration' gives the radius of gyration of each atom. 'spring' prints the per-DOF spring term in the PI Hamiltonian. The others are self-explanatory."})
#
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

      # If we are running a classical simulation, hide the "beads" machinery in the restarts
      if simul.beads.nbeads > 1 :
         self.beads.store(simul.beads)
      else:
         self.atoms.store(simul.beads[0])

      self.cell.store(simul.cell)
      self.prng.store(simul.prng)
      self.step.store(simul.step)
      self.total_steps.store(simul.tsteps)
      self.output.store(simul.outputs)
      self.initialize.store(simul.initlist)
      self.fd_delta.store(simul.properties.fd_delta)

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

      ilist = self.initialize.fetch()
      if not self.initialize._explicit:
         ilist = None

      rsim = engine.simulation.Simulation(nbeads, ncell, self.force.fetch(),
                     self.ensemble.fetch(), nprng, self.output.fetch(), self.step.fetch(),
                     tsteps=self.total_steps.fetch(),  initlist=ilist)

      if self.fd_delta._explicit:
         rsim.properties.fd_delta = self.fd_delta.fetch()

      # binds and inits the simulation object just before returning
      rsim.bind()
      rsim.init()

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

      if self.beads._explicit:
         # nothing to be done here! user/restart provides a beads object
         pass
      elif self.atoms._explicit:
         # user is providing atoms: assume a classical simulation
         atoms = self.atoms.fetch()
         nbeads = 1
         rbeads = Beads(atoms.natoms, nbeads)
         rbeads[0] = atoms.copy()
         # we create a dummy beads storage so that fetch can proceed as if a
         # beads object had been specified
         self.beads.store(rbeads)
      else:
         raise TypeError("Either a <beads> or a <atoms> block must be provided")

      if self.total_steps.fetch() <= self.step.fetch():
         raise ValueError("Current step greater than total steps, no dynamics will be done.")

      for init in self.initialize.fetch():
         if not init in ["velocities", "normal_modes", "cell_velocities"]:
            raise ValueError("Initialization parameter " + init + " is not a valid keyword for initialize.")

      #TODO do something about fd_delta...
