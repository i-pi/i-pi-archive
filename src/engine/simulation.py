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
import math, random
import os.path, sys
from utils.depend import *
from utils.units  import *
from utils.prng   import *
from utils.io     import *
from utils.io.io_xml import *
from atoms import *
import time
from cell import *
from forces import ForceBeads
from beads import Beads
from properties import Properties, Trajectories
import inputs.simulation

_DEFAULT_STRIDES = {"checkpoint": 1000, "properties": 10, "progress": 100, "centroid": 20,  "trajectory": 100}
_DEFAULT_OUTPUT = [ "time", "conserved", "kinetic_cv", "potential" ]
_DEFAULT_TRAJ = [ "positions" ]

class Simulation(dobject):
   """Main simulation object. 

   Contains all the references and the main dynamics loop. Also handles the
   initialisation and output.
   
   Attributes:
      nbeads: The number of the replicas of the system.
      beads: A beads object giving the atom positions.
      cell: A cell object giving the system box.
      prng: A random number generator object.
      _forcemodel: A forcefield object giving the force calculator for each 
         replica of the system.
      forces: A ForceBeads object for calculating the forces for all the 
         replicas.
      ensemble: An ensemble object giving the objects necessary for producing
         the correct ensemble.
      tsteps: The total number of steps.
      prefix: A prefix for all the output files.
      format: A string specifying both the format and the extension of traj 
         output.
      dstride: A dictionary giving number of steps between printing out 
         data for the different types of data. Defaults to _DEFAULT_STRIDES.
      outlist: An array of strings giving the different properties to output.
      initlist: A dictionary of the properties that should be initialised with
         their values. Set to zero after the initialisation, so that the 
         checkpoints don't specify any properties to be initialised after the 
         simulation is restarted.
      properties: A properties object.
      fout: File to output the properties to.
      tout: File to output the full trajectory to.
      ichk: A number keeping track of all the restart files generated so far,
         so that old files are not overwritten.
      status: A InputSimulation object used to deal gracefully with soft exit
         in the middle of a force calculation

   Depend objects:
      step: The current simulation step.
   """   

   def __init__(self, beads, cell, force, ensemble, prng, step=0, tsteps=1000, 
               stride=None, prefix="prefix", outlist=None, trajlist=None, 
               initlist=None):
      """Initialises Simulation class.

      Args:
         beads: A beads object giving the atom positions.
         cell: A cell object giving the system box.
         force: A forcefield object giving the force calculator for each replica
            of the system.
         ensemble: An ensemble object giving the objects necessary for 
            producing the correct ensemble.
         prng: A random number object.
         step: An optional integer giving the current simulation time step. 
            Defaults to 0.
         tsteps: An optional integer giving the total number of steps. Defaults
            to 1000.
         stride: An optional dictionary giving the number of steps between 
            printing out data for the different types of data. 
         prefix: An optional string giving the prefix for all the output files. 
            Defaults to 'prefix'.
         outlist: An array of strings giving all the properties that should 
            be output space separated.
         trajlist: An array of strings giving all the trajectories that should
            be output. 
         initlist: A dictionary of keys giving all the quantities that should
            be initialized with values giving their initial value.
      """

      print " # Initializing simulation object "
      self.nbeads = len(beads)
      self.beads = beads
      self.cell = cell
      self.prng = prng
      self._forcemodel = force
      self.forces = ForceBeads()
      self.status = None
      
      self.ensemble = ensemble
      
      dset(self, "step", depend_value(name="step", value=step))
      self.tsteps = tsteps
      
      self.prefix = prefix      
      if stride is None:
         self.dstride = dict(_DEFAULT_STRIDES)
      else:
         self.dstride = stride

      if outlist is None:
         self.outlist = np.array(_DEFAULT_OUTPUT, np.dtype('|S12') )
      else:
         self.outlist = outlist                                    

      if trajlist is None:
         self.trajlist = np.array(_DEFAULT_TRAJ, np.dtype('|S12') )
      else:
         self.trajlist = trajlist                                    

      self.properties = Properties()
      self.trajs = Trajectories()
         
      if initlist is None:
         self.initlist = {}
      else:
         self.initlist = initlist

                        
   def bind(self):
      """Calls the bind routines for all the objects in the simulation.

      Raises:
         KeyError: Raised if one of the requested output properties is not 
            defined in the property list.
      """

      # binds important computation engines
      self.forces.bind(self.beads, self.cell,  self._forcemodel, softexit=self.soft_exit)
      self.ensemble.bind(self.beads, self.cell, self.forces, self.prng)

      # binds output management objects
      self.properties.bind(self)
      self.trajs.bind(self)
      
      self.status = inputs.simulation.InputSimulation()
      self.status.store(self)
      
      # Checks as soon as possible if some asked-for properties are missing or mispelled
      for what in self.outlist:
         if '(' in what:
            argstart = what.find('(')
            key = what[0:argstart]
            if not key in self.properties.property_dict.keys():
               print "Computable properties list: ", self.properties.property_dict.keys()
               raise KeyError(key + " is not a recognized property")

         elif not what in self.properties.property_dict.keys():
            print "Computable properties list: ", self.properties.property_dict.keys()
            raise KeyError(what + " is not a recognized property")
   
   def soft_exit(self, rollback=True):
      """Deals with a soft exit request. 
      
      Tries to ensure that a consistent restart checkpoint is 
      written out. 
      """   
      
      if self.step < self.tsteps:
         self.step += 1         
      if not rollback:
         self.status.store(self)         
      self.write_chk()

      self._forcemodel.socket.end_thread()      
      sys.exit()

   def run(self):
      """Runs the simulation.

      Does all the simulation steps, and outputs data to the appropriate files
      when necessary. Also deals with starting and cleaning up the threads used
      in the communication between the driver and the PIMD code.
      """

      self._forcemodel.socket.start_thread()
   
      # prints inital configuration -- only if we are not restarting
      if (self.step == 0):
         self.step = -1
         self.write_output()
         self.write_traj()
         self.step = 0
                 
      # main MD loop
      for self.step in range(self.step,self.tsteps):   
         # stores the state before doing a step. 
         # this is a bit time-consuming but makes sure that we can honor soft 
         # exit requests without screwing the trajectory
         self.status.store(self) 
         
         self.ensemble.step()
         print " # MD step % 7d complete. Timings --> p-step: %10.5f  q-step: %10.5f  t-step: %10.5f" % (self.step, self.ensemble.ptime, self.ensemble.qtime, self.ensemble.ttime )
         print " # MD diagnostics: V: %10.5e    Kcv: %10.5e   Ecns: %10.5e" % (self.properties["potential"], self.properties["kinetic_cv"], self.properties["conserved"] )

         if ((self.step+1) % self.dstride["checkpoint"] == 0):
            self.write_chk()      
         if ((self.step+1) % self.dstride["properties"] == 0):
            self.write_output()
         if ((self.step+1) % self.dstride["trajectory"] == 0):
            self.write_traj()
         if os.path.exists("EXIT"): # soft-exit
            self.soft_exit(rollback=False)
            
      self.soft_exit(rollback=False)
            
   def init(self):
      """Deals with the file initialization.

      Opens the different output files. Also initialises the cell and
      atom velocities if required, and then removes the list of quantities to
      be initialized, so that if the simulation is restarted these quantities
      are not re-initialized.
      """

      self.fout = open(self.prefix + ".out", "a")
      ohead = "# "
      for l in self.outlist:
         ohead += "%16s"%(l) + " "
      self.fout.write(ohead + "\n")
      
      
      # self.tcout = open(self.prefix + ".centroid." + self.trajs.format, "a")  
      self.tout = {}    
      for what in self.trajlist:
         if what in [ "positions", "velocities", "forces" ]:
            self.tout[what] = [ open(self.prefix + "." + what[0:3] + "_" + str(b) + "." + self.trajs.format, "a") for b in range(self.beads.nbeads) ]
         else:
            self.tout[what] = open(self.prefix + "." + what[0:3] + "." + self.trajs.format, "a")

      self.ichk = 0      
      if "velocities" in self.initlist:
         init_temp = float(self.initlist["velocities"])*self.beads.nbeads
         if init_temp == 0:
            self.beads.p = math.sqrt(self.ensemble.ntemp*Constants.kb)*self.beads.sm3*self.prng.gvec((self.beads.nbeads, 3*self.beads.natoms))
         else:
            self.beads.p = math.sqrt(init_temp*Constants.kb)*self.beads.sm3*self.prng.gvec((self.beads.nbeads, 3*self.beads.natoms))

      if self.ensemble.fixcom:
         self.ensemble.rmcom()
      
      # Zeroes out the initlist, such that in restarts no initialization will be required
      self.initlist = {}

   def write_traj(self):
      """Writes out the required trajectories."""
      
      for what in self.trajlist:
         # quick-and-dirty way to check whether a trajectory is "global" or per-bead
         if hasattr(self.tout[what], "__getitem__"):   
            for b in range(self.beads.nbeads):
               self.trajs.print_traj(what, self.tout[what][b], b)
         else:
            self.trajs.print_traj(what, self.tout[what])
   
   def write_output(self):
      """Outputs the required properties of the system.

      Note that properties are outputted using the same format as for the 
      output to the xml checkpoint files, as specified in io_xml.

      Raises:
         KeyError: Raised if one of the properties specified in the output list
            are not contained in the property_dict member of properties.
      """

      self.fout.write("  ")
      for what in self.outlist:
         try:
            quantity = self.properties[what]
         except KeyError:
            raise KeyError(what + " is not a recognized property")
         self.fout.write(write_type(float, quantity) + " ")
         
      self.fout.write("\n")   
      self.fout.flush()   
      
   def write_chk(self):
      """Outputs the xml checkpoint files.

      Note that each checkpoint is written to a different file, so that old 
      checkpoint files can still be used to restart the system if a problem
      is found.
      """

      new = False
      self.ichk += 1
      while (not new):
         try:
            check_file = open(self.prefix + ".restart" + str(self.ichk), "r")
            check_file.close()
            self.ichk += 1
         except IOError:
            check_file = open(self.prefix + ".restart" + str(self.ichk), "w")
            new = True
      
      check_file.write(self.status.write(name="simulation"))
      check_file.close()
      
