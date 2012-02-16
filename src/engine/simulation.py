"""Contains the class that deals with the running of the simulation and 
outputting the results.

The root class for the whole simulation. Contains references to all the top
level objects used in the simulation, and controls all the steps that are
not inherently system dependent, like the running of each time step, 
choosing which properties to initialise, and which properties to output.

Classes:
   RestartSimulation: Deals with creating the simulation object from a file, 
      and writing the checkpoints.
   Simulation: Deals with running the simulation and outputting the results.
"""

__all__ = ['RestartSimulation', 'Simulation']

import numpy as np
import math, random
from utils.depend import *
from utils.restart import *
from utils.units  import *
from utils.prng   import *
from utils.io     import *
from atoms import *
from cell import *
from ensembles import RestartEnsemble
from forces import RestartForce, ForceBeads
from beads import Beads, RestartBeads
from properties import Properties

_DEFAULT_STRIDES={"checkpoint": 1000, "properties": 10, "progress": 100, "trajectory": 20,  "trajectory_full": 100}
_DEFAULT_OUTPUT=[ "time", "conserved", "kinetic", "potential" ]
class RestartSimulation(Restart):
   """Simulation restart class.

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
      nbeads: A float giving the number of beads.
      step: An integer giving the current simulation step. Defaults to 0.
      total_steps: The total number of steps. Defaults to 0.
      stride: A dictionary giving the number of steps between printing out 
         data for the different types of data. Defaults to _DEFAULT_STRIDES.
      prefix: A string giving the prefix for all the output files. Defaults to
         'prefix'.
      traj_format: A string giving the format of the trajectory output files. 
         Defaults to 'pdb'.
      properties: An array of strings giving all the properties that should 
         be output space separated. Defaults to _DEFAULT_OUTPUT.
      initialize: An array of strings giving all the quantities that should
         be output.
      fd_delta: A float giving the size of the finite difference
         parameter used in the Yamamoto kinetic energy estimator. Defaults 
         to 0.
   """

   fields= { "force" : (RestartForce, ()),  "ensemble": (RestartEnsemble, ()), 
             "atoms" : (RestartAtoms, ()), "beads" : (RestartBeads, ()), 
             "cell" : (RestartCell, ()), "prng" : (RestartRandom, ()),
             "nbeads": (RestartValue, (int, 0 ) ),              
             "step" : ( RestartValue, (int, 0)), 
             "total_steps": (RestartValue, (int, 1000) ), 
             "stride" : ( RestartValue, (dict, {})),
             "prefix": (RestartValue, (str, "prefix")), 
             "traj_format": (RestartValue, (str, "pdb")),              
             "properties": (RestartArray, (str,np.zeros(0, np.dtype('|S12'))) ),
             "initialize": (RestartArray, (str,np.zeros(0, np.dtype('|S12'))) ),
             "fd_delta":   ( RestartValue, (float, 0.0)) 
            }

   def store(self, simul):
      """Takes a simulation instance and stores a minimal representation of it.

      Args:
         simul: A simulation object.
      """

      self.force.store(simul._forcemodel)
      self.ensemble.store(simul.ensemble)
      self.beads.store(simul.beads)
      self.cell.store(simul.cell)
      self.prng.store(simul.prng)
      self.step.store(simul.step)
      self.total_steps.store(simul.tsteps)
      self.stride.store(simul.dstride)
      self.prefix.store(simul.prefix)
      self.traj_format.store(simul.format)      
      self.properties.store(simul.outlist)
      self.initialize.store(simul.initlist)
      self.fd_delta.store(simul.properties.fd_delta)
            
   def fetch(self):
      """Creates a simulation object.

      Returns:
         A simulation object of the appropriate type and with the appropriate
         properties and other objects given the attributes of the 
         RestartSimulation object.
      """

      self.check()
      nbeads = self.beads.fetch()
      ncell = self.cell.fetch()
      nprng = self.prng.fetch()

      dstride = dict(_DEFAULT_STRIDES)
      istride = self.stride.fetch()
      vstride = {}
      for k,s in istride.items(): 
         if not k in dstride:
            raise TypeError(k + " is not a valid input for the stride keyword")
         vstride[k] = int(s)
      dstride.update(vstride)
      
      olist = self.properties.fetch()
      if (len(olist) == 0):
         olist = None

      ilist = self.initialize.fetch()
      if (len(ilist) == 0):
         ilist = None
      
      rsim = Simulation(nbeads, ncell, self.force.fetch(), 
                     self.ensemble.fetch(), nprng, self.step.fetch(), 
                     tsteps=self.total_steps.fetch(), stride=dstride,
                     prefix=self.prefix.fetch(), format=self.traj_format.fetch(), outlist=olist, initlist=ilist)

      if (self.fd_delta.fetch() != 0.0):
         rsim.properties.fd_delta = self.fd_delta.fetch()
      return rsim

   def check(self):
      """Function that deals with optional arguments.

      Deals with the difference between classical and PI dynamics. If there is
      no beads argument, the bead positions are generated from the atoms, with 
      the necklace being fixed at the atom position. Similarly, if no nbeads
      argument is specified a classical simulation is done.
      """

      if self.beads.nbeads.fetch() == 0:
         atoms = self.atoms.fetch() 
         if atoms.natoms == 0:
            raise TypeError("Either a <beads> or a <atoms> block must be provided")
         nbeads = self.nbeads.fetch()
         if nbeads == 0:
            nbeads = 1
         rbeads = Beads(atoms.natoms, nbeads)
         for b in range(rbeads.nbeads):
            rbeads[b] = atoms.copy()
         self.beads.store(rbeads)      

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
      initlist: An array of the properties that should be initialised. Set to 
         zero after the initialisation, so that the checkpoints don't specify
         any properties to be initialised.
      properties: A properties object.
      fout: File to output the properties to.
      tcout: File to output the centroid trajectory to.
      tout: File to output the full trajectory to.
      ichk: A number keeping track of all the restart files generated so far,
         so that old files are not overwritten.

   Depend objects:
      step: The current simulation step.
   """   

   def __init__(self, beads, cell, force, ensemble, prng, step=0, tsteps=1000, 
               stride=None, prefix="prefix", format="pdb", outlist=None, initlist=None):
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
         format: An optional string giving the output format for the 
            trajectory files. Defaults to 'pdb'.
         outlist: An array of strings giving all the properties that should 
            be output space separated.
         initlist: An array of strings giving all the quantities that should
            be output.
      """

      print " # Initializing simulation object "
      self.nbeads = len(beads)
      self.beads = beads
      self.cell = cell
      self.prng = prng
      self._forcemodel = force
      self.forces = ForceBeads()
            
      self.ensemble = ensemble
      dset(self, "step", depend_value(name="step", value=step))
      self.tsteps = tsteps
      
      self.prefix = prefix
      self.format = format      
      if stride is None:
         self.dstride = dict(_DEFAULT_STRIDES)
      else:
         self.dstride = stride
      if outlist is None:
         self.outlist = np.array(_DEFAULT_OUTPUT, np.dtype('|S12') )
      else:
         self.outlist = outlist                                    

      self.properties = Properties()
         
      if initlist is None:
         self.initlist = np.zeros(0, np.dtype('|S12'))
      else:
         self.initlist = initlist
                  
      self.bind()
      
   def bind(self):
      """Calls the bind routines for all the objects in the simulation."""

      self.forces.bind(self.beads, self.cell,  self._forcemodel)
      self.ensemble.bind(self.beads, self.cell, self.forces, self.prng)
      self.properties.bind(self)
      self.init()

      # Checks as soon as possible if some asked-for properties are missing or mispelled
      for what in self.outlist:
         if not what in self.properties.property_dict.keys():
            print "Computable properties list: ", self.properties.property_dict.keys()
            raise KeyError(what + " is not a recognized property")


   def run(self):
      """Runs the simulation.

      Does all the simulation steps, and outputs data to the appropriate files
      when necessary. Also deals with starting and cleaning up the threads used
      in the communication between the driver and the PIMD code.
      """

      self._forcemodel.socket.start_thread()
   
      if (self.step == 0):
         self.step = -1
         self.write_output()
         if self.format == "pdb": 
            io_pdb.print_pdb_path(self.beads, self.cell, self.tout)
            io_pdb.print_pdb(self.beads.centroid, self.cell, self.tcout)
         elif self.format == "xyz": 
            io_xyz.print_xyz_path(self.beads, self.cell, self.tout)            
            io_xyz.print_xyz(self.beads.centroid, self.cell, self.tcout)
         self.step = 0
      for self.step in range(self.step,self.tsteps):               
         self.ensemble.step()
         print " # MD step % 7d complete. Timings --> p-step: %10.5f  q-step: %10.5f  t-step: %10.5f" % (self.step, self.ensemble.ptime, self.ensemble.qtime, self.ensemble.ttime )
         print " # MD diagnostics: V: %10.5e    Kcv: %10.5e   Ecns: %10.5e" % (self.properties["potential"], self.properties["kinetic_cv"], self.properties["conserved"] )

         if ((self.step+1) % self.dstride["checkpoint"] == 0):
            self.write_chk()      
         if ((self.step+1) % self.dstride["properties"] == 0):
            self.write_output()
         if ((self.step+1) % self.dstride["trajectory_full"] == 0):            
            if self.format == "pdb":
               io_pdb.print_pdb_path(self.beads, self.cell, self.tout)
            elif self.format == "xyz":
               io_xyz.print_xyz_path(self.beads, self.cell, self.tout)
         if ((self.step+1) % self.dstride["trajectory"] == 0):
            if self.format == "pdb":
               io_pdb.print_pdb(self.beads.centroid, self.cell, self.tcout)
            elif self.format == "xyz":
               io_xyz.print_xyz(self.beads.centroid, self.cell, self.tcout)

      if self.step < self.tsteps:
         self.step += 1         
      self.write_chk()

      self._forcemodel.socket.end_thread()      
   
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
         ohead += "%16s"%(l)
      self.fout.write(ohead + "\n")
      self.tcout = open(self.prefix + "." + self.format, "a")
      self.tout = open(self.prefix + ".full." + self.format, "a")      
      self.ichk = 0
      
      if "velocities" in self.initlist:
         self.beads.p = math.sqrt(self.ensemble.ntemp*Constants.kb)*self.beads.sm3*self.prng.gvec((self.beads.nbeads, 3*self.beads.natoms))

      if self.ensemble.fixcom:
         self.ensemble.rmcom()
      
      # Zeroes out the initlist, such that in restarts no initialization will be required
      self.initlist = np.zeros(0, np.dtype('|S12'))
   
   def write_output(self):
      """Outputs the required properties of the system.

      Note that properties are outputted using the same format as for the 
      output to the xml checkpoint files, as specified in io_xml.

      Raises:
         KeyError: Raised if one of the properties specified in the output list
            are not contained in the property_dict member of properties.
      """

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
      
      r = RestartSimulation()
      r.store(self)
      check_file.write(r.write(name="simulation"))
      check_file.close()
      
