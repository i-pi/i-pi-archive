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
from utils.restart import Restart
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
      prng: A restart random number instance.
      nbeads: A float giving the number of beads.
      step: An integer giving the current simulation step. Defaults to 0.
      total_steps: The total number of steps. Defaults to 0.
      stride: A dictionary giving the number of steps between printing out 
         data for the different types of data. Defaults to _DEFAULT_STRIDES.
      prefix: A string giving the prefix for all the output files. Defaults to
         'prefix'.
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
                     prefix=self.prefix.fetch(), outlist=olist, initlist=ilist)

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
      beads: A beads object giving the atoms positions.
      cell: A cell object giving the system box.
      prng: A random number object.
      _forcemodel: A forcefield object giving the force calculator for each 
         replica of the system.
      forces: A ForceBeads object for calculating the forces for all the 
         replicas.
      ensemble: An ensemble object giving the objects necessary for producing
         the correct ensemble.
      tsteps: The total number of steps.
      prefix: A prefix for all the output files.
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
               stride=None, prefix="prefix", outlist=None, initlist=None):
      """Initialises Simulation class.

      Args:
      """

      print "@@@@@ initializing with total_Steps", tsteps, "and dt", ensemble.dt
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
      if stride is None:
         self.dstride = dict(_DEFAULT_STRIDES)
      else:
         self.dstride = stride
      if outlist is None:
         self.outlist = np.array(_DEFAULT_OUTPUT,np.dtype('|S12') )
      else:
         self.outlist = outlist
      if initlist is None:
         self.initlist = np.zeros(0, np.dtype('|S12'))
      else:
         self.initlist = initlist

      self.properties = Properties()
      
      self.bind()
      
   def bind(self):
      self.forces.bind(self.beads, self.cell,  self._forcemodel)
      self.ensemble.bind(self.beads, self.cell, self.forces, self.prng)
      self.properties.bind(self)
      self.init()

   def run(self):
      self._forcemodel.socket.start_thread()
   
      if (self.step == 0):
         self.step = -1
         self.write_output()
         io_pdb.print_pdb_path(self.beads, self.cell, self.tout)
         io_pdb.print_pdb(self.beads.centroid, self.cell, self.tcout)
         self.step = 0
      for self.step in range(self.step,self.tsteps):               
         self.ensemble.step()
         if ((self.step+1) % self.dstride["checkpoint"] == 0):
            self.write_chk()      
         if ((self.step+1) % self.dstride["properties"] == 0):
            self.write_output()
         if ((self.step+1) % self.dstride["trajectory_full"] == 0):
            io_pdb.print_pdb_path(self.beads, self.cell, self.tout)
         if ((self.step+1) % self.dstride["trajectory"] == 0):
            io_pdb.print_pdb(self.beads.centroid, self.cell, self.tcout)

         print "times:  #p", self.ensemble.ptime/(self.step+1), "  #q",  self.ensemble.qtime/(self.step+1), "  #t",  self.ensemble.ttime/(self.step+1)

      if self.step < self.tsteps:
         self.step += 1         
      self.write_chk()

      self._forcemodel.socket.end_thread()      
   
   def init(self):
      self.fout = open(self.prefix + ".out", "a")
      ohead="# "
      for l in self.outlist:
         ohead += "%16s"%(l)
      self.fout.write(ohead + "\n")      
      self.tcout = open(self.prefix + ".pdb", "a")
      self.tout = open(self.prefix + ".full.pdb", "a")      
      self.ichk = 0
      
      if "velocities" in self.initlist:
         self.beads.p = math.sqrt(self.ensemble.ntemp*Constants.kb)*self.beads.sm3*self.prng.gvec((self.beads.nbeads, 3*self.beads.natoms))
      if self.ensemble.fixcom:
         self.ensemble.rmcom()
         #self.ensemble.rmcom()
      
      self.initlist = np.zeros(0, np.dtype('|S12'))
   
   def write_output(self):
      for what in self.outlist:
         try:
            quantity = self.properties[what]
         except KeyError:
            raise KeyError(what + " is not a recognized property")
         self.fout.write(write_type(float, quantity) + " ")
         
      self.fout.write("\n")   
      self.fout.flush()   
      
   def write_chk(self):
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
      
