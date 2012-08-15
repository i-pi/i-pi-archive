"""Contains the class that deals with the running of the simulation and 
outputting the results.

The root class for the whole simulation. Contains references to all the top
level objects used in the simulation, and controls all the steps that are
not inherently system dependent, like the running of each time step, 
choosing which properties to initialise, and which properties to output.

Classes:
   Simulation: Deals with running the simulation and outputting the results.
"""

__all__ = ['Simulation', 'PropertyOutput' ]

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

class PropertyOutput(dobject):
   """ Mini-class dealing with outputting a set of properties to file.
   
      Does not do any calculation, just manages opening a file, getting data 
      from a Properties object and outputting with the desired stride. 
   """
   
   def __init__(self, filename="out", stride=1, outlist=[]):
      """ Initializes a property output stream opening the corresponding file name. 
         Also writes out headers.  """
      self.filename=filename
      self.outlist=outlist
      self.stride=stride
      self.out=None
      
      
   def bind(self, simul):
      """ Binds output proxy to simulation object. """

      self.simul=simul
      
      # Checks as soon as possible if some asked-for properties are missing or mispelled
      for what in self.outlist:
         if '(' in what:
            argstart = what.find('(')
            key = what[0:argstart]
            if not key in self.simul.properties.property_dict.keys():
               print "Computable properties list: ", self.properties.property_dict.keys()
               raise KeyError(key + " is not a recognized property")

         elif not what in self.simul.properties.property_dict.keys():
            print "Computable properties list: ", self.properties.property_dict.keys()
            raise KeyError(what + " is not a recognized property")
         
      self.open_stream()
      
   def open_stream(self):
      """ Opens the output stream. """

      try:
         self.out=open(self.filename, "a")
      except:
         raise ValueError("Could not open file "+self.filename+" for output")

      ohead = "# "
      for l in self.outlist:
         ohead += "%16s"%(l) + " "
      self.out.write(ohead + "\n")
         

   def close_stream():
      """ Closes the output stream. """ 
      
      self.out.close()
            
   def write(self):
      """Outputs the required properties of the system. 

      Note that properties are outputted using the same format as for the 
      output to the xml checkpoint files, as specified in io_xml.

      Raises:
         KeyError: Raised if one of the properties specified in the output list
            are not contained in the property_dict member of properties.
      """

   
      if not (self.simul.step+1) % self.stride == 0:   return
      self.out.write("  ")      
      for what in self.outlist:
         try:
            quantity = self.simul.properties[what]
         except KeyError:
            raise KeyError(what + " is not a recognized property")
         if not hasattr(quantity,"__len__") :
            self.out.write(write_type(float, quantity) + "   ")
         else:
            self.out.write(" [ ")
            for el in quantity:
               self.fout.write(write_type(float, el) + " ")
            self.out.write(" ] ")
                     
      self.out.write("\n")   
      self.out.flush()   


class TrajectoryOutput(dobject):
   """ Mini-class dealing with outputting atom-based properties as a trajectory file.

      Does not do any calculation, just manages opening a file, getting data 
      from a Trajectories object and outputting with the desired stride. 
   """
   
   
   def __init__(self, filename="out", stride=1, what="", format="xyz"):
      """ Initializes a property output stream opening the corresponding file name. 
         Also writes out headers.  """
      self.filename=filename
      self.what=what
      self.stride=stride
      self.format=format
      self.out=None
      
      
   def bind(self, simul):
      """ Binds output proxy to simulation object. """

      self.simul=simul      
      self.open_stream()
      
   def open_stream(self):
      """ Opens the output stream(s). """

      if self.what in [ "positions", "velocities", "forces" ]: 
         # must write out trajectories for each bead, so must create b streams
         self.out=[]
         for b in range(self.simul.beads.nbeads):
            padb=( ("%0"+str(int(1+np.floor(np.log(self.simul.beads.nbeads)/np.log(10))))+"d") % (b) )
            try:
               self.out.append( open(self.filename+"_"+padb+"."+self.format, "a") )
            except:
               raise ValueError("Could not open file "+self.filename+"_"+padb+"."+self.format+" for output")         

   def close_stream():
      """ Closes the output stream. """ 
      
      if hasattr(self.out, "__getitem__"):   
         for o in self.out: o.close()
      else:
         self.out.close()

   def write(self):
      """Writes out the required trajectories."""
      
      if not (self.simul.step+1) % self.stride == 0:   return
      
      # quick-and-dirty way to check if a trajectory is "global" or per-bead
      # Checks to see if there is a list of files or just a single file.
      if hasattr(self.out, "__getitem__"):   
         for b in range(len(self.out)):
            self.simul.trajs.print_traj(self.what, self.out[b], b, format=self.format)
      else:
         self.simul.trajs.print_traj(self.what, self.out, b=0, format=self.format)

class CheckpointOutput(dobject):
   """ Mini-class dealing with outputting checkpoints. """
   
   
   def __init__(self, filename="restart", stride=1000, overwrite=True, step=0):
      """ Initializes a property output stream opening the corresponding file name. 
         Also writes out headers.  """
      self.filename=filename
      self.step=step
      self.stride=stride
      self.overwrite=overwrite
      
      
   def bind(self, simul):
      """ Binds output proxy to simulation object. """

      self.simul=simul
      self.status=inputs.simulation.InputSimulation(); self.status.store(simul)
      
   def store(self):
      self.status.store(self.simul)
      
   def write(self):
      """Writes out the required trajectories."""
      
      if not (self.simul.step+1) % self.stride == 0:   return
      
      if self.overwrite: filename=self.filename
      else: filename=self.filename+"_"+str(self.step)      
      
      self.store()
      check_file = open(filename, "w")            
      check_file.write(self.status.write(name="simulation"))
      check_file.close()      
      
      self.step+=1
         
                 
                 

class Simulation(dobject):
   """Main simulation object. 

   Contains all the references and the main dynamics loop. Also handles the
   initialisation and output.
   
   Attributes:
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
      format: A string specifying both the format and the extension of traj 
         output.
      outputs: A list of output objects that should 
      initlist: A dictionary of the properties that should be initialised with
         their values. Set to zero after the initialisation, so that the 
         checkpoints don't specify any properties to be initialised after the 
         simulation is restarted.
      properties: A properties object.
      chk: A CheckoutOutput object used to deal gracefully with soft exit
         in the middle of a force calculation

   Depend objects:
      step: The current simulation step.
   """   

   def __init__(self, beads, cell, force, ensemble, prng, outputs, step=0, tsteps=1000,  initlist=None):
      """Initialises Simulation class.

      Args:
         beads: A beads object giving the atom positions.
         cell: A cell object giving the system box.
         force: A forcefield object giving the force calculator for each replica
            of the system.
         ensemble: An ensemble object giving the objects necessary for 
            producing the correct ensemble.
         prng: A random number object.
         outputs: A list of output objects. 
         step: An optional integer giving the current simulation time step. 
            Defaults to 0.
         tsteps: An optional integer giving the total number of steps. Defaults
            to 1000.
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
      self.outputs = outputs
      self.chk = None
      
      self.ensemble = ensemble
      
      dset(self, "step", depend_value(name="step", value=step))
      self.tsteps = tsteps                             

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
      for o in self.outputs:    o.bind(self)

      self.chk=CheckpointOutput("RESTART", 1, True, 0)
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
      self.chk.write()

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
         for o in self.outputs:  o.write()
         self.step = 0
                 
      # main MD loop
      for self.step in range(self.step,self.tsteps):   
         # stores the state before doing a step. 
         # this is a bit time-consuming but makes sure that we can honor soft 
         # exit requests without screwing the trajectory
         self.chk.store() 
         
         self.ensemble.step()
         print " # MD step % 7d complete. Timings --> p-step: %10.5f  q-step: %10.5f  t-step: %10.5f" % (self.step, self.ensemble.ptime, self.ensemble.qtime, self.ensemble.ttime )
         print " # MD diagnostics: V: %10.5e    Kcv: %10.5e   Ecns: %10.5e" % (self.properties["potential"], self.properties["kinetic_cv"], self.properties["conserved"] )
            
         for o in self.outputs:  o.write()

         if os.path.exists("EXIT"): # soft-exit
            self.soft_exit(rollback=False)

      self.soft_exit(rollback=False)
            
   def init(self):
      """Deals with initialization.

      Opens the different output files. Also initialises the 
      atom velocities, and the higher frequency normal modes if required.
      It then removes the list of quantities to
      be initialized, so that if the simulation is restarted these quantities
      are not re-initialized.
      """  

      if "normal_modes" in self.initlist:
         init_temp = float(self.initlist["normal_modes"])*self.beads.nbeads
         for b in range(1,self.beads.nbeads):
            if (self.beads.qnm[b] == 0.0).all:
               if init_temp == 0: 
                  self.beads.qnm[b] = math.sqrt(self.ensemble.ntemp*Constants.kb)/(self.ensemble.omegak[b]*self.beads.sm3[b])*np.prng.gvec(3*self.beads.natoms)
                  self.beads.pnm[b] = math.sqrt(self.ensemble.ntemp*Constants.kb)*self.beads.sm3[b]*self.prng.gvec(3*self.beads.natoms)
               else:
                  self.beads.qnm[b] = math.sqrt(init_temp*Constants.kb)/(self.ensemble.omegak[b]*self.beads.sm3[b])*self.prng.gvec(3*self.beads.natoms)
                  self.beads.pnm[b] = math.sqrt(init_temp*Constants.kb)*self.beads.sm3[b]*self.prng.gvec(3*self.beads.natoms)

      if "velocities" in self.initlist:
         init_temp = float(self.initlist["velocities"])*self.beads.nbeads
         if init_temp == 0:
            self.beads.p = math.sqrt(self.ensemble.ntemp*Constants.kb)*self.beads.sm3*self.prng.gvec((self.beads.nbeads, 3*self.beads.natoms))
         else:
            self.beads.p = math.sqrt(init_temp*Constants.kb)*self.beads.sm3*self.prng.gvec((self.beads.nbeads, 3*self.beads.natoms))

      if "cell_velocities" in self.initlist:
         init_temp = float(self.initlist["cell_velocities"])*self.beads.nbeads
         if init_temp == 0:
            init_temp = math.sqrt(self.ensemble.ntemp)
         if hasattr(self.cell,"p6"):
            self.cell.p6 = math.sqrt(init_temp*Constants.kb*self.cell.m)*self.prng.gvec(6)
         else:
            self.cell.P = math.sqrt(init_temp*Constants.kb*self.cell.m)*self.prng.gvec(1)

      if self.ensemble.fixcom:
         self.ensemble.rmcom()
      
      # Zeroes out the initlist, such that in restarts no initialization will be required
      self.initlist = {}
