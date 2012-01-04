import numpy as np
import math, random
import gc, objgraph
from utils.depend import *
from utils.restart import Restart
from utils.units  import *
from utils.prng   import *
from utils.io     import *
from atoms import *
from cell import *
from ensembles import RestartEnsemble
from forces import RestartForce, ForceBeads
#from output import RestartOutput
from beads import Beads, RestartBeads
from properties import Properties
#from forces import *

_DEFAULT_STRIDES={"checkpoint": 1000, "properties": 10, "progress": 100, "trajectory": 20,  "trajectory_full": 100}
_DEFAULT_OUTPUT=[ "time", "conserved", "kinetic", "potential" ]
class RestartSimulation(Restart):
   fields= { "force" : (RestartForce, ()),  "ensemble": (RestartEnsemble, ()), 
             "atoms" : (RestartAtoms, ()), "beads" : (RestartBeads, ()), 
             "cell" : (RestartCell, ()), "prng" : (RestartRandom, ()),
             "nbeads": (RestartValue, (int, 0 ) ),              
             "step" : ( RestartValue, (int, 0)), 
             "total_steps": (RestartValue, (int, 1000) ), 
             "stride" : ( RestartValue, (dict, {})),
             "prefix": (RestartValue, (str, "prefix")), 
             "properties": (RestartArray, (str,np.zeros(0, np.dtype('|S12'))) ),
             "initialize": (RestartArray, (str,np.zeros(0, np.dtype('|S12'))) )             
            }

   def store(self, simul):
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
            
   def fetch(self):
      self.check()
      nbeads=self.beads.fetch();  ncell= self.cell.fetch(); nprng=self.prng.fetch()
#      if self.init_velocities.fetch():
#         nbeads.p=np.sqrt(units.Constants.kb*self.ensemble.fetch().temp*nbeads.m3)*nprng.gvec((nbeads.nbeads,nbeads.natoms*3))
#      if self.init_cell.fetch():
#         ncell.p6=np.sqrt(units.Constants.kb*self.ensemble.fetch().temp*ncell.m6)*nprng.gvec(6)
      
      # gets a list of strides
      dstride=dict(_DEFAULT_STRIDES)
      istride=self.stride.fetch()
      vstride={}
      for k,s in istride.items(): 
         if not k in dstride: raise TypeError(k+" is not a valid input for the stride keyword")
         vstride[k]=int(s)
      dstride.update(vstride)
      
      olist=self.properties.fetch()
      if (len(olist)==0) : olist=None

      ilist=self.initialize.fetch()
      if (len(ilist)==0) : ilist=None
      
      # list of properties to be printed out
      return Simulation(nbeads, ncell, self.force.fetch(), self.ensemble.fetch(), 
                     nprng, self.step.fetch(), tsteps=self.total_steps.fetch(), stride=dstride,
                     prefix=self.prefix.fetch(), outlist=olist, initlist=ilist )

   def check(self):
      if self.beads.nbeads.fetch() == 0:  # no beads keyword provided, default to creating a necklace centered at a single configuration -- provided by the atoms keyword
         atoms=self.atoms.fetch() # fetches the template atoms block and creates a beads object
         if atoms.natoms == 0:  raise TypeError("Either a <beads> or a <atoms> block must be provided")
         nbeads=self.nbeads.fetch()
         if nbeads==0: nbeads=1   # defaults to classical simulation
         rbeads=Beads(atoms.natoms, nbeads)
         for b in range(rbeads.nbeads): rbeads[b]=atoms.copy()
         self.beads.store(rbeads)
      
class Simulation(dobject):
   """Represents a simulation cell. Includes the cell parameters, 
      the atoms and the like."""   

   def __init__(self, beads, cell, force, ensemble, prng,
               step=0, tsteps=1000, stride=None, prefix="prefix", outlist=None, initlist=None):
      print "@@@@@ initializing with total_Steps", tsteps, "and dt", ensemble.dt
      self.nbeads=len(beads)
      self.beads=beads
      self.cell=cell
      self.prng=prng
      self._forcemodel=force
      self.forces=ForceBeads()
            
      self.ensemble=ensemble
      dset(self,"step",depend_value(name="step",value=step) )
      self.tsteps=tsteps
      
      # output & properties
      self.prefix=prefix
      if stride is None: self.dstride=dict(_DEFAULT_STRIDES)
      else:              self.dstride=stride
      if outlist is None: self.outlist=np.array(_DEFAULT_OUTPUT,np.dtype('|S12') )
      else:               self.outlist=outlist
      if initlist is None: self.initlist=np.zeros(0, np.dtype('|S12'))
      else:               self.initlist=initlist

      
      self.properties=Properties()
      
      self.bind()
      
   def bind(self):
      self.forces.bind(self.beads, self.cell,  self._forcemodel)
      self.ensemble.bind(self.beads, self.cell, self.forces, self.prng)
      self.properties.bind(self)
      self.init()

   def run(self):
      self._forcemodel.socket.start_thread()
   
      if (self.step == 0):
         self.write_output();  io_pdb.print_pdb_path(self.beads, self.cell, self.tout);  io_pdb.print_pdb(self.beads.centroid, self.cell, self.tcout)
      for self.step in range(self.step,self.tsteps):               
         self.ensemble.step()
         if ((self.step+1) % self.dstride["checkpoint"] ==0) : self.write_chk()      
         if ((self.step+1) % self.dstride["properties"] ==0) : self.write_output()
         if ((self.step+1) % self.dstride["trajectory_full"] ==0) : io_pdb.print_pdb_path(self.beads, self.cell, self.tout)
         if ((self.step+1) % self.dstride["trajectory"] ==0) : io_pdb.print_pdb(self.beads.centroid, self.cell, self.tcout)         

         print "times:  #p", self.ensemble.ptime/(self.step+1), "  #q",  self.ensemble.qtime/(self.step+1), "  #t",  self.ensemble.ttime/(self.step+1)
         print  "objvalue count: ", objgraph.count('depend_value'), objgraph.count('depend_array'), objgraph.count('Atom'),

      if self.step < self.tsteps: self.step+=1         
      self.write_chk()

      self._forcemodel.socket.end_thread()      
   
   def init(self):
      # initialize output
      self.fout=open(self.prefix+".out", "a")
      ohead="# "
      for l in self.outlist: ohead+="%16s"%(l)
      self.fout.write(ohead+"\n")      
      self.tcout=open(self.prefix+".pdb", "a")
      self.tout=open(self.prefix+".full.pdb", "a")      
      self.ichk=0
      
      # initialize velocities
      if "velocities" in self.initlist:  self.beads.p=math.sqrt(self.ensemble.ntemp*Constants.kb)* self.beads.sm3*self.prng.gvec((self.beads.nbeads, 3*self.beads.natoms))
      
      self.initlist=np.zeros(0, np.dtype('|S12'))   # sets this to nothing so in the checkpoint initialization won't be required
   
   def write_output(self):
      for what in self.outlist:
         try:
            quantity = self.properties[what]
         except: raise TypeError(what+" is not a recognized property")
         self.fout.write(write_type(float, quantity) + " ")
         
      self.fout.write("\n")   
      self.fout.flush()   
      
   def write_chk(self):
      # tries to open a new restart file
      new = False;  self.ichk += 1
      while (not new):
         try:
            check_file = open(self.prefix + ".restart" + str(self.ichk), "r")
            check_file.close()
            self.ichk += 1
         except IOError:
            check_file = open(self.prefix + ".restart" + str(self.ichk), "w")
            new = True
      
      r=RestartSimulation()
      r.store(self)
      check_file.write(r.write(name="simulation"))
      check_file.close()
      
