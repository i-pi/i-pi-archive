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
#from properties import *
#from forces import *

class RestartSimulation(Restart):
   fields= { "force" : (RestartForce, ()),  "ensemble": (RestartEnsemble, ()), 
             "atoms" : (RestartAtoms, ()), "beads" : (RestartBeads, ()), 
             "cell" : (RestartCell, ()), "prng" : (RestartRandom, ()),
             "init_velocities" : (RestartValue, (bool, False)), 
             "init_cell" : (RestartValue, (bool, False)), 
             "nbeads": (RestartValue, (int, 0 ) ),              
             "step" : ( RestartValue, (int, 0)), 
             "total_steps": (RestartValue, (int, 1000) ) }

   def store(self, simul):
      self.force.store(simul.force)
      self.ensemble.store(simul.ensemble)
      self.beads.store(simul.beads)
      self.cell.store(simul.cell)
      self.prng.store(simul.prng)
      self.step.store(simul.step)
      self.total_steps.store(simul.tsteps)
            
   def fetch(self):
      self.check()
      nbeads=self.beads.fetch();  ncell= self.cell.fetch(); nprng=self.prng.fetch()
      if self.init_velocities.fetch():
         nbeads.p=np.sqrt(units.Constants.kb*self.ensemble.fetch().temp*nbeads.m3)*nprng.gvec((nbeads.nbeads,nbeads.natoms*3))
      if self.init_cell.fetch():
         ncell.p6=np.sqrt(units.Constants.kb*self.ensemble.fetch().temp*ncell.m6)*nprng.gvec(6)
      return Simulation(nbeads, ncell, self.force.fetch(), self.ensemble.fetch(), 
                     nprng, self.step.fetch(), tsteps=self.total_steps.fetch() )

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
               step=0, tsteps=1000):
      self.nbeads=len(beads)
      self.beads=beads
      self.cell=cell
      self.prng=prng
      self._forcemodel=force
      self.forces=ForceBeads()
         
      self.ensemble=ensemble
      self.step=step
      self.tsteps=tsteps
      self.bind()
      
   def bind(self):
      self.forces.bind(self.beads, self.cell,  self._forcemodel)
      self.ensemble.bind(self.beads, self.cell, self.forces, self.prng)

   def run(self):

      fout = open("output","w")
      tout = open("pitraj.pdb","w")      
      for self.step in range(self.step,self.tsteps):
         self.ensemble.step()
         fout.write( str(self.step)+" "+str(self.ensemble.econs)+" "+str(self.beads.kin.sum())+" "+
                  str(self.forces.pot.sum())+" "+str(self.beads.vpath*self.ensemble.omegan2)+" "+"\n");
         fout.flush()
         print "times:  #p", self.ensemble.ptime/(self.step+1), "  #q",  self.ensemble.qtime/(self.step+1), "  #t",  self.ensemble.ttime/(self.step+1)
         print  "objvalue count: ", objgraph.count('depend_value'), objgraph.count('depend_array'), objgraph.count('Atom'),

         if self.step % 20 ==0:
            io_pdb.print_pdb(self.beads, self.cell, tout)
            print "Cleaning up", gc.collect()
      self.step+=1
#         
#class RestartSimulation(Restart):
#   fields= { "force" : (RestartForce, ()), "atoms" : (RestartAtoms, ()), "cell" : (RestartCell, ()),
#             "ensemble": (RestartEnsemble, ()), "steps_done" : ( RestartValue, (int, 0)), 
#             "total_steps": (RestartValue, (int, 1000)), "output": (RestartOutput, ()) }

#   def store(self, simul):
#      self.force.store(simul.force)
#      self.atoms.store(simul.atoms)
#      self.cell.store(simul.cell)
#      self.ensemble.store(simul.ensemble)
#      self.steps_done.store(simul.step)
#      self.total_steps.store(simul.tsteps)
#      self.output.store(simul.output)
#            
#   def fetch(self):
#      return Simulation(self.atoms.fetch(), self.cell.fetch(), self.force.fetch(), self.ensemble.fetch(), 
#                        self.output.fetch(), step = self.steps_done.fetch(), tsteps=self.total_steps.fetch())
#   
#class Simulation(dobject):
#   """Represents a simulation cell. Includes the cell parameters, 
#      the atoms and the like."""   

#   def __init__(self, atoms, cell, force, ensemble, output, step=0, tsteps=1000):
#      self.atoms=atoms
#      self.cell=cell
#      self.force=force
#      self.ensemble=ensemble
#      self.output=output
#      self.step=step
#      self.tsteps=tsteps
#      self.properties = Properties()
#      self.bind()
#      
#   def bind(self):
#      self.force.bind(self.atoms, self.cell)
#      self.ensemble.bind(self.atoms, self.cell, self.force)
#      self.properties.bind(self.ensemble, self.atoms, self.cell, self.force)
#      self.output.bind(self.properties, self)

#   def run(self):      
#      for self.step in range(self.step,self.tsteps):
#         self.ensemble.step()
#        # print str(self.step+1)+" "+str(self.ensemble.econs)+" "+str(self.atoms.kin)+" "+str(self.force.pot)+" "+str(self.ensemble.thermostat.ethermo)+" "+"\n",
#         self.output.step(self.step)
#      self.step+=1
#         
