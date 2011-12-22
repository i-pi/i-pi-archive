import numpy as np
import math, random
from utils.depend import *
from utils.restart import Restart
from utils.units  import *
from utils.io     import *
from atoms import *
from cell import *
#from ensembles import RestartEnsemble
from forces import RestartForce, ForceBeads
#from output import RestartOutput
from pi_beads import Beads, RestartBeads
#from properties import *
#from forces import *

class RestartSimulation(Restart):
   fields= { "force" : (RestartForce, ()),  #"ensemble": (RestartEnsemble, ()), 
             "atoms" : (RestartAtoms, ()), "beads" : (RestartBeads, ()), 
             "cell" : (RestartCell, ()), "nbeads": (RestartValue, (int, 0 ) ),              
             "step" : ( RestartValue, (int, 0)), 
             "total_steps": (RestartValue, (int, 1000) ) }

   def store(self, simul):
      self.force.store(simul.force)
#      self.ensemble.store(simul.ensemble)

      self.beads.store(simul.beads)
      self.cell.store(simul.cell)
      self.step.store(simul.step)
      self.total_steps.store(simul.tsteps)
            
   def fetch(self):
      self.check()
      return Simulation(self.beads.fetch(), self.cell.fetch(), self.force.fetch(), #self.ensemble.fetch(), 
                     self.step.fetch(), tsteps=self.total_steps.fetch() )

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

   def __init__(self, beads, cell, force, #ensemble, 
               step=0, tsteps=1000):
      self.nbeads=len(beads)
      self.beads=beads
      self.cell=cell
      self._forcemodel=force
      self.forces=ForceBeads()
         
#      self.ensemble=ensemble
      self.step=step
      self.tsteps=tsteps
      self.bind()
      
   def bind(self):
      pass
#      for b in range(self.beads.nbeads): 
#         self.force[b].bind(self.beads[b], self.cell)
      self.forces.bind(self.beads, self.cell,  self._forcemodel)
#      self.ensemble.bind(self.atoms, self.cell, self.force)

   def run(self):
      for self.step in range(self.step,self.tsteps): pass
#         self.ensemble.step()
#         print str(self.step)+" "+str(self.ensemble.econs)+" "+str(self.atoms.kin)+" "+str(self.force.pot)+" "+str(self.ensemble.thermostat.ethermo)+" "+"\n",
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
