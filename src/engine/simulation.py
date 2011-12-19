import numpy as np
import math, random
from utils.depend import *
from utils.restart import Restart
from utils.units  import *
from utils.io     import *
from atoms import *
from cell import *
from ensembles import RestartEnsemble
from forces import RestartForce
from properties import *
#from forces import *

class RestartSimulation(Restart):
   fields= { "force" : (RestartForce, ()), "atoms" : (RestartAtoms, ()), "cell" : (RestartCell, ()),
             "ensemble": (RestartEnsemble, ()), "steps_done" : ( RestartValue, (int, 0)), 
             "total_steps": (RestartValue, (int, 1000) ) }

   def store(self, simul):
      self.force.store(simul.force)
      self.atoms.store(simul.atoms)
      self.cell.store(simul.cell)
      self.ensemble.store(simul.ensemble)
      self.steps_done.store(simul.step)
      self.total_steps.store(simul.tsteps)
            
   def fetch(self):
      return Simulation(self.atoms.fetch(), self.cell.fetch(), self.force.fetch(), self.ensemble.fetch(), 
                     self.steps_done.fetch(), tsteps=self.total_steps.fetch() )
   
class Simulation(dobject):
   """Represents a simulation cell. Includes the cell parameters, 
      the atoms and the like."""   

   def __init__(self, atoms, cell, force, ensemble, step=0, tsteps=1000):
      self.atoms=atoms
      self.cell=cell
      self.force=force
      self.ensemble=ensemble
      self.step=step
      self.tsteps=tsteps
      self.bind()
      
   def bind(self):
      self.force.bind(self.atoms, self.cell)
      self.ensemble.bind(self.atoms, self.cell, self.force)

   def run(self):      
      for self.step in range(self.step,self.tsteps):
         self.ensemble.step()
         print str(self.step+1)+" "+str(self.ensemble.econs)+" "+str(self.atoms.kin)+" "+str(self.force.pot)+" "+str(self.ensemble.thermostat.ethermo)+" "+"\n",
      self.step+=1
         
