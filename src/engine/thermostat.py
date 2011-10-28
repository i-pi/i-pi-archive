import numpy, math
from utils.depend import *
class thermostat(object): 
   """Represent a thermostat for constant T simulations.
      Contains: temp = temperature, dt = time step, econs =
      change in the kinetic energy due to the thermostat
      Initialised by: thermo = thermostat(temp, dt, econs)
      temp = temperature, default = 1.0
      dt = time step, default = 1.0
      econs = conserved energy quantity, default = 0.0"""

   p=None; atoms=None; cell=None; 
   
   def __init__(self, temp = 1.0, dt = 1.0, econs=0.0):
      self.temp=depend(name='temp', value=temp)
      self.dt=depend(name='dt', value=dt)
      self.econs = depend(name='econs',value=econs)

   def bind(self, atoms, p, cell):
      """Dummy binding function"""

      pass
      
   def step(self):                
      """Dummy atoms thermostat step""" 
      
      pass

   def NST_cell_step(self):
      """Dummy NST cell thermostat step"""   
      
      pass

   def NPT_cell_step(self):
      """Dummy NPT cell thermostat step"""   
      
      pass
