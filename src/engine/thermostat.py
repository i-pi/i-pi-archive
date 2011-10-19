import numpy, math
from utils.depend import *
class thermostat(object): 
   p=None; atoms=None; cell=None; 
   
   def __init__(self, temp = 1.0, dt = 1.0, econs=0.0):
      self.temp=depend(name='temp', value=temp)
      self.dt=depend(name='dt', value=dt)
      self.econs = depend(name='econs',value=econs)

   def bind(self, atoms, p, cell): pass
      
   def step(self):                 pass

   def cell_step(self):            pass
