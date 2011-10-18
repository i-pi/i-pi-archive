import numpy, math
from utils.depend import *
class Thermostat(dobject): 
   def __init__(self, temp = 1.0, dt = 1.0, econs=0.0):
      self.temp=depend(name='temp', value=temp)
      self.dt=depend(name='dt', value=dt)
      
      self.econs = econs

   def step(self):
      pass

   def cell_step(self):
      pass
