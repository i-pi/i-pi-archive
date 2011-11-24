import numpy, math
from utils.depend import *

class barostat(): 

   def __init__(self, pext = numpy.zeros((3,3)), dt = 1.0):
      self.pext=depend(name='pext', value=pext)
      self.dt=depend(name='dt', value=dt)

   def bind(self, syst):
      """Dummy binding function"""

      pass
      
   def pstep(self):                
      """Dummy atoms thermostat step""" 
      
      pass

   def rstep(self):
      pass

