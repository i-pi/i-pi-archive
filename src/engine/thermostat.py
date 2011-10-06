import numpy
import math

class Thermostat(object): 
   @property
   def temp(self):
      return self._temp

   @property
   def dt(self):
      return self._dt

   @dt.setter
   def dt(self,new):
      print "Thermo dt setter called"
      self._dt = new
     
   def __init__(self):
      self._temp=1.0
      self._dt=1.0
