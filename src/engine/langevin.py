import numpy
import math
import thermostat

class Thermo_Langevin(thermostat.Thermostat):     
   
   def compute_TS(self):
      print "Re-computing propagator"
      self.__T=math.exp(-self.__dt/self.__tau)
      self.__S=math.sqrt(self.temp*(1-self.__T**2))
   
   #@property
   #def temp(self):
   #   print "langevin temp getter called"
   #   return self.__temp

   @property
   def dt(self):
      print "langevin getter called"
      return self.__dt
     
   @dt.setter
   def dt(self,new):
      print "Thermo_Langevin setter called"
      self.__dt = new
      self.compute_TS()
   
   def __init__(self, temp = 1.0, dt = 1.0):
      self.__dt = dt
      self.__tau=1.0
      self.__T=1.0
      self.__S=1.0
      self.temp = temp

