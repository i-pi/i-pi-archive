import numpy
import math
import thermostat

class Thermo_test(thermostat.Thermostat):     
   
   @property
   def dt(self):
      print "Thermo_test getter called"
      return self.dt_Base
     
   @dt.setter
   def dt(self,new):
      print "Thermo_test setter called"
      self.dt_Base = new
   
   def __init__(self):
      thermostat.Thermostat.__init__(self)
      self.__tau=1.0
      self.__T=1.0
      self.__S=1.0

