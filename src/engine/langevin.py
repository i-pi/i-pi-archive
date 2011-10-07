import numpy
import math
import thermostat

class Thermo_Langevin(thermostat.Thermostat):     
   
   def compute_TS(self):
      print "Re-computing propagator"
      self.__T=math.exp(-self.__dt/self.__tau)
      self.__S=math.sqrt(self.__temp*(1-self.__T**2))
    
   @property
   def dt(self):
      return self.__dt
     
   @dt.setter
   def dt(self,new):
      print "Thermo_Langevin setter called"
      self.__dt = new
      self.compute_TS()
   
   def __init__(self):
      thermostat.Thermostat.__init__(self)
      self.__temp=self.temp
      self.__tau=1.0
      self.__T=1.0
      self.__S=1.0
      self.__dt = self.dt_Base

