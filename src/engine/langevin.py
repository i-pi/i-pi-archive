import numpy
import math
import thermostat

class Thermo_Langevin(thermostat.Thermostat):     
   
   def compute_TS(self):
      print "Re-computing propagator"
      self.__T=math.exp(-self.dt_Base/self.__tau)
      self.__S=math.sqrt(self.temp*(1-self.__T**2))
    
   @property
   def dt(self):
      print "langevin getter called"
      return self.dt_Base
     
   @dt.setter
   def dt(self,new):
      print "Thermo_Langevin setter called"
      self.dt_Base = new
      self.compute_TS()
   
   def __init__(self):
      thermostat.Thermostat.__init__(self)
      self.__tau=1.0
      self.__T=1.0
      self.__S=1.0

