import numpy
import math
import thermostat

class Thermo_Langevin(thermostat.Thermostat):     
   
   def compute_TS(self):
      print "Re-computing propagator"
      self.T=math.exp(-self._dt/self._tau)
      self.S=math.sqrt(self._temp*(1-self._T**2))
         
   @thermostat.Thermostat.dt.setter
   def dt(self,new):
      self._dt = new
      self.compute_TS()
   
   def __init__(self):
      thermostat.Thermostat.__init__(self)
      print self._temp
      self._tau=1.0
      self._T=1.0
      self._S=1.0

