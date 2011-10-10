import numpy
import math
import thermostat
import random

class Thermo_Langevin(thermostat.Thermostat):     
   
   def compute_TS(self):
      print "Re-computing propagator"
      self.__T=math.exp(-self.__dt/self.__tau)
      self.__S=math.sqrt(self.temp*(1-self.__T**2))
   
   @property
   def temp(self):
      print "langevin temp getter called"
      return self.__temp

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
      self.__tau=1.0
      self.k_Boltz = 1.0
      self.__temp = temp
      self.dt = dt

   def step(self, atom):
      sigma = 1.0/(4*math.pi*self.__tau*self.k_Boltz*self.temp*atom.mass)
      for i in range(3):
         atom.p[i] = self.__T*atom.p[i] + self.__S*random.gauss(0.0, sigma)

