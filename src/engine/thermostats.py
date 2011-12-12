import numpy as np
import math
from utils.depend import *
from utils.units  import *

class Thermostat(dobject): 
   """Represent a thermostat for constant T simulations.
      Contains: temp = temperature, dt = time step, econs =
      change in the kinetic energy due to the thermostat
      Initialised by: thermo = thermostat(temp, dt, econs)
      temp = temperature, default = 1.0
      dt = time step, default = 1.0
      econs = conserved energy quantity, default = 0.0"""

   def __init__(self, temp = 1.0, dt = 1.0, ethermo=0.0):
      dset(self,"temp",   depend_value(name='temp', value=temp))
      dset(self,"dt",     depend_value(name='dt', value=dt))      
      dset(self,"ethermo",depend_value(name='ethermo',value=ethermo))

   def bind(self, atoms=None, cell=None, pm=None):
      """Dummy binding function"""     
      
      if not atoms is None:
         dset(self,"p",dget(atoms, "p"))
         dset(self,"m",dget(atoms, "m3"))      
      elif not cell is None:   
         dset(self,"p",dget(cell, "p6"))
         dset(self,"m",dget(cell, "m6"))      
      elif not pm is None:   
         dset(self,"p",pm[0])
         dset(self,"m",pm[1])               
      else: 
         raise TypeError("Thermostat.bind expects either an Atoms, a Cell, or a (p,m) tuple to bind to")
      
      dset(self,"sqrtm",depend_array(name="sqrtm", 
                                     value=np.zeros(len(dget(self,"m"))), 
                                     deps=depend_func(func=self.get_sqrtm, dependencies=[depget(self,"m")] ) ) )      
      
   def get_sqrtm(self):
      return np.sqrt(self.m)
   
   def step(self):                
      """Dummy atoms thermostat step"""       
      pass

class ThermoLangevin(Thermostat):     
   """Represent a langevin thermostat for constant T simulations.
      Contains: temp = temperature, dt = time step, econs =
      change in the kinetic energy due to the thermostat,
      tau = thermostat mass, sm = sqrt(mass), (T,S) = thermostat parameters
      Initialised by: thermo = langevin(temp, dt, tau, econs)
      temp = temperature, default = 1.0
      dt = time step, default = 1.0
      tau = thermostat mass, default = 1.0
      econs = conserved energy quantity, default = 0.0"""

   def get_T(self):
      """Calculates T in p(0) = T*p(dt) + S*random.gauss()"""
      return math.exp(-self.dt/self.tau)
      
   def get_S(self):      
      """Calculates S in p(0) = T*p(dt) + S*random.gauss()"""
      return math.sqrt(Constants.kb*self.temp*(1-self.T**2))   
  
   def __init__(self, temp = 1.0, dt = 1.0, tau = 1.0, ethermo=0.0):
      super(ThermoLangevin,self).__init__(temp,dt,ethermo)
      
      dset(self,"tau",depend_value(value=tau,name='tau'))
      dset(self,"T",depend_value(name="T",deps=depend_func(func=self.get_T, dependencies=[depget(self,"tau"), depget(self,"dt")])))
      dset(self,"S",depend_value(name="S",deps=depend_func(func=self.get_S, dependencies=[depget(self,"temp"), depget(self,"T")])))
      
   def step(self):
      """Updates the atom velocities with a langevin thermostat"""
      
      p=numpy.array(self.p,copy=True)
      
      p/=self.sqrtm
      self.ethermo+=numpy.dot(p,p)*0.5
      p*=self.T
      p+=self.S*np.random.standard_normal(len(p))
      self.ethermo-=numpy.dot(p,p)*0.5      
      p*=self.sqrtm
      
      self.p=p

