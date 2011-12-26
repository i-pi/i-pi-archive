import numpy as np
import math
from utils.depend   import *
from utils.units    import *
from utils.restart  import *
from utils.prng import Random


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

   def bind(self, atoms=None, cell=None, pm=None, prng=Random()):
      self.prng=prng  
      
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
      
      dset(self,"sm",depend_array(name="sm", value=np.zeros(len(dget(self,"m"))), 
                                     func=self.get_sm, dependencies=[dget(self,"m")] ) )
      
   def get_sm(self):  return np.sqrt(self.m)
   
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
      return math.exp(-0.5*self.dt/self.tau)
      
   def get_S(self):      
      """Calculates S in p(0) = T*p(dt) + S*random.gauss()"""
      return math.sqrt(Constants.kb*self.temp*(1-self.T**2))   
  
   def __init__(self, temp = 1.0, dt = 1.0, tau = 1.0, ethermo=0.0):
      super(ThermoLangevin,self).__init__(temp,dt,ethermo)
      
      dset(self,"tau",depend_value(value=tau,name='tau'))
      dset(self,"T",depend_value(name="T",func=self.get_T, dependencies=[dget(self,"tau"), dget(self,"dt")]))
      dset(self,"S",depend_value(name="S",func=self.get_S, dependencies=[dget(self,"temp"), dget(self,"T")]))
      
   def step(self):
      """Updates the atom velocities with a langevin thermostat"""
      
      p=self.p.view(np.ndarray).copy()
      
      p/=self.sm

      self.ethermo+=np.dot(p,p)*0.5
      p*=self.T
      p+=self.S*self.prng.gvec(len(p))
      self.ethermo-=np.dot(p,p)*0.5      

      p*=self.sm      
            
      self.p=p
      
      
class RestartThermo(Restart):
   attribs={ "kind": (RestartValue, (str, "langevin")) }
   fields={ "ethermo" : (RestartValue, (float, 0.0)), 
            "tau" : (RestartValue, (float, 1.0))  }
   
   def store(self, thermo):
      if type(thermo) is ThermoLangevin: 
         self.kind.store("langevin")
         self.tau.store(thermo.tau)
      else: self.kind.store("unknown")      
      self.ethermo.store(thermo.ethermo)
      
   def fetch(self):
      if self.kind.fetch() == "langevin" : thermo=ThermoLangevin(tau=self.tau.fetch())
      else: thermo=Thermostat()
      thermo.ethermo=self.ethermo.fetch()
      return thermo      
     
      
