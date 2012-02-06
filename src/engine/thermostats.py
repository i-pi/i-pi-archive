import numpy as np
import math
from utils.depend   import *
from utils.units    import *
from utils.restart  import *
from utils.prng import Random
from beads import Beads

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

   def bind(self, beads=None, atoms=None, cell=None, pm=None, prng=None, ndof=None):

      if prng is None:  self.prng=Random()
      else: self.prng=prng  
      
      if not beads is None:
         dset(self,"p",beads.p.flatten())
         dset(self,"m",beads.m3.flatten())     
      elif not atoms is None:
         dset(self,"p",dget(atoms, "p"))
         dset(self,"m",dget(atoms, "m3"))               
      elif not cell is None:   
         dset(self,"p",dget(cell, "p6"))
         dset(self,"m",dget(cell, "m6"))      
      elif not pm is None:   
         dset(self,"p",pm[0])
         dset(self,"m",pm[1])               
      else: 
         raise TypeError("Thermostat.bind expects either Beads, Atoms, a Cell, or a (p,m) tuple to bind to")
      
      if ndof is None: self.ndof=len(self.p)
      else: self.ndof=ndof
      
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
      tau = thermostat time scale, sm = sqrt(mass), (T,S) = thermostat parameters
      Initialised by: thermo = langevin(temp, dt, tau, econs)
      temp = temperature, default = 1.0
      dt = time step, default = 1.0
      tau = thermostat relaxation time, default = 1.0
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
      dset(self,"T",  depend_value(name="T",func=self.get_T, dependencies=[dget(self,"tau"), dget(self,"dt")]))
      dset(self,"S",  depend_value(name="S",func=self.get_S, dependencies=[dget(self,"temp"), dget(self,"T")]))
      
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

class ThermoPILE_L(Thermostat):    
   def __init__(self, temp = 1.0, dt = 1.0, tau = 1.0, ethermo=0.0):
      super(ThermoPILE_L,self).__init__(temp,dt,ethermo)
      dset(self,"tau",depend_value(value=tau,name='tau'))

   # optionally does not bind the centroid, so we can re-use all of this in the PILE_G case
   def bind(self, beads=None, prng=None, bindcentroid=True, ndof=None):  
      if beads is None or not type(beads) is Beads: raise TypeError("ThermoPILE_L.bind expects a Beads argument to bind to")
      if prng is None:  self.prng=Random()
      else: self.prng=prng  
      
      # creates a set of thermostats to be applied to individual normal modes
      self._thermos=[ ThermoLangevin(temp=1, dt=1, tau=1) for b in range(beads.nbeads) ]
      if not bindcentroid: self._thermos[0]=None
      
      dset(self,"tauk",depend_array(name="tauk", value=np.zeros(beads.nbeads-1,float), func=self.get_tauk, dependencies=[dget(self,"temp")]) )
      
      # must pipe all the dependencies in such a way that values for the nm thermostats
      # are automatically updated based on the "master" thermostat
      def make_taugetter(k): return lambda: self.tauk[k-1]
      it=0
      for t in self._thermos:
         if t is None: 
            it+=1; continue   
         t.bind(pm=(beads.pnm[it,:],beads.m3[0,:]),prng=self.prng, ndof=(ndof if it==0 else None)) # bind thermostat t to the it-th bead
         # pipes temp and dt
         deppipe(self,"temp", t, "temp")
         deppipe(self,"dt", t, "dt")

         # for tau it is slightly more complex
         if it==0:
            deppipe(self,"tau", t, "tau")
         else:
            # here we manually connect _thermos[i].tau to tauk[i]. simple and clear.
            dget(t,"tau").add_dependency(dget(self,"tauk"))
            dget(t,"tau")._func=make_taugetter(it)
         dget(self,"ethermo").add_dependency(dget(t,"ethermo"))
         it+=1     
               
      dget(self,"ethermo")._func=self.get_ethermo;
         
   def get_tauk(self):  
      return np.array([1.0/(4*self.temp*Constants.kb/Constants.hbar*math.sin(k*math.pi/len(self._thermos))) for k in range(1,len(self._thermos)) ])

   def get_ethermo(self):
      et=0.0;
      for t in self._thermos: et+=t.ethermo
      return et
            
   def step(self):
      # super-cool! just loop over the thermostats! it's as easy as that! 
      for t in self._thermos: t.step()        

class ThermoSVR(Thermostat):     
   """Represent a stochastic velocity rescaling thermostat for constant T simulations.
      Contains: temp = temperature, dt = time step, econs =
      change in the kinetic energy due to the thermostat,
      tau = thermostat relaxation time, sm = sqrt(mass), (T,S) = thermostat parameters
      Initialised by: thermo = langevin(temp, dt, tau, econs)
      temp = temperature, default = 1.0
      dt = time step, default = 1.0
      tau = thermostat relaxation time, default = 1.0
      econs = conserved energy quantity, default = 0.0"""

   def get_et(self):
      return math.exp(-0.5*self.dt/self.tau)

   def get_K(self):
      return Constants.kb*self.temp*0.5
      
   def __init__(self, temp = 1.0, dt = 1.0, tau = 1.0, ethermo=0.0):
      super(ThermoSVR,self).__init__(temp,dt,ethermo)
      
      dset(self,"tau",depend_value(value=tau,name='tau'))
      dset(self,"et",  depend_value(name="et",func=self.get_et, dependencies=[dget(self,"tau"), dget(self,"dt")]))
      dset(self,"K",  depend_value(name="K",func=self.get_K, dependencies=[dget(self,"temp")]))
      
   def step(self):
      """Updates the atom velocities with a langevin thermostat"""
      
      K=np.dot(depstrip(self.p),depstrip(self.p)/depstrip(self.m))*0.5
      
      r1=self.prng.g
      if (self.ndof-1)%2==0: rg=2.0*self.prng.gamma((self.ndof-1)/2)
      else: rg=2.0*self.prng.gamma((self.ndof-2)/2)+self.prng.g**2
            
      alpha2=self.et+self.K/K*(1-self.et)*(r1**2+rg)+2.0*r1*math.sqrt(self.K/K*self.et*(1-self.et))
      alpha=math.sqrt(alpha2)
      if (r1+math.sqrt(2*K/self.K*self.et/(1-self.et))) <0: alpha*=-1

      self.ethermo+=K*(1-alpha2)
      self.p*=alpha

class ThermoPILE_G(ThermoPILE_L):    
   def __init__(self, temp = 1.0, dt = 1.0, tau = 1.0, ethermo=0.0):
      super(ThermoPILE_G,self).__init__(temp,dt,tau,ethermo)
      
   def bind(self, beads=None, prng=None, ndof=None):   
      # first binds as a local PILE, then substitutes the thermostat on the centroid
      super(ThermoPILE_G,self).bind(beads=beads,prng=prng,bindcentroid=False, ndof=ndof) 
            
      self._thermos[0]=ThermoSVR(temp=1, dt=1, tau=1) 
      
      t=self._thermos[0]      
      t.bind(pm=(beads.pnm[0,:],beads.m3[0,:]),prng=self.prng, ndof=ndof)
      deppipe(self,"temp", t, "temp")
      deppipe(self,"dt", t, "dt")
      deppipe(self,"tau", t, "tau")
      dget(self,"ethermo").add_dependency(dget(t,"ethermo"))
  
class ThermoGLE(Thermostat):     
   """Represent a GLE thermostat.
      Contains: dt = time step, econs = change in the kinetic energy due to the thermostat,
      tau = thermostat relaxation time, sm = sqrt(mass), A = friction matrix, C = static covariance,
      temp = temperature,  (T,S) = GLE propagator matrices
      Initialised by: thermo = ThermoGLE(temp, dt, A, C, econs)
      temp = temperature, default = 1.0
      A = friction, default 1x1 matrix value 1
      C = static covariance, default identity (sizeof(A)) x temp
      dt = time step, default = 1.0
      tau = thermostat relaxation time, default = 1.0
      econs = conserved energy quantity, default = 0.0"""

   def get_T(self):
      """Calculates T in p(0) = T*p(dt) + S*random.gauss()"""
      return math.exp(-0.5*self.dt/self.tau)
      
   def get_S(self):      
      """Calculates S in p(0) = T*p(dt) + S*random.gauss()"""
      return math.sqrt(Constants.kb*self.temp*(1-self.T**2))   
  
   def __init__(self, temp = 1.0, dt = 1.0, A = None, C = None, ethermo=0.0):
      super(ThermoLangevin,self).__init__(temp,dt,ethermo)
      
      dset(self,"tau",depend_value(value=tau,name='tau'))
      dset(self,"T",  depend_value(name="T",func=self.get_T, dependencies=[dget(self,"tau"), dget(self,"dt")]))
      dset(self,"S",  depend_value(name="S",func=self.get_S, dependencies=[dget(self,"temp"), dget(self,"T")]))
      
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
      elif type(thermo) is ThermoSVR: 
         self.kind.store("svr")
         self.tau.store(thermo.tau)         
      elif type(thermo) is ThermoPILE_L: 
         self.kind.store("pile_l")
         self.tau.store(thermo.tau)
      elif type(thermo) is ThermoPILE_G: 
         self.kind.store("pile_g")
         self.tau.store(thermo.tau)                  
      else: self.kind.store("unknown")      
      self.ethermo.store(thermo.ethermo)
      
   def fetch(self):
      if self.kind.fetch() == "langevin" : thermo=ThermoLangevin(tau=self.tau.fetch())
      elif self.kind.fetch() == "svr" : thermo=ThermoSVR(tau=self.tau.fetch())
      elif self.kind.fetch() == "pile_l" : thermo=ThermoPILE_L(tau=self.tau.fetch())
      elif self.kind.fetch() == "pile_g" : thermo=ThermoPILE_G(tau=self.tau.fetch())            
      else: thermo=Thermostat()
      thermo.ethermo=self.ethermo.fetch()
      return thermo
           
