import numpy as np
from utils.depend import *
from utils import units
from thermostats import *
from barostats import *
import time

#class EnsembleRestart(object):  
#   def __init__(self, ensemble=None): 
#      if not ensemble is
#   
#   def write(self, ensemble):
#      if   type(ensemble) is Ensemble:    self.kind="none"
#      elif type(ensemble) is NVEEnsemble: self.kind="nve"
#      else: raise TypeError("Unknown ensemble type")
#      
#      self.atoms=AtomsRestart(); 
#      self.atoms.write(simul.atoms)
#      
#      self.ensemble=EnsembleRestart();
#      self.ensemble.write(simul.ensemble)
#      
#   def read(self):
#      atoms=self.atoms.read()
#      sim=Simulation(atoms,cell=None,None,None)
#      return sim
#      

class Ensemble(dobject): 
   """General ensemble object, with no particle motion
      Contains: syst = System object, containing the atom and cell coordinates, 
      econs = conserved energy quantity.
      Initialised by: ens = ensemble(system)
      system is a System object, containing the atom and cell coordinates"""

   def __init__(self):
      dset(self,"econs",depend_value(name='econs',deps=depend_func(func=self.get_econs)) )
      self.timer=0.0
      
   def bind(self, atoms, cell, force):
      self.atoms=atoms
      self.cell=cell
      self.force=force
      depget(self,"econs").add_dependency(depget(self.atoms, "kin"))
      depget(self,"econs").add_dependency(depget(self.force, "pot"))      
      
   def step(self): 
      """Dummy routine which does nothing"""
      pass

   def get_econs(self):
      """Calculates the conserved energy quantity for constant E ensembles"""
      return self.atoms.kin+self.force.pot
      
      
class NVEEnsemble(Ensemble):
   def __init__(self, dt=1.0):
      super(NVEEnsemble,self).__init__()
      dset(self,"dt",depend_value(name='dt',value=dt))
   
   def step(self):
      """Velocity Verlet time step"""

      self.atoms.p += self.force.f * (self.dt*0.5)
      self.atoms.q += self.atoms.p/self.atoms.m3*self.dt
      self.atoms.p += self.force.f * (self.dt*0.5)


class NVTEnsemble(NVEEnsemble):
   def __init__(self, dt=None, temp=None, thermostat=Thermostat(), fixcom=False):
      super(NVTEnsemble,self).__init__(dt=dt)
      self.fixcom=fixcom
      self.thermostat=thermostat
      
      # binds options for dt and temperature of the thermostat to those in the ensemble
      # dt
      if not dt   is None: self.dt=dt      
      else: self.dt=2.0*self.thermostat.dt   
      dset(self.thermostat,"dt",   #this involves A LOT of piping
         depend_value(name="dt",deps=depend_func(func=self.get_halfdt,dependencies=[depget(self,"dt")],dependants=depget(self.thermostat,"dt")._dependants) ) )
      
      #temp
      dset(self, "temp", dget(self.thermostat,"temp"))
      if not temp is None: self.temp=temp
      
   def get_halfdt(self):  return self.dt*0.5
   
   def bind(self, atoms, cell, force):
      super(NVTEnsemble,self).bind(atoms,cell,force)
      self.thermostat.bind(atoms)

   def get_econs(self):
      """Calculates the conserved energy quantity for constant T ensembles"""
      return NVEEnsemble.get_econs(self)+self.thermostat.ethermo  
      
   def rm_com(self):
      if (self.fixcom):
         pcom=np.zeros(3,float);
         
         pcom[0]=self.atoms.px.sum();          pcom[1]=self.atoms.py.sum();          pcom[2]=self.atoms.pz.sum()
         self.thermostat.ethermo+=np.dot(pcom,pcom)/(2.0*self.atoms.M)
         
         pcom*=1.0/self.atoms.M
         self.atoms.px-=self.atoms.m*pcom[0];  self.atoms.py-=self.atoms.m*pcom[1];  self.atoms.pz-=self.atoms.m*pcom[2]; 
   
   def step(self):
      """Velocity Verlet time step"""

      self.timer-=time.clock()    
      self.thermostat.step()      
      self.timer+=time.clock()    
      self.rm_com()
      super(NVTEnsemble,self).step()
      self.timer-=time.clock()          
      self.thermostat.step()
      self.timer+=time.clock()          
      self.rm_com()
            
class NPTEnsemble(NVTEnsemble):
   def __init__(self, dt=None, temp=None, pext=None, thermostat=Thermostat(), barostat=Barostat(), fixcom=False ):
      super(NPTEnsemble,self).__init__(dt=dt, temp=temp, thermostat=thermostat, fixcom=fixcom)
      self.barostat=barostat
            
      # "binds" the dt, temp and p properties of the barostat to those of the ensemble. 
      # the thermostat has been already bound in NVTEnsemble init so we must make sure to copy all dependencies as well
      depcopy(self,"temp", self.barostat,"temp")
      depcopy(self,"temp", self.barostat.thermostat,"temp")
      depcopy(self, "dt", self.barostat, "dt")
            
      #dulicate pext & sext from the barostat 
      dset(self, "pext", dget(self.barostat,"pext"))
      if not pext is None: self.pext=pext
      
   def bind(self, atoms, cell, force):
      super(NPTEnsemble,self).bind(atoms, cell, force)
      self.barostat.bind(atoms, cell, force)

   def get_econs(self):
      """Calculates the conserved energy quantity for constant T ensembles"""
      return NVTEnsemble.get_econs(self)+self.barostat.thermostat.ethermo+self.barostat.pot+ self.cell.kin - 2.0*Constants.kb*self.thermostat.temp*math.log(self.cell.V)
      
   def step(self):
      """Velocity Verlet time step"""
      self.thermostat.step()
      self.rm_com()            
      self.barostat.step()
      self.thermostat.step()      
      self.rm_com()
                        
class NSTEnsemble(NVTEnsemble):
   def __init__(self, dt=None, temp=None, pext=None, sext=None, thermostat=Thermostat(), barostat=Barostat(), fixcom=False ):
      super(NSTEnsemble,self).__init__(dt=dt, temp=temp, thermostat=thermostat,fixcom=fixcom)
      self.barostat=barostat
            
      # "binds" the dt, temp and p properties of the barostat to those of the ensemble. 
      # the thermostat has been already bound in NVTEnsemble init so we must make sure to copy all dependencies as well
      depcopy(self,"temp", self.barostat,"temp")
      depcopy(self,"temp", self.barostat.thermostat,"temp")
      depcopy(self, "dt", self.barostat, "dt")
            
      #dulicate pext & sext from the barostat 
      dset(self, "pext", dget(self.barostat,"pext"))
      dset(self, "sext", dget(self.barostat,"sext"))      
      if not pext is None: self.pext=pext
      if not sext is None: self.sext=sext
      
   def bind(self, atoms, cell, force):
      super(NSTEnsemble,self).bind(atoms, cell, force)
      self.barostat.bind(atoms, cell, force)

   def get_econs(self):
      """Calculates the conserved energy quantity for constant T ensembles"""
      xv=0.0 # extra term stemming from the Jacobian
      for i in range(3):  xv+=math.log(self.cell.h[i,i])*(3-i)
      return NVTEnsemble.get_econs(self)++self.barostat.thermostat.ethermo+self.barostat.pot+self.cell.kin-2.0*Constants.kb*self.thermostat.temp*xv
      
   def step(self):
      """Velocity Verlet time step"""
      self.thermostat.step()         
      self.timer-=time.clock()       
      self.rm_com()            
      self.timer+=time.clock()        
      self.barostat.step()
      self.thermostat.step()      
      self.rm_com()      
