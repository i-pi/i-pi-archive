import numpy as np
from utils.depend import *
from utils.restart import *
from utils import units
from thermostats import *
from barostats import *
import time

class RestartEnsemble(Restart):
   attribs={"type"  : (RestartValue, (str, "unknown")) }
   fields={"thermostat" : (RestartThermo, () ), "barostat" : (RestartThermo, () ), 
           "timestep": (RestartValue, (float,"1.0")) ,
           "temperature" : (RestartValue, (float, 1.0)), "pressure" : (RestartValue, (float,"1.0")) ,
           "stress" : (RestartArray, (float, np.identity(3))), 
           "nm_propagator": (RestartValue, (bool, True)), "fixcom": (RestartValue, (bool, False)) }
   
   def store(self, ens):
      if type(ens) is NVEEnsemble:    
         self.type.store("nve"); tens=0
      elif type(ens) is NVTEnsemble:  
         self.type.store("nvt"); tens=1
      elif type(ens) is NPTEnsemble:  
         self.type.store("npt"); tens=2
      elif type(ens) is NSTEnsemble:  
         self.type.store("nst"); tens=3
      
      self.timestep.store(ens.dt)
      self.temperature.store(ens.temp)
      self.nm_propagator.store(ens.nmstep)
      
      if tens>0: 
         self.thermostat.store(ens.thermostat)
         self.fixcom.store(ens.fixcom)
      if tens>1:
         self.barostat.store(ens.barostat)
      if tens == 2:
         self.pressure.store(ens.pext)
      if tens == 3:
         self.stress.store(ens.pext)

   def fetch(self):
      if   self.type.fetch().upper() == "NVE" :
         ens=NVEEnsemble(dt=self.timestep.fetch(), temp=self.temperature.fetch(), nmstep=self.nm_propagator.fetch())
      elif self.type.fetch().upper() == "NVT" : 
         ens=NVTEnsemble(dt=self.timestep.fetch(), temp=self.temperature.fetch(), thermostat=self.thermostat.fetch(),  nmstep=self.nm_propagator.fetch(),
                        fixcom=self.fixcom.fetch() )
      elif self.type.fetch().upper() == "NPT" : 
         ens=NPTEnsemble(dt=self.timestep.fetch(), temp=self.temperature.fetch(), thermostat=self.thermostat.fetch(),   nmstep=self.nm_propagator.fetch(),
                        fixcom=self.fixcom.fetch(), pext=self.pressure.fetch(), barostat=self.barostat.fetch() )
      elif self.type.fetch().upper() == "NST" : 
         ens=NSTEnsemble(dt=self.timestep.fetch(), temp=self.temperature.fetch(), thermostat=self.thermostat.fetch(),  nmstep=self.nm_propagator.fetch(),
                        fixcom=self.fixcom.fetch(), sext=self.stress.fetch(), barostat=self.barostat.fetch() )
      return ens
      
class Ensemble(dobject): 
   """General ensemble object, with no particle motion
      Contains: syst = System object, containing the atom and cell coordinates, 
      econs = conserved energy quantity.
      Initialised by: ens = ensemble(system)
      system is a System object, containing the atom and cell coordinates"""

   def __init__(self, dt=1.0, temp=1.0, nmstep=True ):
      dset(self,"econs",depend_value(name='econs',func=self.get_econs) )
      dset(self, "temp", depend_value(name='temp',value=temp))       
      dset(self, "dt",   depend_value(name='dt',value=dt))
      self.nmstep=nmstep      
      self.timer=0.0
      
   def bind(self, beads, cell, bforce, prng):
      self.beads=beads
      self.cell=cell
      self.forces=bforce
      self.prng=prng

      dset(self, "ntemp", depend_value(name='ntemp',func=self.get_ntemp))           
      dget(self,"econs").add_dependency(dget(self.beads, "kin"))
      dget(self,"econs").add_dependency(dget(self.forces, "pot"))
      dget(self,"econs").add_dependency(dget(self.beads, "vpath"))
      
      #also computes frequencies and matrices needed for the propagator
      dset(self,"omegan",depend_value(name='omegan',func=self.get_omegan, dependencies=[dget(self,"ntemp")]) )
      dset(self,"omegan2",depend_value(name='omegan2',func=self.get_omegan2, dependencies=[dget(self,"omegan")]) )
      dset(self,"omegak",depend_array(name='omegak',value=np.zeros(self.beads.nbeads,float),func=self.get_omegak, dependencies=[dget(self,"omegan")]) )
      dset(self,"prop_pq",depend_array(name='prop_pq',value=np.zeros((self.beads.nbeads,2,2),float),func=self.get_prop_pq, 
                                      dependencies=[dget(self,"omegak"), dget(self,"dt")]) )
   def get_ntemp(self):   return self.temp*self.beads.nbeads
   def get_omegan(self):  return self.ntemp*units.Constants.kb/units.Constants.hbar
   def get_omegan2(self): return self.omegan**2
   def get_omegak(self):  return 2*self.omegan*np.array([math.sin(k*math.pi/self.beads.nbeads) for k in range(self.beads.nbeads) ])
   def get_prop_pq(self): 
      pqk=np.zeros((self.beads.nbeads,2,2), float)
      pqk[0]=np.array([[1,0],[self.dt,1]])
      for b in range(1,self.beads.nbeads):
         dtomegak=self.omegak[b]*self.dt
         c=math.cos(dtomegak); s=math.sin(dtomegak)
         pqk[b,0,0]=c;         pqk[b,1,1]=c
         pqk[b,0,1]=-s*self.omegak[b];         pqk[b,1,0]=s/self.omegak[b];
      return pqk

         
   def pstep(self): 
      """Dummy routine which does nothing"""
      pass

   def qstep(self): 
      """Dummy routine which does nothing"""
      pass

   def step(self): 
      """Dummy routine which does nothing"""
      pass

   def get_econs(self):
      """Calculates the conserved energy quantity for constant E ensembles"""
      return self.beads.kin+self.beads.vpath*self.omegan2+self.forces.pot
      
      
class NVEEnsemble(Ensemble):
   def __init__(self, dt=1.0, temp=1.0, nmstep=True, fixcom=False):  # PIMD requires a temperature even for NVE (to define the spring constant)
      super(NVEEnsemble,self).__init__(dt=dt,temp=temp, nmstep=nmstep)
      self.fixcom=fixcom

   def pstep(self): 
      self.beads.p += self.forces.f * (self.dt*0.5)
      pass
   
   def qstep(self):
      """Velocity Verlet time step"""
      if self.beads.nbeads==1:
         # just classical propagator
         self.beads.q += self.beads.p/self.beads.m3*self.dt
      elif not self.nmstep:
         # does cartesian propagation
         self.beads.p += self.beads.fpath* (self.omegan2*self.dt*0.5)
         self.beads.q += self.beads.p/self.beads.m3*self.dt
         self.beads.p += self.beads.fpath* (self.omegan2*self.dt*0.5)
      else:
         # does normal modes propagation
         pq=np.zeros((2,self.beads.natoms*3),float)     #allocates a temporary to store pnm[i] and qnm[i]
         sm=depstrip(self.beads.sm3[0])                 #all the beads have the same masses
         for k in range(self.beads.nbeads):
            # works in mass-scaled coordinates so that all free-particle propagators are equal
            pq[0,:]=depstrip(self.beads.pnm[k])/sm
            pq[1,:]=depstrip(self.beads.qnm[k])*sm
            pq=np.dot(self.prop_pq[k],pq)
            self.beads.qnm[k]=pq[1,:]/sm         
            self.beads.pnm[k]=pq[0,:]*sm                    
      
   def step(self): 
      self.pstep()
      self.qstep()
      self.pstep()
      
class NVTEnsemble(NVEEnsemble):
   def __init__(self, dt=None, temp=None, nmstep=True, thermostat=Thermostat(), fixcom=False):
      super(NVTEnsemble,self).__init__(dt=dt,temp=temp, nmstep=nmstep)
      self.thermostat=thermostat
      
      self.temp=self.thermostat.temp;       
      self.dt=self.thermostat.dt;
      if not dt is None:   self.dt=dt
      if not temp is None: self.temp=temp
      self.ptime=0
      self.qtime=0
      self.ttime=0


   def bind(self,beads, cell, bforce,prng):
      super(NVTEnsemble,self).bind(beads, cell, bforce, prng)
      self.thermostat.bind(pm=(self.beads.p.flatten(),self.beads.m3.flatten()),prng=prng)
      #merges the definitions of dt and temp
      for d in dget(self.thermostat,"temp")._dependants:  dget(self, "ntemp").add_dependant(d)
      dset(self.thermostat, "temp", dget(self,"ntemp"))
      for d in dget(self.thermostat,"dt")._dependants:    dget(self, "dt").add_dependant(d)
      dset(self.thermostat, "dt", dget(self,"dt"))

      
   def step(self): 
      self.ttime-=time.time()
      self.thermostat.step()
      self.ttime+=time.time()
      self.ptime-=time.time()
      self.pstep()
      self.ptime+=time.time()
      
      self.qtime-=time.time()
      self.qstep()
      self.qtime+=time.time()

      self.ptime-=time.time()
      self.pstep()
      self.ptime+=time.time()
      self.ttime-=time.time()
      self.thermostat.step()
      self.ttime+=time.time()


   def get_econs(self):
      """Calculates the conserved energy quantity for constant T ensembles"""
      return NVEEnsemble.get_econs(self)+self.thermostat.ethermo 

#class NVTEnsemble(NVEEnsemble):
#   def __init__(self, dt=None, temp=None, thermostat=Thermostat(), fixcom=False):
#      super(NVTEnsemble,self).__init__(dt=dt)
#      self.fixcom=fixcom
#      self.thermostat=thermostat
#      
#      # binds options for dt and temperature of the thermostat to those in the ensemble
#      # dt
#      if not dt   is None: self.dt=dt      
#      else: self.dt=2.0*self.thermostat.dt   
#      dset(self.thermostat,"dt",   #this involves A LOT of piping
#         depend_value(name="dt",deps=depend_func(func=self.get_halfdt,dependencies=[depget(self,"dt")],dependants=depget(self.thermostat,"dt")._dependants) ) )
#      
#      #temp
#      dset(self, "temp", dget(self.thermostat,"temp"))
#      if not temp is None: self.temp=temp
#      
#   def get_halfdt(self):  return self.dt*0.5
#   
#   def bind(self, atoms, cell, force):
#      super(NVTEnsemble,self).bind(atoms,cell,force)
#      self.thermostat.bind(atoms)

#   def get_econs(self):
#      """Calculates the conserved energy quantity for constant T ensembles"""
#      return NVEEnsemble.get_econs(self)+self.thermostat.ethermo  
#      
#   def rm_com(self):
#      if (self.fixcom):
#         pcom=np.zeros(3,float);
#         
#         pcom[0]=self.atoms.px.sum();          pcom[1]=self.atoms.py.sum();          pcom[2]=self.atoms.pz.sum()
#         self.thermostat.ethermo+=np.dot(pcom,pcom)/(2.0*self.atoms.M)
#         
#         pcom*=1.0/self.atoms.M
#         self.atoms.px-=self.atoms.m*pcom[0];  self.atoms.py-=self.atoms.m*pcom[1];  self.atoms.pz-=self.atoms.m*pcom[2]; 
#   
#   def step(self):
#      """Velocity Verlet time step"""

#      self.timer-=time.clock()    
#      self.thermostat.step()      
#      self.timer+=time.clock()    
#      self.rm_com()
#      super(NVTEnsemble,self).step()
#      self.timer-=time.clock()          
#      self.thermostat.step()
#      self.timer+=time.clock()          
#      self.rm_com()
#            
#class NPTEnsemble(NVTEnsemble):
#   def __init__(self, dt=None, temp=None, pext=None, thermostat=Thermostat(), barostat=Barostat(), fixcom=False ):
#      super(NPTEnsemble,self).__init__(dt=dt, temp=temp, thermostat=thermostat, fixcom=fixcom)
#      self.barostat=barostat
#            
#      # "binds" the dt, temp and p properties of the barostat to those of the ensemble. 
#      # the thermostat has been already bound in NVTEnsemble init so we must make sure to copy all dependencies as well
#      depcopy(self,"temp", self.barostat,"temp")
#      depcopy(self,"temp", self.barostat.thermostat,"temp")
#      depcopy(self, "dt", self.barostat, "dt")
#            
#      #dulicate pext & sext from the barostat 
#      dset(self, "pext", dget(self.barostat,"pext"))
#      if not pext is None: self.pext=pext
#      
#   def bind(self, atoms, cell, force):
#      super(NPTEnsemble,self).bind(atoms, cell, force)
#      self.barostat.bind(atoms, cell, force)

#   def get_econs(self):
#      """Calculates the conserved energy quantity for constant T ensembles"""
#      return NVTEnsemble.get_econs(self)+self.barostat.thermostat.ethermo+self.barostat.pot+ self.cell.kin - 2.0*Constants.kb*self.thermostat.temp*math.log(self.cell.V)
#      
#   def step(self):
#      """Velocity Verlet time step"""
#      self.thermostat.step()
#      self.rm_com()            
#      self.barostat.step()
#      self.thermostat.step()      
#      self.rm_com()
#                        
#class NSTEnsemble(NVTEnsemble):
#   def __init__(self, dt=None, temp=None, pext=None, sext=None, thermostat=Thermostat(), barostat=Barostat(), fixcom=False ):
#      super(NSTEnsemble,self).__init__(dt=dt, temp=temp, thermostat=thermostat,fixcom=fixcom)
#      self.barostat=barostat
#            
#      # "binds" the dt, temp and p properties of the barostat to those of the ensemble. 
#      # the thermostat has been already bound in NVTEnsemble init so we must make sure to copy all dependencies as well
#      depcopy(self,"temp", self.barostat,"temp")
#      depcopy(self,"temp", self.barostat.thermostat,"temp")
#      depcopy(self, "dt", self.barostat, "dt")
#            
#      #dulicate pext & sext from the barostat 
#      dset(self, "pext", dget(self.barostat,"pext"))
#      dset(self, "sext", dget(self.barostat,"sext"))      
#      if not pext is None: self.pext=pext
#      if not sext is None: self.sext=sext
#      
#   def bind(self, atoms, cell, force):
#      super(NSTEnsemble,self).bind(atoms, cell, force)
#      self.barostat.bind(atoms, cell, force)

#   def get_econs(self):
#      """Calculates the conserved energy quantity for constant T ensembles"""
#      xv=0.0 # extra term stemming from the Jacobian
#      for i in range(3):  xv+=math.log(self.cell.h[i,i])*(3-i)
#      return NVTEnsemble.get_econs(self)++self.barostat.thermostat.ethermo+self.barostat.pot+self.cell.kin-2.0*Constants.kb*self.thermostat.temp*xv
#      
#   def step(self):
#      """Velocity Verlet time step"""
#      self.thermostat.step()         
#      self.timer-=time.clock()       
#      self.rm_com()            
#      self.timer+=time.clock()        
#      self.barostat.step()
#      self.thermostat.step()      
#      self.rm_com()      
