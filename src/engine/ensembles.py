from utils.depend import *
from utils.restart import *
from utils import units
from thermostats import *
from barostats import *
import time

class RestartEnsemble(Restart):
   attribs={"type"  : (RestartValue, (str, "unknown")) }
   fields={"thermostat" : (RestartThermo, () ), "barostat" : (RestartBaro, () ), 
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
         ens=NVEEnsemble(dt=self.timestep.fetch(), temp=self.temperature.fetch(), nmstep=self.nm_propagator.fetch(), fixcom=self.fixcom.fetch())
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

   def __init__(self, dt, temp, nmstep=True):
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

      dset(self,"ntemp", depend_value(name='ntemp',func=self.get_ntemp,dependencies=[dget(self,"temp")]))
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
   def __init__(self, dt, temp, nmstep=True, fixcom=False):  # PIMD requires a temperature even for NVE (to define the spring constant)
      super(NVEEnsemble,self).__init__(dt=dt,temp=temp, nmstep=nmstep)
      self.fixcom=fixcom

   def rmcom(self):
      if (self.fixcom):
         pcom=np.zeros(3,float);
         
         p=depstrip(self.beads.p);
         na3=self.beads.natoms*3; nb=self.beads.nbeads;         
         m=depstrip(self.beads.m3)[:,0:na3:3];                  
         # computes center of mass
         for i in range(3): pcom[i]=p[:,i:na3:3].sum();

         if hasattr(self,"thermostat"):
            if hasattr(self.thermostat, "_thermos") : # in case of multiple thermostats, mangle with the centroid thermostat only
               self.thermostat._thermos[0].ethermo+=np.dot(pcom,pcom)/(2.0*self.beads[0].M*nb)
            else:
               self.thermostat.ethermo+=np.dot(pcom,pcom)/(2.0*self.beads[0].M*nb)

         # subtracts COM _velocity_
         pcom*=1.0/(nb*self.beads[0].M)
         print "COM VELOCITY: ", pcom
         for i in range(3): self.beads.p[:,i:na3:3]-=m*pcom[i]
         
   

   def pstep(self): 
      self.beads.p += self.forces.f * (self.dt*0.5)
   
   def qcstep(self):
      self.beads.qc += self.beads.pc/self.beads.atoms.m3*self.dt

   def qstep(self):
      """Velocity Verlet time step"""
      if self.beads.nbeads==1:
         pass
      else:
         # does normal modes propagation
         pq=np.zeros((2,self.beads.natoms*3),float)     #allocates a temporary to store pnm[i] and qnm[i]
         sm=depstrip(self.beads.sm3[0])                 #all the beads have the same masses
         for k in range(1,self.beads.nbeads):
            # works in mass-scaled coordinates so that all free-particle propagators are equal
            pq[0,:]=depstrip(self.beads.pnm[k])/sm
            pq[1,:]=depstrip(self.beads.qnm[k])*sm
            pq=np.dot(self.prop_pq[k],pq)
            self.beads.qnm[k]=pq[1,:]/sm         
            self.beads.pnm[k]=pq[0,:]*sm                    
      
   def step(self): 
      self.pstep()
      self.qcstep()
      self.qstep()
      self.pstep()
      self.rmcom()
      
class NVTEnsemble(NVEEnsemble):
   def __init__(self, dt, temp, nmstep=True, thermostat=Thermostat(), fixcom=False):
      super(NVTEnsemble,self).__init__(dt=dt,temp=temp, nmstep=nmstep, fixcom=fixcom)
      self.thermostat=thermostat
      
      self.dt=dt
      self.temp=temp
      self.ptime=0
      self.qtime=0
      self.ttime=0


   def bind(self,beads, cell, bforce,prng):
      super(NVTEnsemble,self).bind(beads, cell, bforce, prng)
      self.thermostat.bind(beads=self.beads,prng=prng,ndof=(3*(self.beads.natoms-1) if self.fixcom else None) )
      #merges the definitions of dt and temp
      deppipe(self,"ntemp", self.thermostat,"temp")
      deppipe(self,"dt", self.thermostat, "dt")

      dget(self,"econs").add_dependency(dget(self.thermostat, "ethermo"))
      
   def step(self): 
      self.ttime-=time.time()
      self.thermostat.step()
      self.rmcom()
      self.ttime+=time.time()
      self.ptime-=time.time()
      self.pstep()
      self.ptime+=time.time()
      
      self.qtime-=time.time()
      self.qcstep()
      self.qstep()
      self.qtime+=time.time()

      self.ptime-=time.time()
      self.pstep()
      self.ptime+=time.time()
      self.ttime-=time.time()
      self.thermostat.step()
      self.rmcom()      
      self.ttime+=time.time()


   def get_econs(self):
      """Calculates the conserved energy quantity for constant T ensembles"""
      return NVEEnsemble.get_econs(self)+self.thermostat.ethermo 

class NPTEnsemble(NVTEnsemble):
   def __init__(self, dt, temp, pext, nmstep = True, thermostat=Thermostat(), barostat=Barostat(), fixcom=False):
      super(NPTEnsemble,self).__init__(dt=dt, temp=temp, nmstep=nmstep, thermostat=thermostat, fixcom=fixcom)
      self.barostat=barostat
      
      #set the barostat pressure
      dset(self,"pext",depend_value(name="pext", value=pext) )
      deppipe(self, "pext", self.barostat, "pext")

   def bind(self, beads, cell, bforce, prng):
      super(NPTEnsemble,self).bind(beads, cell, bforce, prng)
      self.barostat.bind(beads, cell, bforce)

      # "binds" the dt, temp and p properties of the barostat to those of the ensemble. 
      # the thermostat has been already bound in NVTEnsemble init so we must make sure to copy all dependencies as well
      deppipe(self,"ntemp", self.barostat,"temp")
      deppipe(self,"ntemp", self.barostat.thermostat,"temp")
      deppipe(self, "dt", self.barostat, "dt")
            
      dget(self,"econs").add_dependency(dget(self.barostat.thermostat, "ethermo"))
      dget(self,"econs").add_dependency(dget(self.barostat, "pot"))
      dget(self,"econs").add_dependency(dget(self.thermostat, "temp"))
      dget(self,"econs").add_dependency(dget(self.cell, "kin"))
      dget(self,"econs").add_dependency(dget(self.cell, "V"))

   def get_econs(self):
      """Calculates the conserved energy quantity for constant T ensembles"""
      return NVTEnsemble.get_econs(self)+self.barostat.thermostat.ethermo+self.barostat.pot+ self.cell.kin - 2.0*Constants.kb*self.thermostat.temp*math.log(self.cell.V)
      
   def step(self):
      """Velocity Verlet time step"""
      self.thermostat.step()
      self.barostat.thermostat.step()
      self.rmcom()            
      self.barostat.pstep()
      self.barostat.qcstep()
      self.qstep()
      self.barostat.pstep()
      self.barostat.thermostat.step()
      self.thermostat.step()      
      self.rmcom()
                        
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
