import numpy as np
import math, random
from utils.depend import *
from utils.units import Constants
from utils.mathtools import h2abc
from atoms import *
from cell import *
from ensembles import *
from forces import *

_DEFAULT_FINDIFF=1e-5
_DEFAULT_FDERROR=1e-9
_DEFAULT_MINFID=1e-12
class Properties(dobject):
   """A proxy to compute and output properties of the system, which may require putting together bits and pieces
      from different elements of the simulation."""   

   def __init__(self):
      self.property_dict = {}
      self.fd_delta=-_DEFAULT_FINDIFF
      
   def bind(self, simul):
      self.ensemble = simul.ensemble
      self.beads = simul.beads
      self.cell = simul.cell
      self.forces = simul.forces
      self.simul=simul      

      dset(self, "time", depend_value(name="time",  func=self.get_time, dependencies=[dget(self.simul,"step"), dget(self.ensemble,"dt")]))
      self.property_dict["time"] = dget(self,"time")
      
      dset(self, "econs", depend_value(name="econs", func=self.get_econs, dependencies=[dget(self.ensemble,"econs")]))
      self.property_dict["conserved"] = dget(self,"econs")
      
      dset(self, "kin", depend_value(name="kin", func=self.get_kin, dependencies=[dget(self.beads,"kin"),dget(self.cell,"kin")]))
      self.property_dict["kinetic_md"] = dget(self,"kin")
      
      dset(self, "pot", depend_value(name="pot", func=self.get_pot, dependencies=[dget(self.forces,"pot")]))      
      self.property_dict["potential"] = dget(self,"pot")

      dset(self, "temp", depend_value(name="temp", func=self.get_temp, dependencies=[dget(self.beads,"kin")]))      
      self.property_dict["temperature"] = dget(self,"temp")     
     
      self.property_dict["V"] = dget(self.cell,"V")
      dset(self, "cell_params", depend_value(name="cell_params", func=self.get_cell_params, dependencies=[dget(self.cell, "h")]))
      self.property_dict["cell_parameters"] = dget(self,"cell_params")

      dset(self, "stress", depend_value(name="stress", func=self.get_stress, dependencies=[dget(self.beads, "kstress"), dget(self.forces, "vir"), dget(self.cell, "V")]))
      self.property_dict["stress_md.xx"] = depend_value(name="scl_xx", dependencies=[dget(self, "stress")], func=(lambda : self.stress[0,0]) ) 
      
      dset(self, "press", depend_value(name="press", func=self.get_press, dependencies=[dget(self,"stress")]))
      self.property_dict["pressure_md"] = dget(self,"press")

      dset(self, "kin_cv", depend_value(name="kin_cv", func=self.get_kincv, dependencies=[dget(self.beads,"q"),dget(self.forces,"f"),dget(self.ensemble,"temp")]))
      self.property_dict["kinetic_cv"] = dget(self,"kin_cv")

      dset(self, "kstress_cv", depend_value(name="kstress_cv", func=self.get_kstresscv, dependencies=[dget(self.beads,"q"),dget(self.forces,"f"),dget(self.ensemble,"temp")]))
      dset(self, "stress_cv", depend_value(name="stress_cv", func=self.get_stresscv, dependencies=[dget(self,"kstress_cv"),dget(self.forces,"vir"), dget(self.cell, "V")]))
      self.property_dict["stress_cv.xx"] = depend_value(name="scv_xx", dependencies=[dget(self, "stress_cv")], func=(lambda : self.stress_cv[0,0]) ) 
      dset(self, "press_cv", depend_value(name="press_cv", func=self.get_presscv, dependencies=[dget(self,"stress_cv")]))
      self.property_dict["pressure_cv"] = dget(self,"press_cv")
      
      # creates dummy beads and forces to compute displaced and scaled estimators
      self.dbeads=simul.beads.copy()
      self.dforces=ForceBeads()
      self.dforces.bind(self.dbeads, self.simul.cell,  self.simul._forcemodel)     
      
      dset(self, "kin_yama", depend_value(name="kin_yama", func=self.get_kinyama, dependencies=[dget(self.beads,"q"),dget(self.ensemble,"temp")]))
      self.property_dict["kinetic_yamamoto"] = dget(self,"kin_yama")
      
      
   def get_kin(self):          return self.beads.kin/self.beads.nbeads
   def get_time(self):         return (1+self.simul.step)*self.ensemble.dt
   def __getitem__(self,key):  return self.property_dict[key].get()

   def get_pot(self):          return self.forces.pot/self.beads.nbeads
   def get_temp(self):         return self.beads.kin/(0.5*Constants.kb*(3*self.beads.natoms*self.beads.nbeads-(3 if self.ensemble.fixcom else 0))*self.beads.nbeads)
   def get_econs(self):        return self.ensemble.econs/self.beads.nbeads

   def get_stress(self):  return (self.forces.vir+self.beads.kstress)/self.cell.V                  
   def get_press(self):   return np.trace(self.stress)/3.0
   def get_stresscv(self):  return (self.forces.vir+self.kstress_cv)/self.cell.V                  
   def get_presscv(self):   return np.trace(self.stress_cv)/3.0
   
   def get_kincv(self):        
      kcv=0.0
      for b in range(self.beads.nbeads):
         kcv+=np.dot(depstrip(self.beads.q[b])-depstrip(self.beads.qc),depstrip(self.forces.f[b]))
      kcv*=-0.5/self.beads.nbeads
      kcv+=0.5*Constants.kb*self.ensemble.temp*(3*self.beads.natoms) 
      return kcv

   def get_kstresscv(self):        
      kst=np.zeros((3,3),float)
      q=depstrip(self.beads.q); qc=depstrip(self.beads.qc); na3=3*self.beads.natoms;
      for b in range(self.beads.nbeads):
         for i in range(3):
            for j in range(i,3):
               kst[i,j]+=np.dot(q[b,i:na3:3]-qc[i:na3:3],depstrip(self.forces.f[b])[j:na3:3])

      kst*=-1/self.beads.nbeads
      for i in range(3): kst[i,i]+=Constants.kb*self.ensemble.temp*(3*self.beads.natoms) 
      return kst

   def get_kinyama(self):              
      
      dbeta=abs(self.fd_delta)
      
      v0=self.pot
      while True: 
         splus=math.sqrt(1.0+dbeta); sminus=math.sqrt(1.0-dbeta)
         
         for b in range(self.beads.nbeads):
            self.dbeads[b].q=self.beads.centroid.q*(1.0-splus)+splus*self.beads[b].q      
         vplus=self.dforces.pot/self.beads.nbeads
         
         for b in range(self.beads.nbeads):
            self.dbeads[b].q=self.beads.centroid.q*(1.0-sminus)+sminus*self.beads[b].q      
         vminus=self.dforces.pot/self.beads.nbeads

         kyama=((1.0+dbeta)*vplus-(1.0-dbeta)*vminus)/(2*dbeta)-v0
         kyama+=0.5*Constants.kb*self.ensemble.temp*(3*self.beads.natoms) 
         if (self.fd_delta<0 and abs((vplus+vminus)/(v0*2)-1.0)>_DEFAULT_FDERROR and dbeta> _DEFAULT_MINFID):
            dbeta*=0.5; print "Reducing displacement in Yamamoto kinetic estimator"; continue
         else: break
         
      return kyama
         
   def get_cell_params(self):
      a, b, c, alpha, beta, gamma = h2abc(self.cell.h)
      return [a, b, c, alpha, beta, gamma]

