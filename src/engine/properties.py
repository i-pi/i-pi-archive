import numpy as np
import math, random
from utils.depend import *
from utils.units import Constants
from utils.mathtools import h2abc
from atoms import *
from cell import *
from ensembles import *
from forces import *

class Properties(dobject):
   """A proxy to compute and output properties of the system, which may require putting together bits and pieces
      from different elements of the simulation."""   

   def __init__(self):
      self.property_dict = {}
      
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
      self.property_dict["kinetic"] = dget(self,"kin")
      
      dset(self, "pot", depend_value(name="pot", func=self.get_pot, dependencies=[dget(self.forces,"pot")]))      
      self.property_dict["potential"] = dget(self,"pot")

      dset(self, "temp", depend_value(name="temp", func=self.get_temp, dependencies=[dget(self.beads,"kin")]))      
      self.property_dict["temperature"] = dget(self,"temp")     
     
      self.property_dict["V"] = dget(self.cell,"V")
      dset(self, "cell_params", depend_value(name="cell_params", func=self.get_cell_params, dependencies=[dget(self.cell, "h")]))
      self.property_dict["cell_parameters"] = dget(self,"cell_params")

      dset(self, "press", depend_value(name="press", func=self.get_press, dependencies=[dget(self.beads, "p"), dget(self.beads, "m3"), dget(self.forces, "vir"), dget(self.cell, "V")]))
      self.property_dict["pressure"] = dget(self,"press")

      dset(self, "kin_cv", depend_value(name="kin_cv", func=self.get_kincv, dependencies=[dget(self.beads,"q"),dget(self.forces,"f"),dget(self.ensemble,"temp")]))
      self.property_dict["kinetic_cvirial"] = dget(self,"kin_cv")

      
   def get_kin(self):          return self.beads.kin/self.beads.nbeads
   def get_time(self):         return (1+self.simul.step)*self.ensemble.dt
   def __getitem__(self,key):  return self.property_dict[key].get()

   def get_pot(self):          return self.forces.pot/self.beads.nbeads
   def get_temp(self):         return self.beads.kin/(0.5*Constants.kb*(3*self.beads.natoms*self.beads.nbeads-(3 if self.ensemble.fixcom else 0))*self.beads.nbeads)
   def get_econs(self):        return self.ensemble.econs/self.beads.nbeads

   def get_kincv(self):        
      kcv=0.0
      for b in range(self.beads.nbeads):
         kcv+=np.dot(depstrip(self.beads.q[b])-depstrip(self.beads.qc),depstrip(self.forces.f[b]))
      kcv*=-0.5/self.beads.nbeads
      kcv+=0.5*Constants.kb*self.ensemble.temp*(3*self.beads.natoms) 
      return kcv

      
   def get_press(self):
      p = depstrip(self.atoms.p)
      m = depstrip(self.atoms.m3)
      vir = depstrip(self.forces.vir)
      return (np.dot(p,p/m) + np.trace(vir))/(3.0*self.cell.V)
   
   def get_cell_params(self):
      a, b, c, alpha, beta, gamma = h2abc(self.cell.h)
      return [a, b, c, alpha, beta, gamma]

