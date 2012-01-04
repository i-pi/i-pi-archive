import numpy as np
import math, random
from utils.depend import *
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
      
      self.property_dict["conserved"] = dget(self.ensemble,"econs")
      
      dset(self, "kin", depend_value(name="kin", func=self.get_kin, dependencies=[dget(self.beads,"kin"),dget(self.cell,"kin")]))
      self.property_dict["kinetic"] = dget(self,"kin")
      
      self.property_dict["potential"] = dget(self.forces,"pot")
     
      self.property_dict["V"] = dget(self.cell,"V")
      dset(self, "cell_params", depend_value(name="cell_params", func=self.get_cell_params, dependencies=[dget(self.cell, "h")]))
      self.property_dict["cell_parameters"] = dget(self,"cell_params")

      dset(self, "press", depend_value(name="press", func=self.get_press, dependencies=[dget(self.beads, "p"), dget(self.beads, "m3"), dget(self.forces, "vir"), dget(self.cell, "V")]))
      self.property_dict["pressure"] = dget(self,"press")

      
   def get_kin(self):          return self.beads.kin + self.cell.kin
   def get_time(self):         return self.simul.step * self.ensemble.dt
   def __getitem__(self,key):  return self.property_dict[key].get()
      
   def get_press(self):
      p = depstrip(self.atoms.p)
      m = depstrip(self.atoms.m3)
      vir = depstrip(self.forces.vir)
      return (np.dot(p,p/m) + np.trace(vir))/(3.0*self.cell.V)
   
   def get_cell_params(self):
      a, b, c, alpha, beta, gamma = h2abc(self.cell.h)
      return [a, b, c, alpha, beta, gamma]

