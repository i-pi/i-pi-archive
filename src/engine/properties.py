import numpy as np
import math, random
from utils.depend import *
from utils.mathtools import h2abc
from atoms import *
from cell import *
from ensembles import *
from forces import *

class Properties(dobject):
   """Represents a simulation cell. Includes the cell parameters, 
      the atoms and the like."""   

   def __init__(self):
      self.property_dict = {}
      
   def bind(self, ensemble, atoms, cell, force):
      self.ensemble = ensemble
      self.atoms = atoms
      self.cell = cell
      self.force = force
      self.property_dict["econs"] = dget(self.ensemble,"econs")
      dset(self, "kin", depend_value(name="kin", deps=depend_func(func=self.get_kin, dependencies=[depget(self.atoms,"kin"),depget(self.cell,"kin")])))
      self.property_dict["kin"] = dget(self,"kin")
      self.property_dict["V"] = dget(self.cell,"V")
      self.property_dict["pot"] = dget(self.force,"pot")
      dset(self, "cell_params", depend_value(name="cell_params", deps=depend_func(func=self.get_cell_params, dependencies=[depget(self.cell, "h")])))
      self.property_dict["cell_parameters"] = dget(self,"cell_params")
      dset(self, "press", depend_value(name="press", deps=depend_func(func=self.get_press, dependencies=[depget(self.atoms, "p"), depget(self.atoms, "m3"), depget(self.force, "vir"), depget(self.cell, "V")])))
      self.property_dict["pressure"] = dget(self,"press")

   def get_kin(self):
      return self.atoms.kin + self.cell.kin
   
   def get_press(self):
      p = depstrip(self.atoms.p)
      m = depstrip(self.atoms.m3)
      vir = depstrip(self.force.vir)
      return (np.dot(p,p/m) + np.trace(vir))/(3.0*self.cell.V)
   
   def get_cell_params(self):
      a, b, c, alpha, beta, gamma = h2abc(self.cell.h)
      return [a, b, c, alpha, beta, gamma]

