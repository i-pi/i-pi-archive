import numpy as np
from utils.depend import *
from utils import units

class Atom(dobject):
   """Represent an atom, with position, velocity, mass and related properties."""
      
   def get_kin(self):
      return np.dot(self.p,self.p)/(2.0*self.m)
      
   def __init__(self, system, index, label="X"):
      self.name=label
      self.p=system.p[3*index:3*index+3]
      self.q=system.q[3*index:3*index+3]
      self.f=system.f[3*index:3*index+3]
      self.__dict__["m"]=system.m3[3*index:3*index+3]            
      self.kin=depend_value(name="kin", deps=depend_func(func=self.get_kin, dependants=[depget(self,"p"),depget(system,"m3")]) )

   def __getattribute__(self, name):
      if (name == "m"):
         return dget(self,"m")[0]
      else:
         return super(Atom,self).__getattribute__(name)

   def __setattr__(self, name, value):
      if (name == "m"):
         dget(self,"m")[0:3]=value
         return value
      else:
         return super(Atom,self).__setattr__(name,value)

