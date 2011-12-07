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

class Atoms(dobject):
   """Represents a simulation cell. Includes the cell parameters, 
      the atoms and the like."""   

   def __init__(self, natoms):
      self.natoms=natoms
      dset(self,"q",depend_array(name="q",value=np.zeros(3*natoms, float)) ) 
      dset(self,"p",depend_array(name="p",value=np.zeros(3*natoms, float)) )
      dset(self,"m3",depend_array(name="m3",value=np.zeros(3*natoms, float)) )
      dset(self,"m",self.m3[0:3*self.natoms:3])
      
      self._alist=[ Atom(self, i) for i in range(natoms) ]
      dset(self,"kin",depend_value(name="kin",deps=depend_func(func=self.get_kin,dependencies=[depget(self,"p"),depget(self,"m3")])) )
      
   def __getitem__(self,index):
      return self._alist[index]

   def __setitem__(self,index,value):
      self._alist[index].p=value.p
      self._alist[index].q=value.q 
      self._alist[index].m=value.m

   def __setattr__(self, name, value):
      if (name == "m"):
         dget(self,"m3")[0:3*self.natoms:3]=value
         dget(self,"m3")[1:3*self.natoms:3]=value         
         dget(self,"m3")[2:3*self.natoms:3]=value         
         return value
      else:
         return super(Atoms,self).__setattr__(name,value)
         
         
   def get_kin(self):
      """Calculates the total kinetic energy of the system,
      by summing the atomic contributions"""
      return 0.5*np.dot(self.p,self.p/self.m3)         
