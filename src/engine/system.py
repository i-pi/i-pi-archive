import numpy as np
import math, random
from utils.depend import *
from utils.units  import *
from utils.io     import *
from atoms import *
from cell import *
#from forces import *

class System(dobject):
   """Represents a simulation cell. Includes the cell parameters, 
      the atoms and the like."""   

   def __init__(self, natoms):
      self.natoms=natoms
      self.q=depend_array(name="q",value=np.zeros(3*natoms, float))
      self.p=depend_array(name="p",value=np.zeros(3*natoms, float))
      self.f=depend_array(name="f",value=np.zeros(3*natoms, float))
      self.m3=depend_array(name="m3",value=np.zeros(3*natoms, float))      
      self.__dict__["m"]=self.m3[0:3*self.natoms:3]
      
      self.cell=Cell()
      self.atoms=[ Atom(self, i) for i in range(natoms) ]
      
      self.pcom=depend_array(name="pcom",value=np.zeros(3,float),
         deps=depend_func(func=self.get_pcom,dependencies=[depget(self,"p"),depget(self,"m3")]) )
      self.kin=depend_value(name="kin",deps=depend_func(func=self.get_kin,dependencies=[depget(self,"p"),depget(self,"m3")]) )         
      
      self.pot = depend_value(value=0.0,name="pot",deps=depend_func(func=(lambda : 0.0)))
      self.vir = depend_array(value=numpy.zeros((3,3),float),name='vir', deps=depend_func(func=(lambda : 0.0), dependencies = [depget(self.cell,"V")]))

#      self.stress = depend(name='stress',deps=depend_func(func=self.get_stress,dependencies=[depget(self,"vir"), depget(self,"kstress")]))
#      self.press = depend(name='press',deps=depend_func(func=self.get_press,dependencies=[depget(self,stress)]))

   def __setattr__(self, name, value):
      if (name == "m"):
         dget(self,"m3")[0:3*self.natoms:3]=value
         dget(self,"m3")[1:3*self.natoms:3]=value         
         dget(self,"m3")[2:3*self.natoms:3]=value         
         return value
      else:
         return super(System,self).__setattr__(name,value)

   def get_pcom(self):
      p = numpy.zeros(3)      
      for i in range(3):
         p[i]=self.p[i:3*self.natoms:3].sum()
      return p/self.m.sum()

   def get_kin(self):
      """Calculates the total kinetic energy of the system,
      by summing the atomic contributions"""
      return 0.5*np.dot(self.p,self.p/self.m3)

#   def get_tempk(self):  TODO put this somewhere more appropriate (COM removal, etc)
#      """Calculates an estimate of the temperature from the kinetic energy"""
#      return 2.0*self.kin/(*units.Constants.kb*len(self.atoms))

#   def get_stress(self):
#      """Calculates the internal stress tensor"""
#   
#      return self.kstress.get()+self.vir.get()
#   
#   def get_press(self):
#      """Calculates the internal pressure scalar"""

#      return numpy.trace(self.stress.get())/3.0
#      
#   def get_kstress(self):
#      """Calculates the kinetic stress tensor"""

#      ks=numpy.zeros((3,3),float)
#      for at in self.atoms:
#         ks += at.kstress.get()
#      ks/=self.cell.V.get()
#      return ks      
