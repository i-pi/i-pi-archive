import numpy as np
import math, random
from utils.depend import *
from utils.units  import *
from utils.io     import *
from atoms import *
from cell import *
#from forces import *

class Simulation(dobject):
   """Represents a simulation cell. Includes the cell parameters, 
      the atoms and the like."""   

   def __init__(self, atoms, cell, force, ensemble):
      self.atoms=atoms
      self.cell=cell
      self.force=force
      self.ensemble=ensemble
      self.bind()
      
   def bind(self):
      self.force.bind(self.atoms, self.cell)
      self.ensemble.bind(self.atoms, self.cell, self.force)
   
   
#      self.properties
      
#      
#      self.pcom=depend_array(name="pcom",value=np.zeros(3,float),
#         deps=depend_func(func=self.get_pcom,dependencies=[depget(self.atoms,"p"),depget(self.atoms,"m3")]) )
#          
#      
#      self.pot = depend_value(value=0.0,name="pot",deps=depend_func(func=(lambda : 0.0)))
#      self.vir = depend_array(value=numpy.zeros((3,3),float),name='vir', deps=depend_func(func=(lambda : 0.0), dependencies = [depget(self.cell,"V")]))

#      self.stress = depend(name='stress',deps=depend_func(func=self.get_stress,dependencies=[depget(self,"vir"), depget(self,"kstress")]))
#      self.press = depend(name='press',deps=depend_func(func=self.get_press,dependencies=[depget(self,stress)]))



#   def get_pcom(self):
#      p = numpy.zeros(3)      
#      for i in range(3):
#         p[i]=self.atoms.p[i:3*self.atoms.natoms:3].sum()
#      return p/self.atoms.m.sum()

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
