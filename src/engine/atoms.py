import numpy, math
from utils.depend import *

class Atom(object):
   """Represent an atom, with position, velocity, mass and related properties"""
    
#qpslice holds position and momentum data. 
#qslice holds a reference to the position data, pslice to the momentum data.
      
   def __init__(self, qpfslice, name="X", mass=1.0):
      # un-dependant properties
      self.q = depend(value=qpfslice[:,0],name='q')
      self.p = depend(value=qpfslice[:,1],name='p')
      self.f = depend(value=qpfslice[:,2],name='f')
      self.mass = depend(value=mass,name='mass')
      self.name = depend(value=name,name='name')
      
      # dependant properties (w/computation function)
      self.kin = depend(name='kin', func=self.get_kin, deplist=[self.p,self.mass])
      self.kstress = depend(name='kstress', func=self.get_kstress, deplist=[self.p,self.mass])      

   def __str__(self):
      return "    Name = %s\n    q = %s\n    p = %s\n    f = %s " % (self.name.get(), self.q.get(), self.p.get(), self.f.get())

   def get_kin(self):
      """Calculates the kinetic energy of the particle from the particle 
         momentum"""
      #print " [ upd. atom.kin ]", 
      ke = 0.0
      p=self.p.get()
      for i in range(3):
         ke += p[i]**2 
      return ke/(2.0*self.mass.get())

   def get_kstress(self):
      """Calculates the kinetic energy of the particle from the particle 
         momentum"""
      #print " [ upd. atom.kin ]", 
      p=self.p.get()
      return numpy.outer(p,p)/self.mass.get()
