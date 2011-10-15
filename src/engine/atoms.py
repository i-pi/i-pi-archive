import numpy, math
from utils.depend import *

class Atom(object):
   """Represent an atom, with position, velocity, mass and related properties"""
    
#qpslice holds position and momentum data. 
#qslice holds a reference to the position data, pslice to the momentum data.

   @depend()
   def q(self): pass
   @depend()
   def p(self): pass
   @depend()
   def f(self): pass
   @depend()
   def name(self): pass
   @depend()
   def mass(self): pass
      
   def __init__(self, qpfslice, name="X", mass=1.0):
      self.q = qpfslice[:,0]      
      self.p = qpfslice[:,1]
      self.f = qpfslice[:,2]
      self.mass = mass
      self.name = name

   def __str__(self):
      return "    Name = %s\n    q = %s\n    p = %s\n    f = %s " % (self.name, self.q, self.p, self.f)

   def pot(self):
      pot = 0.0
      return pot

   @depend([p,mass])
   def kinetic(self):
      """Calculates the kinetic energy of the particle from the particle 
         momentum"""
      ke = 0.0
      for i in range(3):
         ke += self.p[i]**2 / (2.0*self.mass)
      return ke
      
   @depend([p]) 
   def palias(self): 
      return self.p
   
