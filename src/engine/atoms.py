import numpy
import math 

class Atom:
   """Represent an atom, with position, velocity, mass and related properties"""
    
#qpslice holds position and momentum data. 
#qslice holds a reference to the position data, pslice to the momentum data.

   def __init__(self, qpslice, name="X", mass=1.0):
      self.q = qpslice[:,0]
      self.p = qpslice[:,1]
      self.mass = mass
      self.name = name

   def __str__(self):
      return "    Name = %s, q = %s, p = %s " % (self.name, self.q, self.p)

   def pot(self):
      pot = 0.0
      return pot

   def kinetic(self):
      ke = 0.0
      for i in range(3):
         ke += self.p[i]**2 / (2.0*self.mass)
      return ke
   
