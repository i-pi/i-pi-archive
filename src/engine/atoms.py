import numpy, math
from utils.depend import *

class Atom(object):
   """Represent an atom, with position, velocity, mass and related properties.
      Contains: q = positions, p = momenta, f = forces, mass = mass,
      name = atom label, kin = kinetic energy, kstress = kinetic 
      contribution to the stress tensor.
      Initialised by: atom = Atom(qpfslice, name, mass)
      qpfslice[:,0] = positions
      qpfslice[:,1] = momenta
      qpfslice[:,2] = forces
      name is a small string giving the atom label, default = "X"
      mass is the mass, default = 1.0"""
      
   def __init__(self, qpfslice, name="X", mass=1.0):
      # un-dependent properties
      self.q = depend(value=qpfslice[:,0],name='q')
      self.p = depend(value=qpfslice[:,1],name='p')
      self.f = depend(value=qpfslice[:,2],name='f')
      self.mass = depend(value=mass,name='mass')
      self.name = depend(value=name,name='name')
      
      # dependent properties (w/computation function)
      self.kin = depend(name='kin', func=self.get_kin, deplist=[self.p,self.mass])
      self.kstress = depend(name='kstress', func=self.get_kstress, deplist=[self.p,self.mass])      

   def __str__(self):
      return "    Name = %s\n    q = %s\n    p = %s\n    f = %s " % (self.name.get(), self.q.get(), self.p.get(), self.f.get())

   def get_kin(self):
      """Calculates the kinetic energy of the particle from the particle 
         momentum"""

      p=self.p.get()
      return numpy.dot(p,p)/(2.0*self.mass.get())

   def get_kstress(self):
      """Calculates the contribution of the atom to the kinetic stress tensor"""

      p=self.p.get()
      kstress = numpy.zeros((3,3),float)
      for i in range(3):
         for j in range(i,3):
            kstress[i,j] = p[i]*p[j]
      return kstress/self.mass.get()

