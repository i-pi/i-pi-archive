import numpy, math, random
import upper_T
from utils.depend import *
from utils import units

class Cell(object):
   """Represents the simulation cell in a periodic system
      Contains: h = lattice vector matrix, p = lattice momentum matrix (NST),
      pc = lattice scalar momentum (NPT), w = barostat mass, pext = external 
      pressure tensor, V = volume, ih = inverse lattice matrix, 
      (h0, ih0, V0) = as above, but for reference cell (pext = 0), 
      kin = kinetic energy, pot = strain potential, 
      strain = strain tensor, piext = external stress tensor
      Initialised by: cell = Cell(h, w, h0, pext)
      h is the lattice vector matrix
      w = barostat mass, default = 1.0
      h0 = reference cell, default = h
      pext = external pressure tensor, default = 0"""

   @classmethod 
   def fromSidesAngles(cls, mycell=[1,1,1,math.pi/2.0,math.pi/2.0,math.pi/2.0], w=1.0, h0=None, pext = numpy.zeros((3,3),float)):
      """A different initialiser, which uses lattice vector lengths and angles
         to calculate the lattice vector matrix.
         Initialised by: cell = Cell.fromSidesAngles(mycell, w, h0, pext)
         All variables as for the standard initialisation, except
         mycell = [a, b, c, alpha, beta, gamma], default = [1,1,1,90,90,90]"""

      a, b, c, alpha, beta, gamma = mycell
      return cls(h=abc2h(a,b,c,alpha,beta,gamma), w = w, h0 = h0)

   def __init__(self, h = numpy.identity(3, float), w = 1.0, h0=None, pext = numpy.zeros((3,3),float) ):      
      #un-dependent properties
      self.h = depend(value = h, name = 'h')
      self.p = depend(value = numpy.zeros((3,3), float), name = 'p')
      self.pc = depend(value = 0.0, name = 'pc')
      self.w = depend(value = w, name='w')
      
      if (h0 is None): 
         h0=numpy.copy(h)
      self.h0 = depend(name='h0', value = h0)
      
      self.pext = depend(name = 'pext', value = pext)

      #dependent properties
      self.V = depend(name='V', func=self.get_volume)
      self.ih = depend(name='ih', func=self.get_ih)
      self.h.add_dependant(self.V);  self.h.add_dependant(self.ih)

      self.V0 = depend(name='V0', func=self.get_vol0)
      self.ih0 = depend(name='ih', func=self.get_ih0)
      self.h0.add_dependant(self.V0); self.h0.add_dependant(self.ih0)

      self.kin = depend(name='kin', func=self.get_kin) 
      self.p.add_dependant(self.kin);  self.w.add_dependant(self.kin);  self.pc.add_dependant(self.kin);
      
      self.strain = depend(name = 'strain', func=self.get_strain)
      self.h.add_dependant(self.strain); self.ih0.add_dependant(self.strain);
      self.pot = depend(name = 'pot', func = self.get_pot)
      self.pext.add_dependant(self.pot); self.V0.add_dependant(self.pot); self.strain.add_dependant(self.pot);
      
      self.piext = depend(name='piext', func=self.get_piext, deplist=[self.pext, self.V, self.V0, self.ih0, self.h])
   
   def init_velocity(self, temp):
      """Initialises the cell velocity matrix according to the 
         Maxwell-Boltzmann distribution"""

      sigma = math.sqrt(self.w.get()*units.kb*temp)
      for i in range(3):
         for j in range(i,3):
            self.p.get()[i,j] = random.gauss(0.0, sigma)
      self.p.taint(taintme=False)

   def get_volume(self):
      """Calculates the volume of the unit cell, assuming an upper-triangular
         lattice vector matrix"""

      h = self.h.get()
      return ut_det(h)
        
   def get_vol0(self):
      """Calculates the volume of the reference cell, assuming an 
         upper-triangular lattice vector matrix"""

      h0 = self.h0.get()
      return ut_det(h0)      

   def get_ih(self):
      """Inverts a 3*3 (upper-triangular) cell matrix"""

      h = self.h.get()      
      return ut_inverse(h)

   def get_ih0(self):
      """Inverts a 3*3 (upper-triangular) reference cell matrix"""

      h0 = self.h0.get()      
      return ut_inverse(h0)

   def get_kin(self):
      """Calculates the kinetic energy of the cell from the cell parameters"""

      p=self.p.get()
      ke = 0.0
      for i in range(3):
         for j in range(i,3):
            ke += p[i, j]**2
      ke += self.pc.get()**2            
            
      ke /= 2*self.w.get()
      return ke

   def get_strain(self):
      """Computes the strain tensor from the unit cell and reference cell"""

      root = numpy.dot(self.h.get(), self.ih0.get())
      eps = numpy.dot(numpy.transpose(root), root) - numpy.identity(3, float)
      eps *= 0.5
      return eps

   def get_pot(self):
      """Calculates the elastic strain energy of the cell"""

      return self.V0.get()*numpy.trace(numpy.dot(self.pext.get(), self.strain.get()))

   def get_piext(self):
      """Calculates the external stress tensor"""
      root = numpy.dot(self.h.get(), self.ih0.get())
      pi = numpy.dot(root, self.pext.get())
      pi = numpy.dot(pi, numpy.transpose(root))
      pi *= self.V0.get()/self.V.get()
      return pi

   def __str__(self):
      h=self.h.get(); p = self.p.get()
      return "    h1 = %s\n    h2 = %s\n    h3 = %s\n\n    p1 = %s\n    p2 = %s\n    p3 = %s\n\n    w = %s, volume = %s\n" % (
         h[:,0], h[:,1], h[:,2], p[:,0], p[:,1], p[:,2], self.w.get(), self.V.get() )
      
   def apply_pbc(self, atom):
      """Uses the minimum image convention to return a particle to the
         unit cell"""

      s=numpy.dot(self.ih,atom.q.get())
      for i in range(3):
         s[i] = s[i] - round(s[i])
      new_pos = numpy.dot(self.h,s)
      return new_pos

   def minimum_distance(self, atom1, atom2):
      """Takes two atoms and tries to find the smallest vector between two 
         images. This is only rigorously accurate in the case of a cubic cell,
         but gives the correct results as long as the cut-off radius is defined
         as smaller than the smallest width between parallel faces."""

      s = numpy.dot(self.ih.get(),atom1.q.get() - atom2.q.get())
      for i in range(3):
         s[i] -= round(s[i])
      return numpy.dot(self.h.get(), s)

      
def h2abc(h):
   """Returns a description of the cell in terms of the length of the 
      lattice vectors and the angles between them."""
   
   a=h[0,0]
   b=math.sqrt(h[0,1]**2+h[1,1]**2)
   c=math.sqrt(h[0,2]**2+h[1,2]**2+h[2,2]**2)
   gamma=math.acos(h[0,1]/b)
   beta=math.acos(h[0,2]/c)
   alpha = math.acos(numpy.dot(h[:,1], h[:,2])/(b*c))

   return a, b, c, alpha, beta, gamma

def abc2h(a, b, c, alpha, beta, gamma):
   """Returns a cell matrix given a description in terms of the vector lengths
      and the angles in between"""

   h = numpy.zeros((3,3) ,float)
   h[0,0] = a
   h[0,1] = b * math.cos(gamma)
   h[0,2] = c * math.cos(beta)
   h[1,1] = b * math.sin(gamma)
   h[1,2] = (b*c*math.cos(alpha) - h[0,1]*h[0,2]) / h[1,1]
   h[2,2] = math.sqrt(c**2 - h[0,2]**2 - h[1,2]**2)
   return h
   
   
def ut_inverse(h):
   """Inverts a 3*3 upper-triangular matrix"""

   ih = numpy.zeros((3,3), float)
   for i in range(3):
      ih[i,i] = 1.0/h[i,i]
   ih[0,1] = -ih[0,0]*h[0,1]*ih[1,1]
   ih[1,2] = -ih[1,1]*h[1,2]*ih[2,2]
   ih[0,2] = -ih[1,2]*h[0,1]*ih[0,0]-ih[0,0]*h[0,2]*ih[2,2]
   return ih
   
def ut_det(h):
   """Finds the determinant of a 3*3 upper-triangular matrix"""

   return h[0,0]*h[1,1]*h[2,2]
