import numpy
import math

def compute_ih(h) :
   """Inverts a (upper-triangular) cell matrix"""
   ih = numpy.zeros((3,3), float)
   for i in range(3):
      ih[i,i] = 1.0/h[i,i]
   ih[0,1] = -ih[0,0]*h[0,1]/h[1,1]
   ih[1,2] = -ih[1,1]*h[1,2]/h[2,2]
   ih[0,2] = -ih[1,2]*h[0,1]*ih[1,1]-ih[0,0]*h[0,2]*ih[2,2]
   return ih

def compute_strain(h, ih_0):
   """Computes the strain tensor from the unit cell and reference cell"""
   root = numpy.dot(h, ih_0)
   eps = numpy.dot(numpy.transpose(root), root) - numpy.identity(3, float)
   eps /= 2
   return eps   

def volume(h):
   """Calculates the volume of the unit cell, assuming an upper-triangular
      unit vector matrix"""
   #return numpy.inner(self.h[0:3,0], numpy.cross(self.h[0:3,1], self.h[0:3,2])) # general matrix
   return h[0,0]*h[1,1]*h[2,2]   # upper-triangular matrix

class Cell(object):
   """Represents the simulation cell in a periodic system"""

#h is the lattice basis matrix, which will hold the basis vectors in column form
#h is going to be upper triangular, i.e. a_1=(a1x,0,0); a_2=(a2x,a2y,0); a_3=(a3x,a3y,a3z)
#p will hold the box momenta, in the same form as for h
#w is probably the barostat mass, or something similar. 
  
   @property
   def h(self): 
      return self.__h
   
   @h.setter
   def h(self, newh): 
      self.__h = newh
      self.__ih = numpy.zeros((3,3),float)
      self.__taint_ih = True
      self.__taint_eps = True

   @property
   def p(self): 
      return self.__p
   
   @p.setter
   def p(self, newp): 
      self.__p=newp

   @property
   def ih(self):       
      if (self.__taint_ih) :
         self.__ih=compute_ih(self.__h)
         self.__taint_ih = False
      return self.__ih

   @property
   def strain(self):
      if (self.__taint_eps):
         print "New eps formed"
         self.__eps = compute_strain(self.__h, self.__ih_0)
         self.__taint_eps = False
      return self.__eps
      
   
   def __init__(self, cell = [ 1, 1, 1, math.pi/2, math.pi/2, math.pi/2], P_ext = numpy.zeros(3, float) ):
      
      a, b, c, alpha, beta, gamma = cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]
      self.__h = abc2h(a, b, c, alpha, beta, gamma)
      self.__p = numpy.zeros((3,3) ,float)
      self.__taint_ih = True
      self.w = 1.0
      self.__P_ext = P_ext
      self.__h_0 = numpy.identity(3, float) #needs to be the unstrained cell
      self.__ih_0 = compute_ih(self.__h_0)
      self.__V_0 = volume(self.__h_0)

   def __str__(self):
      return "    h1 = %s, h2 = %s, h3 = %s \n    p1 = %s, p2 = %s, p3 = %s \n    w = %s, volume = %s" % (self.h[:,0], self.h[:,1], self.h[:,2], self.p[:,0], self.p[:,1], self.p[:,2], self.w, volume(self.__h))
      
   def pot(self):
      """Calculates the elastic strain energy of the cell"""
      pe = self.__V_0*numpy.trace(numpy.dot(self.__P_ext, self.strain))
      return pe

   def kinetic(self):
      """Calculates the kinetic energy of the cell from the cell parameters"""

      ke = 0.0
      for i in range(3):
         for j in range(3):
            ke += self.__p[i, j]**2
      ke /= 2*self.w
      return ke
      
   def apply_pbc(self, atom):
      """Uses the minimum image convention to return a particle to the
         unit cell"""

      s=numpy.dot(self.ih,atom.q)
      for i in range(3):
         s[i] = s[i] - math.floor(s[i])
      new_pos = numpy.dot(self.h,s)
      return new_pos

   def minimum_distance(self, atom1, atom2):
      """Takes two atoms and tries to find the smallest distance between two 
         images. This is only rigorously accurate in the case of a cubic cell,
         but gives the correct results as long as the cut-off radius is defined
         as smaller than the smallest width between parallel faces."""

      s1 = numpy.dot(self.ih, atom1.q)
      s2 = numpy.dot(self.ih, atom2.q)
      s = numpy.empty(3)
      for i in range(3):
         s[i] = s1[i] - s2[i] - round(s1[i] - s2[i])
      r = numpy.dot(self.h, s)
      return r
      
def h2abc(h):
   """
   Returns a description of the cell in terms of the length of the 
   lattice vectors and the angles between them."""
   
   a=h[0,0]
   b=math.sqrt(h[0,1]**2+h[1,1]**2)
   c=math.sqrt(h[0,2]**2+h[1,2]**2+h[2,2]**2)
   gamma=math.acos(h[0,1]/b)
   beta=math.acos(h[0,2]/c)
   alpha = math.acos(numpy.dot(h[:,1], h[:,2])/(b*c))

   return a, b, c, alpha, beta, gamma

def abc2h(a, b, c, alpha, beta, gamma):
   """
   Returns a cell matrix given a description in terms of the vector lengths and
   the angles in between"""

   h = numpy.zeros((3,3) ,float)
   h[0,0] = a
   h[0,1] = b * math.cos(gamma)
   h[0,2] = c * math.cos(beta)
   h[1,1] = b * math.sin(gamma)
   h[1,2] = (b*c*math.cos(alpha) - h[0,1]*h[0,2]) / h[1,1]
   h[2,2] = math.sqrt(c**2 - h[0,2]**2 - h[1,2]**2)
   return h
   
