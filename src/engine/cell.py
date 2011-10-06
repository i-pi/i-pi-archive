import numpy
import math

def compute_ih(h) :
   ih = numpy.zeros((3,3), float)
   for i in range(3):
      ih[i,i] = 1.0/h[i,i]
   ih[0,1] = -h[0,1]/(h[0,0] * h[1,1])
   ih[1,2] = -h[1,2]/(h[1,1] * h[2,2])
   ih[0,2] = (h[0,1]*h[1,0] - h[0,2]*h[1,1]) / (h[0,0]*h[1,1]*h[2,2])
   return h


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
   
   def __init__(self, cell):
      self.__h = numpy.zeros((3,3), float)
      a, b, c, alpha, beta, gamma = cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]
      self.abc2h(a, b, c, alpha, beta, gamma)
      self.__p = numpy.zeros((3,3) ,float)
      self.__taint_ih = True
      self.w = 1.0

   def __str__(self):
      return "    Unit vectors: %s %s %s \n    Momenta: %s %s %s \n    w = %s, volume = %s" % (self.h[:,0], self.h[:,1], self.h[:,2], self.p[:,0], self.p[:,1], self.p[:,2], self.w, self.volume())
      
   def volume(self):
      #return numpy.inner(self.h[0:3,0], numpy.cross(self.h[0:3,1], self.h[0:3,2])) # general matrix
      return self.__h[0,0]*self.__h[1,1]*self.__h[2,2]   # upper-triangular matrix

   def kinetic(self):
      ke = 0.0
      for i in range(3):
         for j in range(3):
            ke += self.__p[i, j]**2
      ke /= 2*self.w
      return ke
      
   def apply_pbc(self, atom):
      s=numpy.dot(self.ih,atom.q)
      for i in range(0,3):
         s[i] = s[i] - round(s[i])
      print s
      atom = numpy.dot(self.h,s)
      print atom
      
   def h2abc(self):
      """
      Returns a description of the cell in terms of the length of the 
      lattice vectors and the angles between them."""
      
      a=self.__h[0,0]; b=math.sqrt(self.__h[0,1]**2+self.__h[1,1]**2);  c=math.sqrt(self.__h[0,2]**2+self.__h[1,2]**2+self.__h[2,2]**2);
      gamma=math.acos(self.__h[0,1]/b); beta=math.acos(self.__h[0,2]/c); 
      alpha = math.acos(numpy.dot(self.__h[:,1], self.__h[:,2])/(b*c))

      return a, b, c, alpha, beta, gamma

   def abc2h(self, a, b, c, alpha, beta, gamma):
      self.__h[0,0] = a
      self.__h[0,1] = b * math.cos(gamma)
      self.__h[0,2] = c * math.cos(beta)
      self.__h[1,1] = b * math.sin(gamma)
      self.__h[1,2] = (b*c*math.cos(alpha) - self.__h[0,1]*self.__h[0,2]) / self.__h[1,1]
      self.__h[2,2] = math.sqrt(c**2 - self.__h[0,2]**2 - self.__h[1,2]**2)
