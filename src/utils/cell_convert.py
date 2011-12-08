import numpy as np
import math

def h2abc(h):
   """Returns a description of the cell in terms of the length of the 
      lattice vectors and the angles between them."""
   
   a=h[0,0]
   b=math.sqrt(h[0,1]**2+h[1,1]**2)
   c=math.sqrt(h[0,2]**2+h[1,2]**2+h[2,2]**2)
   gamma=math.acos(h[0,1]/b)
   beta=math.acos(h[0,2]/c)
   alpha = math.acos(np.dot(h[:,1], h[:,2])/(b*c))

   return a, b, c, alpha, beta, gamma

def abc2h(a, b, c, alpha, beta, gamma):
   """Returns a cell matrix given a description in terms of the vector lengths
      and the angles in between"""

   h = np.zeros((3,3) ,float)
   h[0,0] = a
   h[0,1] = b * math.cos(gamma)
   h[0,2] = c * math.cos(beta)
   h[1,1] = b * math.sin(gamma)
   h[1,2] = (b*c*math.cos(alpha) - h[0,1]*h[0,2]) / h[1,1]
   h[2,2] = math.sqrt(c**2 - h[0,2]**2 - h[1,2]**2)
   return h
