import numpy as np
import math

def matrix_exp(M, ntaylor=8, nsquare=8):
   n=M.shape[1]
   # computes the coefficients of a Taylor expansion of the exponential
   tc=np.zeros(ntaylor+1)
   tc[0]=1.0
   for i in range(ntaylor): tc[i+1]=tc[i]/(i+1)
   
   # scale
   SM=np.copy(M)/2.0**nsquare
   
   # truncated Taylor exp.
   EM=np.identity(n,float)*tc[ntaylor]
   for i in range(ntaylor-1,-1,-1):
      EM=np.dot(SM,EM)
      EM+=np.identity(n)*tc[i]
   
   # squaring
   for i in range(nsquare): EM=np.dot(EM,EM)
   return EM

def stab_cholesky(M):
   n=M.shape[1]
   D=np.zeros(n,float)
   L=np.zeros(M.shape,float)
   for i in range(n):
      L[i,i]=1.
      D[i]=M[i,i]
      for j in range(i):
         L[i,j]=M[i,j]
         for k in range(j): L[i,j]-=L[i,k]*L[j,k]*D[k]
         if (not D[j]==0.0): L[i,j]=L[i,j]/D[j]

      for k in range(i): D[i]-=L[i,k]*L[i,k]*D[k]

   S=np.zeros(M.shape,float)
   for i in range(n): 
      if (D[i]>0): D[i]=math.sqrt(D[i])
      else: D[i]=0.0
      for j in range(i+1):
         S[i,j]+=L[i,j]*D[j]
   return S

def h2abc(h):
   """Returns a description of the cell in terms of the length of the 
      lattice vectors and the angles between them."""
   
   a=float(h[0,0])
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
   
def invert_ut3x3(h):
   """Inverts a 3*3 (upper-triangular) cell matrix"""
   ih = np.zeros((3,3), float)
   for i in range(3):
      ih[i,i] = 1.0/h[i,i]
   ih[0,1] = -ih[0,0]*h[0,1]*ih[1,1]
   ih[1,2] = -ih[1,1]*h[1,2]*ih[2,2]
   ih[0,2] = -ih[1,2]*h[0,1]*ih[0,0]-ih[0,0]*h[0,2]*ih[2,2]
   return ih

def eigensystem_ut3x3(p):
   """Finds the eigenvector matrix of a 3*3 upper-triangular matrix"""
   eigp = np.zeros((3,3), float)
   eigvals = np.zeros(3, float)

   for i in range(3):
      eigp[i,i] = 1
   eigp[0,1] = -p[0,1]/(p[0,0] - p[1,1])
   eigp[1,2] = -p[1,2]/(p[1,1] - p[2,2])
   eigp[0,2] = -(p[0,1]*p[1,2] - p[0,2]*p[1,1] + p[0,2]*p[2,2])/((p[0,0] - p[2,2])*(p[2,2] - p[1,1]))

   for i in range(3):
      eigvals[i] = p[i,i]
   return eigp, eigvals

def det_ut3x3(h):
   """Calculates the volume of the unit cell, assuming an upper-triangular
      unit vector matrix"""
   return h[0,0]*h[1,1]*h[2,2]
   
MINSERIES=1e-8
def exp_ut3x3(h):
   """Computes the matrix exponential for a 3x3 upper-triangular matrix"""
   eh=np.zeros((3,3), float)
   e00=math.exp(h[0,0]);    e11=math.exp(h[1,1]);    e22=math.exp(h[2,2])
   eh[0,0]=e00;    eh[1,1]=e11;    eh[2,2]=e22; 
   
   print "minseries is ", MINSERIES
   if (abs((h[0,0]-h[1,1])/h[0,0])>MINSERIES): 
      r01=(e00-e11)/(h[0,0]-h[1,1])
   else:
      r01=e00*(1+(h[0,0]-h[1,1])*(0.5+(h[0,0]-h[1,1])/6.0))
   if (abs((h[1,1]-h[2,2])/h[1,1])>MINSERIES): 
      r12=(e11-e22)/(h[1,1]-h[2,2])   
   else:
      r12=e11*(1+(h[1,1]-h[2,2])*(0.5+(h[1,1]-h[2,2])/6.0))
   if (abs((h[2,2]-h[0,0])/h[2,2])>MINSERIES): 
      r02=(e22-e00)/(h[2,2]-h[0,0])   
   else:
      r02=e22*(1+(h[2,2]-h[0,0])*(0.5+(h[2,2]-h[0,0])/6.0))
      
   eh[0,1]=h[0,1]*r01
   eh[1,2]=h[1,2]*r12

   eh[0,2]=h[0,2]*r02
   if (abs((h[2,2]-h[0,0])/h[2,2])>MINSERIES):    
      eh[0,2]+=h[0,1]*h[0,2]*(r01-r12)/(h[0,0]-h[2,2])
   elif (abs((h[1,1]-h[0,0])/h[1,1])>MINSERIES):    
      eh[0,2]+=h[0,1]*h[0,2]*(r12-r02)/(h[1,1]-h[0,0])
   elif (abs((h[1,1]-h[2,2])/h[1,1])>MINSERIES):    
      eh[0,2]+=h[0,1]*h[0,2]*(r02-r01)/(h[2,2]-h[1,1])
   else:
      eh[0,2]+=h[0,1]*h[0,2]*e00/24.0*(12.0+4*(h[1,1]+h[2,2]-2*h[0,0])+(h[1,1]-h[0,0])*(h[2,2]-h[0,0]))
      
   return eh  
