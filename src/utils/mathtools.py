import numpy
import math

def matrix_exp(M, ntaylor=8, nsquare=8):
   n=M.shape[1]
   # computes the coefficients of a Taylor expansion of the exponential
   tc=numpy.zeros(ntaylor+1)
   tc[0]=1.0
   for i in range(ntaylor): tc[i+1]=tc[i]/(i+1)
   
   # scale
   SM=numpy.copy(M)/2.0**nsquare
   
   # truncated Taylor exp.
   EM=numpy.identity(n,float)*tc[ntaylor]
   for i in range(ntaylor-1,-1,-1):
      EM=numpy.dot(SM,EM)
      EM+=numpy.identity(n)*tc[i]
   
   # squaring
   for i in range(nsquare): EM=numpy.dot(EM,EM)
   return EM

def stab_cholesky(M):
   n=M.shape[1]
   D=numpy.zeros(n,float)
   L=numpy.zeros(M.shape,float)
   for i in range(n):
      L[i,i]=1.
      D[i]=M[i,i]
      for j in range(i):
         L[i,j]=M[i,j]
         for k in range(j): L[i,j]-=L[i,k]*L[j,k]*D[k]
         if (not D[j]==0.0): L[i,j]=L[i,j]/D[j]

      for k in range(i): D[i]-=L[i,k]*L[i,k]*D[k]

   S=numpy.zeros(M.shape,float)
   for i in range(n): 
      if (D[i]>0): D[i]=math.sqrt(D[i])
      else: D[i]=0.0
      for j in range(i+1):
         S[i,j]+=L[i,j]*D[j]
   return S
   
def invert_ut3x3(h):
   """Inverts a 3*3 (upper-triangular) cell matrix"""
   ih = numpy.zeros((3,3), float)
   for i in range(3):
      ih[i,i] = 1.0/h[i,i]
   ih[0,1] = -ih[0,0]*h[0,1]*ih[1,1]
   ih[1,2] = -ih[1,1]*h[1,2]*ih[2,2]
   ih[0,2] = -ih[1,2]*h[0,1]*ih[0,0]-ih[0,0]*h[0,2]*ih[2,2]
   return ih

def eigensystem_ut3x3(p):
   """Finds the eigenvector matrix of a 3*3 upper-triangular matrix"""
   eigp = numpy.zeros((3,3), float)
   eigvals = numpy.zeros(3, float)

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
