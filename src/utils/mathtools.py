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
