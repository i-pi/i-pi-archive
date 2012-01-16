import numpy as np
import math
from restart import Restart, RestartValue, RestartArray
import pdb

class MRG32k3a:
   # a stripped-down version of L'Ecuyer MRG32k3a. Good efficiency and RN quality,
   # it is possible to do arbitrary jumps in the stream (but this is not implemented here)
   anti=False
   m1     =  4294967087.0        
   m2     =  4294944443.0        
   a12    =    1403580.0         
   a13n   =     810728.0         
   a21    =     527612.0         
   a23n   =    1370589.0         
   norm   =  2.328306549295727688e-10
   fact   =  5.9604644775390625e-8

   @staticmethod
   def U01(status):
      p1 = MRG32k3a.a12*status[0,1] - MRG32k3a.a13n*status[0,0]
      k = int(p1/MRG32k3a.m1)

      p1 = p1 - k*MRG32k3a.m1
      if (p1<0.0):
         p1 += MRG32k3a.m1
      status[0,0] = status[0,1]
      status[0,1] = status[0,2]
      status[0,2] = p1

      p2 = MRG32k3a.a21*status[1,2] - MRG32k3a.a23n*status[1,0]
      k = int(p2/MRG32k3a.m2)
      p2 = p2-k*MRG32k3a.m2
      if (p2<0.0):
         p2 += MRG32k3a.m2
      status[1,0] = status[1,1]
      status[1,1] = status[1,2]
      status[1,2] = p2
      
      if (p1>p2):
         ui = (p1 - p2)*MRG32k3a.norm
      else:
         u = (p1 - p2 + MRG32k3a.m1)*MRG32k3a.norm

      if (MRG32k3a.anti):
         return 1 - u
      else:
         return u

class Random(object):
   def __init__(self, seed=12345, state=None):
      self.rng = np.random.mtrand.RandomState(seed=seed)
      self.seed = seed
      if state is None:   
         self.rng.seed(seed)        
      else:
         self.state = state

   @property
   def state(self):
      return self.rng.get_state()

   @state.setter
   def state(self, value):
      return self.rng.set_state(value)

   @property
   def u(self):
      return self.rng.random_sample()

   @property
   def g(self):
      return self.rng.standard_normal()
   
   def gamma(self, k, theta=1.0):
      return self.rng.gamma(k,theta)
      
#      # generates a Gaussian variate with zero mean and unit varianceq
#      # uses a ratio-of-uniforms rather than the usual box-mueller method
#      # as this is (often) faster, and avoids the hassle of storing one of the samples
#      while True:
#         u=self.u; v=self.u 
#         v=1.7156*(v-0.5)  
#         x=u-0.449871
#         y=abs(v)+0.386595
#         q=x*x+y*(0.19600*y-0.25472*x)
#         if (q < 0.27597): break
#         if (q > 0.27846): continue
#         if (v*v<-4.0*math.log(u)*u*u): break

#      return v/u

   def gvec(self, shape):
      return self.rng.standard_normal(shape)
#      v=np.zeros(shape,float)
#      sz=v.size; fv=v.flat
#      for i in range(sz):
#         fv[i]=self.g
#      return v


class RestartRandom(Restart):
   fields = {"seed" : (RestartValue, (int, 123456)), 
             "state": (RestartArray, (np.uint, np.zeros(0, np.uint ))),
             "has_gauss": (RestartValue, (int, 0)),  
             "gauss": (RestartValue, (float,  0.00 )),
             "set_pos": (RestartValue, (int, 0))
            }

   def store(self, prng):
      self.seed.store(prng.seed)
      gstate = prng.state
      self.state.store(gstate[1])
      self.set_pos.store(gstate[2])
      self.has_gauss.store(gstate[3])
      self.gauss.store(gstate[4])

   def fetch(self):
      state = self.state.fetch()
      if state.shape == (0,):
         return Random(seed=self.seed.fetch())
      else:
         return Random(state=('MT19937', self.state.fetch(), self.set_pos.fetch(), self.has_gauss.fetch(), self.gauss.fetch() ))

