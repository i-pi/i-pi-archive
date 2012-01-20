"""Contains the classes used to keep track of the pseudo-random numbers.

Allows the user to specify a seed for the random number generator.
These are used in initialising the velocities and in stochastic thermostats.
The state of the random number generator is kept track of, so that the if the
simulation is restarted from a checkpoint, we will see the same dynamics as if 
it had not been stopped.

Classes:
   MRG32k3a: A class implementing a basic L'Ecuyer MRG32k3a random number
      generator.
   Random: An interface between the numpy.random module and the user.
   RestartRandom: Deals with creating the Random object from a file, and 
      writing the checkpoints.
"""

__all__ = ['MRG32k3a', 'Random', 'RestartRandom']

import numpy as np
import math
from restart import Restart, RestartValue, RestartArray

class MRG32k3a:
   """Class implementing a simple L'Ecuyer MRG32k3a.

   A simple uniform pseudo-random number generator class.
   This type of random number generator has the ability to jump an arbitrary
   distance in the stream, but this has not yet been implemented.

   Attributes:
      anti: A boolean giving whether the distribution is from 0 to 1 or from
         1 to 0.
      m1:
   """
   anti = False
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

   def gvec(self, shape):
      return self.rng.standard_normal(shape)


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

