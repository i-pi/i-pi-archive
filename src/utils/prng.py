"""Contains the classes used to generate pseudo-random numbers.

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
   """Class to interface with the standard pseudo-random number generator.

   Initialises the standard numpy pseudo-random number generator from a seed
   at the beginning of the simulation, and keeps track of the state so that
   it can be output to the checkpoint files throughout the simulation.

   Attributes:
      rng: The random number generator to be used.
      seed: The seed number to start the generator.
      state: A tuple of five objects giving the current state of the random
         number generator. The first is the type of random number generator, 
         here 'MT19937', the second is an array of 624 integers, the third
         is the current position in the array that is being read from, the 
         fourth gives whether it has a gaussian random number stored, and 
         the fifth is this stored Gaussian random number, or else the last
         Gaussian random number returned.
   """

   def __init__(self, seed=12345, state=None):
      """Initialises Random.

      Args:
         seed: An optional seed giving an integer to initialise the state with.
         state: An optional state tuple to initialise the state with.
      """

      self.rng = np.random.mtrand.RandomState(seed=seed)
      self.seed = seed
      if state is None:   
         self.rng.seed(seed)        
      else:
         self.state = state

   @property
   def state(self):
      """Interface to the standard get_state() function."""

      return self.rng.get_state()

   @state.setter
   def state(self, value):
      """Interface to the standard set_state() function.

      Should only be used with states generated from another similar random
      number generator, such as one from a previous run.
      """

      return self.rng.set_state(value)

   @property
   def u(self):
      """Interface to the standard random_sample() function.

      Returns:
         A pseudo-random number from a uniform distribution from 0-1.
      """

      return self.rng.random_sample()

   @property
   def g(self):
      """Interface to the standard standard_normal() function.

      Returns:
         A pseudo-random number from a normal Gaussian distribution.
      """

      return self.rng.standard_normal()
   
   def gamma(self, k, theta=1.0):
      """Interface to the standard gamma() function.

      Args:
         k: Shape parameter for the gamma distribution.
         theta: Mean of the distribution.

      Returns:
         A random number from a gamma distribution with a shape k and a 
         mean value theta.
      """

      return self.rng.gamma(k,theta)

   def gvec(self, shape):
      """Interface to the standard_normal array function.

      Args:
         shape: The shape of the array to be returned.

      Returns:
         An array with the required shape where each element is taken from
         a normal Gaussian distribution.
      """

      return self.rng.standard_normal(shape)


class RestartRandom(Restart):
   """Random restart class.

   Handles generating the appropriate random number class from the xml 
   input file, and generating the xml checkpoint tags and data from an 
   instance of the object.

   Attributes:
      seed: An optional integer giving a seed to initialise the random number
         generator from. Defaults to 123456.
      state: An optional array giving the state of the random number generator.
         Defaults to an empty array.
      has_gauss: An optional integer giving whether there is a stored 
         Gaussian number or not. Defaults to 0.
      gauss: An optional float giving the stored Gaussian number. Defaults to
         0.0.
      set_pos: An optional integer giving the position in the state array
         that is being read from. Defaults to 0.
   """

   fields = {"seed" : (RestartValue, (int, 123456)), 
             "state": (RestartArray, (np.uint, np.zeros(0, np.uint ))),
             "has_gauss": (RestartValue, (int, 0)),  
             "gauss": (RestartValue, (float,  0.00 )),
             "set_pos": (RestartValue, (int, 0))
            }

   def store(self, prng):
      """Takes a random number instance and stores a minimal 
      representation of it.

      Args:
         prng: A random number object from which to initialise from.
      """

      self.seed.store(prng.seed)
      gstate = prng.state
      self.state.store(gstate[1])
      self.set_pos.store(gstate[2])
      self.has_gauss.store(gstate[3])
      self.gauss.store(gstate[4])

   def fetch(self):
      """Creates a random number object.

      Returns:
         An random number object of the appropriate type and with the 
         appropriate properties given the attributes of the RestartRandom
         object.
      """

      state = self.state.fetch()
      if state.shape == (0,):
         return Random(seed=self.seed.fetch())
      else:
         return Random(state=('MT19937', self.state.fetch(), self.set_pos.fetch(), self.has_gauss.fetch(), self.gauss.fetch() ))

