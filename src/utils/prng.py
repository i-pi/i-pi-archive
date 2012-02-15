"""Contains the classes used to generate pseudo-random numbers.

Allows the user to specify a seed for the random number generator.
These are used in initialising the velocities and in stochastic thermostats.
The state of the random number generator is kept track of, so that the if the
simulation is restarted from a checkpoint, we will see the same dynamics as if 
it had not been stopped.

Classes:
   Random: An interface between the numpy.random module and the user.
   RestartRandom: Deals with creating the Random object from a file, and 
      writing the checkpoints.
"""

__all__ = ['Random', 'RestartRandom']

import numpy as np
import math
from restart import Restart, RestartValue, RestartArray

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

   def get_state(self):
      """Interface to the standard get_state() function."""

      return self.rng.get_state()

   def set_state(self, value):
      """Interface to the standard set_state() function.

      Should only be used with states generated from another similar random
      number generator, such as one from a previous run.
      """

      return self.rng.set_state(value)

   state=property(get_state, set_state)

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

