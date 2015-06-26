"""Contains a helper class for parallel tempering simulations.

Manages options, restarts and the actual exchange process for parallel
tempering simulations.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np

from ipi.utils.depend import *
from ipi.utils.messages import info, verbosity
from ipi.utils.units import Constants
from ipi.engine.thermostats import *


__all__ = ['ParaTemp']


class ParaTemp(dobject):
   """Helper class for parallel tempering simulations.

   Contains options and routines for performing parallel tempering runs.
   Contrary to the typical PT implementation, in which the positions are
   exchanged, here we exchange temperatures. This makes post-processing
   a bit more cumbersome, but is much better considering the interaction
   with the first-principles code -- that would not like to see wildly
   different positions after the exchange.

   Attributes:
      stride: How often to do exchanges (on average)
      temp_list: List of reference temperatures
      temp_index: Keeps track of the temperature of each system
      slist: List of systems corresponding to the various replicas

   Depend objects:
      system_temp: The actual temperatures of the various systems
   """

   def __init__(self, tlist=None, ilist=None, stride=0.0):
      """Initializes ParaTemp object.

      Parameters:
         tlist: List of temperatures to be copied in temp_list
         stride: Exchange stride, to be copied in stride
      """

      self.stride = stride
      if tlist is None:
         tlist = []
      if ilist is None or len(ilist)==0:
         ilist = range(len(tlist))   # defaults to normal ordering of temperatures

      if len(ilist)!= len(tlist):
         raise ValueError("Temperature list and index list have mismatching sizes.")

      dset(self, "temp_index",depend_array(name="temp_index", value=np.asarray(ilist, int).copy()) )
      self.temp_list = np.asarray(tlist, float).copy()

      dset(self,"system_temp",depend_array(name="system_temp", value=np.asarray(tlist).copy(), func=self.get_stemp,
                  dependencies=[dget(self,"temp_index")]))

      self.parafile = None

   def bind(self, slist, prng):
      """Wires up the PT setup, by connecting the replicas to the temperature list.

      """

      self.prng=prng
      self.slist = slist
      if len(slist)!=len(self.temp_list):
         raise ValueError("Number of systems does not match size of list of temperatures.")

      # now makes sure that the temperatures of the ensembles of the systems are
      # piped to temp_list
      def make_tempgetter(k):
         return lambda: self.system_temp[k]

      isys=0
      for s in self.slist:
         dget(s.ensemble,"temp").add_dependency(dget(self,"system_temp"))
         dget(s.ensemble,"temp")._func = make_tempgetter(isys)
         isys+=1

      self.parafile=open("PARATEMP", "a")

   def get_stemp(self):
      """ Returns the temperatures of the various systems. """

      return np.asarray([ self.temp_list[self.temp_index[i]] for i in range(len(self.temp_list))])

   def swap(self, step=-1):
      """ Tries a PT swap move. """

      if self.stride <= 0.0: return

      syspot  = [ s.forces.pot for s in self.slist ]
      # spring potential in a form that can be easily used further down (no temperature included!)
      syspath = [ s.beads.vpath/Constants.hbar**2 for s in self.slist ]

      # tries exchanges. note that we don't just exchange neighbouring replicas but try all pairs
      # 1. since this can in principle speed up diffusion by allowing "double jumps"
      # 2. since temp_list is NOT sorted, and so neighbouring temp_list could be actually far off
      #    and re-sorting would be more bookkeeping I can stand.
      for i in range(len(self.slist)):
         for j in range(i):
            if (1.0/self.stride < self.prng.u) : continue  # tries a swap with probability 1/stride
            # ALL SYSTEMS ARE EXPECTED TO HAVE SAME N OF BEADS!
            betai = 1.0/(Constants.kb*self.system_temp[i]*self.slist[i].beads.nbeads); # exchanges are being done, so it is better to re-compute betai in the inner loop
            betaj = 1.0/(Constants.kb*self.system_temp[j]*self.slist[j].beads.nbeads);


            pxc = np.exp(
              (betai * syspot[i] + syspath[i]/betai +
               betaj * syspot[j] + syspath[j]/betaj) -
              (betai * syspot[j] + syspath[j]/betai +
               betaj * syspot[i] + syspath[i]/betaj)
              )

            if (pxc > self.prng.u): # really does the exchange
               info(" @ PT:  SWAPPING replicas % 5d and % 5d." % (i,j), verbosity.low)
               # adjusts the conserved quantities
               # change in kinetic energy
               self.slist[i].ensemble.eens += self.slist[i].nm.kin *(1.0- (betai/betaj))
               self.slist[j].ensemble.eens += self.slist[j].nm.kin *(1.0- (betaj/betai))
               # change in spring energy
               self.slist[i].ensemble.eens += syspath[i]*(1.0/betai**2- 1.0/betaj**2)
               self.slist[j].ensemble.eens += syspath[j]*(1.0/betaj**2- 1.0/betai**2)

               # adjusts the momenta
               self.slist[i].beads.p *= np.sqrt(betai/betaj)
               self.slist[j].beads.p *= np.sqrt(betaj/betai)

               # if there are GLE thermostats around, we must also rescale the s momenta!
               # should also check the barostat thermostat, but we don't do NPT replica exchange yet so whatever.
               if hasattr(self.slist[i].ensemble.thermostat,"s"):
                  self.slist[i].ensemble.thermostat.s *= np.sqrt(betai/betaj)
                  self.slist[j].ensemble.thermostat.s *= np.sqrt(betaj/betai)

               swp=self.temp_index[j];  self.temp_index[j]=self.temp_index[i];  self.temp_index[i]=swp


   def softexit(self):
      if not self.parafile is None:
         self.parafile.close()
      self.parafile = None
