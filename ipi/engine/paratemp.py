"""Contains a helper class for parallel tempering simulations.

Copyright (C) 2013, Michele Ceriotti

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http.//www.gnu.org/licenses/>.

Manages options, restarts and the actual exchange process for parallel
tempering simulations.

Classes:
   ParaTemp: Contains all parallel-tempering related functionalities
"""

__all__ = ['ParaTemp']

import numpy as np
from ipi.utils.depend import *
from ipi.utils.messages import info, verbosity
from ipi.utils.units import Constants
from ipi.engine.thermostats import *

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
      temp_replicas: The temperatures of the various replicas
   """

   def __init__(self, tlist=None, ilist=None, stride=1):
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

      self.temp_list = np.asarray(tlist, float).copy()
      self.temp_index = np.asarray(ilist, int).copy()

      dset(self,"temp_replicas",depend_array(name="temp_tlist", value=np.asarray(tlist).copy()))
      for i in range(len(tlist)):
         self.temp_replicas[i] = self.temp_list[self.temp_index[i]]
      self.parafile = None

   def bind(self, slist, prng):
      """Wirhttp://home.web.cern.ch/about/updates/2013/10/cern-host-advanced-materials-and-surfaces-workshopes up the PT setup, by connecting the replicas to the temperature list.

      """

      self.prng=prng
      self.slist = slist
      if len(slist)!=len(self.temp_list):
         raise ValueError("Number of systems does not match size of list of temperatures.")

      # now makes sure that the temperatures of the ensembles of the systems are
      # piped to temp_list
      def make_tempgetter(k):
         return lambda: self.temp_replicas[k]

      isys=0
      for s in self.slist:
         dget(s.ensemble,"temp").add_dependency(dget(self,"temp_replicas"))
         dget(s.ensemble,"temp")._func = make_tempgetter(isys)
         isys+=1

      self.parafile=open("PARATEMP", "a")

   def swap(self, step=-1):
      """ Tries a PT swap move. """

      # because of where this is in the loop, we must write out BEFORE doing the swaps.
      self.parafile.write("%10d" % (step+1))
      for i in self.temp_index:
         self.parafile.write(" %5d" %i)
      self.parafile.write("\n")
      
      sysham = [ s.forces.pot for s in self.slist ] #TODO: must check how to do with PIMD (I guess just add the springs)

      # tries exchanges. note that we don't just exchange neighbouring replicas but try all pairs
      # 1. since this can in principle speed up diffusion by allowing "double jumps"
      # 2. since temp_list is NOT sorted, and so neighbouring temp_list could be actually far off
      #    and re-sorting would be more bookkeeping I can stand.
      for i in range(len(self.slist)):
         for j in range(i-1):
            if (1.0/float(self.stride) < self.prng.u) : continue  # tries a swap with probability 1/stride
            pxc = np.exp(
              (1.0/(Constants.kb*self.temp_list[self.temp_index[j]]) -
              1/(Constants.kb*self.temp_list[self.temp_index[i]]) ) *
              (sysham[j]-sysham[i]))
            if (pxc > self.prng.u): # really does the exchange
               info(" @ PT:  SWAPPING replicas % 5d and % 5d." % (i,j), verbosity.medium)               
               self.slist[i].beads.p *= np.sqrt(self.temp_replicas[j]/self.temp_replicas[i])
               self.slist[j].beads.p *= np.sqrt(self.temp_replicas[i]/self.temp_replicas[j])
               # if there are GLE thermostats around, we must also rescale the s momenta!
               # should also check the barostat thermostat, but we don't do NPT replica exchange yet so whatever.
               if hasattr(self.slist[i].ensemble.thermostat,"s"):
                  self.slist[i].ensemble.thermostat.s *= np.sqrt(self.temp_replicas[j]/self.temp_replicas[i])
                  self.slist[j].ensemble.thermostat.s *= np.sqrt(self.temp_replicas[i]/self.temp_replicas[j])
               swp=self.temp_index[j];  self.temp_index[j]=self.temp_index[i];  self.temp_index[i]=swp;
               swp=self.temp_replicas[j];  self.temp_replicas[j]=self.temp_replicas[i];  self.temp_replicas[i]=swp;

   def softexit(self):
      if not self.parafile is None:
         self.parafile.close()
      self.parafile = None


