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
      system_temp: The actual temperatures of the various systems
      system_vf: The bias and force scaling of the various system
      system_v: The bias of the various system
      system_f: The force of the various systems
   """

   def __init__(self, tlist=None, ilist=None, stride=0.0, wteml=None, wtesl=None, wtegl=None):
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

      if wteml is None:  wteml = []
      if wtesl is None:  wtesl = []
      if wtegl is None:  wtegl = []
      if len(wteml)> 0 and (len(wtegl)!= len(wteml) or len(wtesl)!=len(wteml) or len(wteml)!=len(tlist)):
         raise ValueError("WTE parameters list must all be provided, and match temperature list size.")
      self.wte_means = np.asarray(wteml,float).copy()
      self.wte_sigmas = np.asarray(wtesl,float).copy()
      self.wte_gammas = np.asarray(wtegl,float).copy()

      dset(self,"system_temp",depend_array(name="system_temp", value=np.asarray(tlist).copy(), func=self.get_stemp,
                  dependencies=[dget(self,"temp_index")]))

      self.parafile = None
      self.wtefile = None

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

      def make_biasgetter(k):
         return lambda: self.system_v[k]

      if (len(self.wte_means)>0):
         dset(self,"system_vf",depend_array(name="system_vf", value=np.zeros((len(self.wte_means),2), float)
               ) )
         self.system_v = self.system_vf[:,0]
         self.system_f = self.system_vf[:,1]
         for s in self.slist:
            dget(s.forces, "pot").add_dependant(dget(self,"system_vf"))
         dget(self,"temp_index").add_dependant(dget(self,"system_vf"))
         dget(self,"system_vf")._func=self.get_wtevf
         dget(self,"system_f")._func=self.get_wtevf
         dget(self,"system_v")._func=self.get_wtevf
      else:
         dset(self,"system_vf",depend_array(name="system_vf", value=np.zeros((len(self.temp_list),2), float)) )
         self.system_v = self.system_vf[:,0]
         self.system_f = self.system_vf[:,1]

      isys=0
      for s in self.slist:
         dget(s.ensemble,"bias").add_dependency(dget(self,"system_v"))
         dget(s.ensemble,"bias")._func = make_biasgetter(isys)
         dget(s.ensemble,"temp").add_dependency(dget(self,"system_temp"))
         dget(s.ensemble,"temp")._func = make_tempgetter(isys)
         isys+=1

      self.parafile=open("PARATEMP", "a")
      self.wtefile=None
      if len(self.wte_means)>0:
         self.wtefile=open("PARAWTE", "a")

   def get_stemp(self):
      """ Returns the temperatures of the various systems. """

      return np.asarray([ self.temp_list[self.temp_index[i]] for i in range(len(self.temp_list))])

   def wtevf(self, i, j):
      """ Gets the actual bias of system i subject to the j-th WTE bias. """

      betaj=1.0/(Constants.kb*self.temp_list[j])
      s = self.slist[i]
      ui = s.forces.pot
      vij = (1-1.0/self.wte_gammas[j])*self.wte_gammas[j]/betaj * np.exp(
          -0.5/(self.wte_gammas[j])*
            ((ui/s.beads.nbeads-self.wte_means[j])/self.wte_sigmas[j])**2
          )
      fij = -vij*((ui/s.beads.nbeads-self.wte_means[j])/
             (self.wte_sigmas[j]**2*self.wte_gammas[j]*s.beads.nbeads))
      return (vij, fij)

   def get_wtevf(self):

      vlist=np.zeros((len(self.temp_list),2), float)
      for i in range(len(self.temp_list)):
         vlist[i] = self.wtevf(i, self.temp_index[i])
      return vlist

   def wtestep(self, step=-1):

      if len(self.wte_means) == 0 or  len(self.wte_sigmas) == 0: return

      # note that this loops over the TEMPERATURES, not over the SYSTEMS
      for i in range(len(self.temp_list)):
         s = self.slist[i]
         ui = s.forces.pot
         f = depstrip(s.forces.f)*self.system_f[i]
         s.beads.p += f*s.ensemble.dt*0.5



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

            if len(self.wte_gammas)>0:
               vii, dummy = self.wtevf(i, self.temp_index[i])
               vij, dummy = self.wtevf(i, self.temp_index[j])
               vjj, dummy = self.wtevf(j, self.temp_index[j])
               vji, dummy = self.wtevf(j, self.temp_index[i])
            else:
               vii = vij = vji = vjj = 0
               
            pxc = np.exp(
              (betai * syspot[i] + syspath[i]/betai + betai * vii +
               betaj * syspot[j] + syspath[j]/betaj + betaj * vjj) -
              (betai * syspot[j] + syspath[j]/betai + betai * vji +
               betaj * syspot[i] + syspath[i]/betaj + betaj * vij)
              )
            print i, j, pxc
            if (pxc > self.prng.u): # really does the exchange
               info(" @ PT:  SWAPPING replicas % 5d and % 5d." % (i,j), verbosity.low)
               # adjusts the conserved quantities
               self.slist[i].ensemble.eens += self.slist[i].nm.kin *(1.0- (betai/betaj))
               self.slist[j].ensemble.eens += self.slist[j].nm.kin *(1.0- (betaj/betai))

               # adjusts the momenta
               self.slist[i].beads.p *= np.sqrt(betai/betaj)
               self.slist[j].beads.p *= np.sqrt(betaj/betai)
               # if there are GLE thermostats around, we must also rescale the s momenta!
               # should also check the barostat thermostat, but we don't do NPT replica exchange yet so whatever.
               if hasattr(self.slist[i].ensemble.thermostat,"s"):
                  self.slist[i].ensemble.thermostat.s *= np.sqrt(betai/betaj)
                  self.slist[j].ensemble.thermostat.s *= np.sqrt(betaj/betai)

               # adjusts the contribution from the WTE bias
               self.slist[i].ensemble.eens += vii-vij
               self.slist[j].ensemble.eens += vjj-vji

               swp=self.temp_index[j];  self.temp_index[j]=self.temp_index[i];  self.temp_index[i]=swp


   def softexit(self):
      if not self.parafile is None:
         self.parafile.close()
      self.parafile = None


