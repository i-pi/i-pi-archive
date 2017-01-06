"""Holds the algorithms to perform replica exchange.

Algorithms implemented by Robert Meissner and Riccardo Petraglia, 2016
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2016 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np
import time

from ipi.engine.smotion import Smotion
from ipi.utils.depend import *
from ipi.utils.softexit import softexit
from ipi.utils.messages import verbosity, info
from ipi.utils.units import Constants


__all__ = ['ReplicaExchange']


# TODO: Do not shout :-)
#       (1) Exchange of Hamiltonians is missing

class ReplicaExchange(Smotion):
    """Replica exchange routine.

    Attributes:
        every: on which steps REMD should be performed
        exchange:
            temperature: activate temperature replica exchange
            hamiltonian: activate hamiltonian replica exchange
            bias: activate hamiltonian replica exchange ***not yet implemented
    """

    def __init__(self, stride=1.0):
        """Initialises REMD.

        Args:

        """

        super(ReplicaExchange, self).__init__()

        # replica exchange options
        self.stride = stride

    def bind(self, syslist, prng):

        super(ReplicaExchange,self).bind(syslist, prng)

    def step(self, step=None):
        """Tries to exchange replica."""
        
        if self.stride <= 0.0: return
        
        syspot  = [ s.forces.pot for s in self.syslist ]
        # spring potential in a form that can be easily used further down (no temperature included!)
        syspath = [ s.beads.vpath/Constants.hbar**2 for s in self.syslist ]

        info("\nTrying to exchange replicas on STEP %d" % step, verbosity.debug)

        self.ptime = self.ttime = 0
        self.qtime = -time.time()

        for i in range(len(self.syslist)):
           for j in range(i):
              if (1.0/self.stride < self.prng.u) : continue  # tries a swap with probability 1/stride
              betai = 1.0/(Constants.kb*self.syslist[i].ensemble.temp*self.syslist[i].beads.nbeads); 
              betaj = 1.0/(Constants.kb*self.syslist[j].ensemble.temp*self.syslist[j].beads.nbeads);


              pxc = np.exp(
                (betai * syspot[i] + syspath[i]/betai +
                 betaj * syspot[j] + syspath[j]/betaj) -
                (betai * syspot[j] + syspath[j]/betai +
                 betaj * syspot[i] + syspath[i]/betaj)
                )

              if (pxc > self.prng.u): # really does the exchange
                 info(" @ PT:  SWAPPING replicas % 5d and % 5d." % (i,j), verbosity.low)
                 info(" @ PT:  % 5d and % 5d." % (self.syslist[i].ensemble.temp,self.syslist[j].ensemble.temp), verbosity.low)
                 
                 #tempi = copy(self.syslist[i].ensemble.temp)
                 #self.syslist[i].ensemble.temp = copy(self.syslist[j].ensemble.temp)
                 # velocities have to be adjusted according to the new temperature
                 


        self.qtime += time.time()
