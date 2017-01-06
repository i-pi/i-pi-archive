"""Holds the algorithms to perform replica exchange.

Algorithms implemented by Robert Meissner and Riccardo Petraglia, 2016
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2016 i-PI developers
# See the "licenses" directory for full license information.


import numpy as np
import time

from ipi.engine.smootion import Smootion
from ipi.utils.depend import *
from ipi.utils.softexit import softexit
from ipi.utils.messages import verbosity, info


__all__ = ['REMover']


# TODO: Do not shout :-)
# NOTE: REMD needs to be done properly

class REMover(Smootion):
    """Replica exchange routine.

    Attributes:
        every: on which steps REMD should be performed
        exchange:
            temperature: activate temperature replica exchange
            hamiltonian: activate hamiltonian replica exchange
            bias: activate hamiltonian replica exchange ***not yet implemented
    """

    def __init__(self,
                 exchange="temperature",
                 every=1.0):
        """Initialises REMD.

        Args:

        """

        super(REMD, self).__init__()

        # REMD options
        self.exchange = exchange
        self.every = every

    def bind(self, ens, beads, nm, cell, bforce, prng):

        super(REMD,self).bind(ens, beads, nm, cell, bforce, prng)

    def step(self, step=None):
        """Does one simulation time step."""

        info("\nTrying to exchange replicas on STEP %d" % step, verbosity.debug)

        self.ptime = self.ttime = 0
        self.qtime = -time.time()

        # <-- insert code here -->

        self.qtime += time.time()
