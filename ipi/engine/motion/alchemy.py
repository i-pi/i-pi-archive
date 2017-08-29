"""Contains the classes that deal with the different dynamics required in
different types of ensembles.

Holds the algorithms required for normal mode propagators, and the objects to
do the constant temperature and pressure algorithms. Also calculates the
appropriate conserved energy quantity for the ensemble of choice.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import time

import numpy as np

from ipi.engine.motion import Motion
from ipi.utils.depend import *


#__all__ = ['AlchemyMC']

class AlchemyMC(Motion):
    """Monte Carlo alchemical exchanges class.

    Attributes:
        names of the isotopes for exchanges

    Depend objects:
        The spring potential energy, names of the atoms.
    """

    def __init__(self, fixcom=False, fixatoms=None, mode=None, names=[], nmc=1):
        """Initialises a "alchemical exchange" motion object.

        Args: 
            names : A list of isotopes
            nmc : frequency of doing exchanges
            
        """

        super(AlchemyMC, self).__init__(fixcom=fixcom, fixatoms=fixatoms)

        self.atomtype = names
        self.nummc = nmc      

    def bind(self, beads, ens, prng):
        """Binds ensemble beads, cell, bforce, and prng to the dynamics.

        This takes a beads object, a cell object, a forcefield object and a
        random number generator object and makes them members of the ensemble.
        It also then creates the objects that will hold the data needed in the
        ensemble algorithms and the dependency network. Note that the conserved
        quantity is defined in the init, but as each ensemble has a different
        conserved quantity the dependencies are defined in bind.

        Args:
            beads: The beads object from whcih the bead positions are taken.
            prng: The random number generator object which controls random number
                generation.
        """

        super(AlchemyMC, self).bind(beads, ens, prng)
        

    def AXlist(self, atomtype):
        """This compile a list of atoms ready for exchanges.
        """
            
        # selects the types of atoms for exchange
        atomexchangelist = []
        for i in range(self.beads.natoms):
            for j in atomtype:
                if (self.beads.names[i] == j): 
                    atomexchangelist.append(i)
                            
        return atomexchangelist

    def step(self, step=None):

        """Does one round of alchemical exchanges."""
        # record the spring energy (divided by mass) for each atom in the exchange chain
        q = depstrip(self.beads.q)
        nb = self.beads.nbeads
        na3 = self.beads.natoms * 3
        axlist = self.AXlist(self.atomtype)
        lenlist = len(axlist)
        atomspring = np.zeros(lenlist)
        i = 0
        for atomnum in axlist:
            spr = 0.0
            for b in range(1,nb):
                for j in range(na3,na3+3):
                    spr += (q[b,j]-q[b-1,j])**2
            for j in range(na3,na3+3):
                spr += (q[nb-1,j]-q[0,j])**2
            # no mass here
            spr *= 0.5*self.nm.omegan2
            atomspring[i] = spr
            i += 1
            
        # do the exchange
        betaP = 1.0/(Constants.kb*self.ensemble.temp*nb)
        for i in range(lenlist):
            for j in range(i):
                # no exchange for same type of atoms
                if (self.beads.names[axlist[i]] == self.beads.names[axlist[j]]) : continue
                # energy increase due to the swap
                difspring = (atomspring[i]-atomspring[j])*(self.beads.m[axlist[j]]-self.beads.m[axlist[i]])
                pexchange = np.exp(-betaP*difspring)
                #print 'exchange probablity: %10.5e  n. exchanges this far: %5d' % ( pexchange, nexch )
                
                # attemps the exchange
                if (pexchange > self.prng.u):
                    nexch += 1
                    # swap names
                    nameswap = self.beads.names[axlist[i]]
                    self.beads.names[axlist[i]] = self.beads.names[axlist[j]]
                    self.beads.names[axlist[j]] = nameswap
                    # change masses
                    massratio = s.beads.m[axlist[i]]/s.beads.m[axlist[j]]
                    self.beads.m[axlist[i]] /= massratio
                    self.beads.m[axlist[j]] *= massratio
                    # adjust the (classical) momenta
                    self.beads.p[3*axlist[i]:3*(axlist[i]+1)] /= np.sqrt(massratio)
                    self.beads.p[3*axlist[j]:3*(axlist[j]+1)] *= np.sqrt(massratio)
                    #print 'exchange atom No.  ', atomexchangelist[i], '  and  ', atomexchangelist[j]
                            
