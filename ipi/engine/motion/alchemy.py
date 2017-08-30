"""Contains the classes that deal with the MC exchanges of isotopes.

Holds the algorithms required for alchemical exchanges. Also calculates the
appropriate conserved energy quantity.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import time

import numpy as np

from ipi.engine.motion import Motion
from ipi.utils.depend import *
from ipi.utils.units import Constants


#__all__ = ['AlchemyMC']

class AlchemyMC(Motion):
    """Monte Carlo alchemical exchanges class.

    Attributes:
        names of the isotopes for exchanges

    Depend objects:
        The spring potential energy, names of the atoms.
    """

    def __init__(self, fixcom=False, fixatoms=None, mode=None, names=[], nmc=1, ealc=None):
        """Initialises a "alchemical exchange" motion object.

        Args: 
            names : A list of isotopes
            nmc : frequency of doing exchanges
            
        """

        super(AlchemyMC, self).__init__(fixcom=fixcom, fixatoms=fixatoms)

        self.names = names
        self.nmc = nmc    
  
        dself = dd(self)
        dself.ealc = depend_value(name='ealc')
        if not ealc is None:
            self.ealc = ealc
        else: self.ealc = 0.0

    def bind(self, ens, beads, cell, bforce, nm, prng):
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

        super(AlchemyMC, self).bind(ens, beads, cell, bforce, nm, prng)
        self.ensemble.add_econs(dget(self, "ealc"))

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

        if (1.0/self.nmc < self.prng.u) : return  # tries a round of exhanges with probability 1/nmc

        """Does one round of alchemical exchanges."""
        # record the spring energy (divided by mass) for each atom in the exchange chain
        q = depstrip(self.beads.q)
        nb = self.beads.nbeads
        axlist = self.AXlist(self.names)
        lenlist = len(axlist)
        atomspring = np.zeros(lenlist)
        i = 0
        for atomnum in axlist:
            na3 = atomnum * 3
            spr = 0.0
            for b in range(1,nb):
                for j in range(na3,na3+3):
                    spr += (q[b,j]-q[b-1,j])**2
            for j in range(na3,na3+3):
                spr += (q[nb-1,j]-q[0,j])**2
            # no mass here
            spr *= 0.5*self.nm.omegan2
            atomspring[i] = spr
            #print axlist[i], self.beads.names[axlist[i]], spr*self.beads.m[axlist[i]]
            i += 1

       
            
        # do the exchange
        betaP = 1.0/(Constants.kb*self.ensemble.temp*nb)
        nexch = 0
        for i in range(lenlist):
            for j in range(i):
                # no exchange for same type of atoms
                if (self.beads.names[axlist[i]] == self.beads.names[axlist[j]]) : continue
                # energy increase due to the swap
                difspring = (atomspring[i]-atomspring[j])*(self.beads.m[axlist[j]]-self.beads.m[axlist[i]])
                pexchange = np.exp(-betaP*difspring)
                #print 'exchange probablity: %10.5e  n. exchanges this far: %5d' % ( pexchange, nexch )
                pexchange=1.0
                # attemps the exchange
                if (pexchange > self.prng.u):
                    nexch += 1
                    oldspringenergy = self.beads.vpath*self.nm.omegan2
                    #print oldspringenergy
                    print 'exchange atom No.  ', axlist[i], '  and  ', axlist[j]
                    #print 'bofore exchange', 'econs:',self.ensemble.econs,'kin:', self.nm.kin
                    # swap names
                    nameswap = self.beads.names[axlist[i]]
                    self.beads.names[axlist[i]] = self.beads.names[axlist[j]]
                    self.beads.names[axlist[j]] = nameswap
                    # change masses
                    massratio = self.beads.m[axlist[i]]/self.beads.m[axlist[j]]
                    self.beads.m[axlist[i]] /= massratio
                    self.beads.m[axlist[j]] *= massratio
                    # adjust the (classical) momenta
                    self.beads.p[:,3*axlist[i]:3*(axlist[i]+1)] /= np.sqrt(massratio)
                    self.beads.p[:,3*axlist[j]:3*(axlist[j]+1)] *= np.sqrt(massratio)
                    #print 'exchange atom No.  ', axlist[i], '  and  ', axlist[j]
                    # adjusts the conserved quantities
                    # change in spring energy
                    newspringenergy = self.beads.vpath*self.nm.omegan2
                    #print oldspringenergy, newspringenergy, newspringenergy-oldspringenergy-difspring
                    #print 'after exchange', 'econs:',self.ensemble.econs,'kin:', self.nm.kin
                    self.ealc += -difspring
                    #print 'after adjust', 'econs:', self.ensemble.econs, 'ealc:', self.ealc

                            
