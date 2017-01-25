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

from ipi.utils.depend import *
from ipi.utils.softexit import softexit
from ipi.utils.io import read_file
from ipi.utils.io.inputs.io_xml import xml_parse_file
from ipi.utils.units import unit_to_internal, Constants
from ipi.engine.thermostats import *
from ipi.engine.barostats import *
from ipi.engine.forces import Forces, ScaledForceComponent


__all__ = ['Ensemble', 'ensemble_swap']

### IMPORTANT - THIS MUST BE KEPT UP-TO-DATE WHEN THE ENSEMBLE CLASS IS CHANGED
def ensemble_swap(ens1, ens2):
    """ Swaps the definitions of the two ensembles, by
    exchanging all of the inner properties. """

    if ens1.temp != ens2.temp :
        ens1.temp, ens2.temp = ens2.temp, ens1.temp
    if ens1.pext != ens2.pext :
        ens1.pext, ens2.pext = ens2.pext, ens1.pext    
    if len(ens1.bweights) != len(ens2.bweights):
        raise ValueError("Cannot exchange ensembles that have different numbers of bias components")
    if len(ens1.hweights) != len(ens2.hweights):
        raise ValueError("Cannot exchange ensembles that are described by different forces")        
    if not np.array_equal(ens1.bweights, ens2.bweights):        
        ens1.bweights, ens2.bweights = ens2.bweights.copy(), ens1.bweights.copy()                
    if not np.array_equal(ens1.hweights, ens2.hweights):        
        ens1.hweights, ens2.hweights = ens2.hweights.copy(), ens1.hweights.copy()   
    
    print "force components for ", ens1.hweights
    for f in ens1.forces.mforces:        
        print f.pot
    print "bias forces", ens1.bias.pot, ens2.bias.pot
    for f in ens1.bias.mforces:        
        print f.pot
    

class Ensemble(dobject):
    """Base ensemble class.

    Defines the thermodynamic state of the system.

    Depend objects:
        temp: The system's temperature.
        pext: The systems's pressure
        stressext: The system's stress tensor
        bias: Explicit bias forces
    """

    def __init__(self, eens=0.0, econs=0.0, temp=None, pext=None, stressext=None, bcomponents=None, bweights=None, hweights=None):
        """Initialises Ensemble.

        Args:
            temp: The temperature.
            fixcom: An optional boolean which decides whether the centre of mass
                motion will be constrained or not. Defaults to False.
        """

        dset(self, "temp", depend_value(name='temp'))
        if temp is not None:
            self.temp = temp
        else:
            self.temp = 0.0

        dset(self, "stressext", depend_array(name='stressext', value=np.zeros((3,3), float)))
        if stressext is not None:
            self.stressext = np.reshape(np.asarray(stressext), (3,3))
        else:
            self.stressext = 0.0

        dset(self, "pext", depend_value(name='pext'))
        if pext is not None:
            self.pext = pext
        else:
            self.pext = 0.0

        dset(self, "eens", depend_value(name='eens'))
        if eens is not None:
            self.eens = eens
        else:
            self.eens = 0.0
        
        # these are the additional bias components
        if bcomponents is None:
            bcomponents = []
        self.bcomp = bcomponents
        self.bias = Forces()
        
        # and their weights
        if bweights is None:
            bweights = np.ones(len(self.bcomp))
        
        dset(self, "bweights", depend_array(name="bweights", value = np.asarray(bweights)) )
        
        # weights of the Hamiltonian scaling
        if hweights is None:
            hweights = np.ones(0)
        self.hweights = np.asarray(hweights)
        

    def bind(self, beads, nm, cell, bforce, fflist, elist=[], xlpot=[], xlkin=[]):
        self.beads = beads
        self.cell = cell
        self.forces = bforce
        self.nm = nm

        self.bias.bind(self.beads, self.cell, self.bcomp, fflist)
        
        dset(self, "econs", depend_value(name='econs', func=self.get_econs))        
        # dependencies of the conserved quantity
        dget(self, "econs").add_dependency(dget(self.nm, "kin"))
        dget(self, "econs").add_dependency(dget(self.forces, "pot"))
        dget(self, "econs").add_dependency(dget(self.bias, "pot"))
        dget(self, "econs").add_dependency(dget(self.beads, "vpath"))
        dget(self, "econs").add_dependency(dget(self, "eens"))

        # pipes the weights to the list of weight vectors
        i = 0
        for fc in self.bias.mforces:            
            if fc.weight != 1:
                warning("The weight given to forces used in an ensemble bias are given a weight determined by bias_weight")
            deppipe(self, "bweights", fc, "weight", i)
            i += 1
        
        # add Hamiltonian REM bias components
        if len(self.hweights) == 0:
            self.hweights = np.ones(len(self.forces.mforces))
            
        dset(self, "hweights", depend_array(name="hweights", value = np.asarray(self.hweights)) )
        
        for ic in xrange(len(self.forces.mforces)):
            sfc=ScaledForceComponent(self.forces.mforces[ic],1.0)
            self.bias.add_component(self.forces.mbeads[ic], self.forces.mrpc[ic], sfc)
            dget(sfc,"scaling")._func = lambda i=ic: self.hweights[i]-1
            dget(sfc,"scaling").add_dependency(dget(self, "hweights"))        
        
        self._elist = []

        for e in elist:
            self.add_econs(e)

        dset(self, "lpens", depend_value(name='lpens', func=self.get_lpens, 
                dependencies=[ dget(self,"temp") ] ))        
        dget(self, "lpens").add_dependency(dget(self.nm, "kin"))
        dget(self, "lpens").add_dependency(dget(self.forces, "pot"))
        dget(self, "lpens").add_dependency(dget(self.bias, "pot"))
        dget(self, "lpens").add_dependency(dget(self.beads, "vpath"))
        
        # extended Lagrangian terms for the ensemble
        self._xlpot = []
        for p in xlpot:
            self.add_xlpot(p)
            
        self._xlkin = []
        for k in xlkin:
            self.add_xlkin(k)
        

    def add_econs(self, e):
        self._elist.append(e)
        dget(self, "econs").add_dependency(e)

    def add_xlpot(self, p):
        self._xlpot.append(p)
        dget(self, "lpens").add_dependency(p)
        
    def add_xlkin(self, k):
        self._xlkin.append(k)
        dget(self, "lpens").add_dependency(k)
        
    def get_econs(self):
        """Calculates the conserved energy quantity for constant energy
        ensembles.
        """
        eham = self.beads.vpath*self.nm.omegan2 + self.nm.kin + self.forces.pot
        eham += self.bias.pot   # bias
        for e in self._elist:
            eham += e.get()

        return eham + self.eens
        
    def get_lpens(self):
        """Returns the ensemble probability (modulo the partition function) 
        for the ensemble. 
        """
        
        
        lpens = (self.forces.pot+self.bias.pot+self.nm.kin+self.beads.vpath*self.nm.omegan2);
        
        # inlcude terms associated with an extended Lagrangian integrator of some sort
        for p in self._xlpot:
            lpens += p.get()
        for k in self._xlkin:
            lpens += k.get()
            
        lpens *= -1.0/(Constants.kb*self.temp*self.beads.nbeads)
        return lpens
