"""Interface to the calculators of the Atomic Simulation Environment.

Atomic Simulation Environment:
https://wiki.fysik.dtu.dk/ase/
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import time

import numpy as np

from ..interfaces.sockets import Client
from ..utils import units


class ClientASE(Client):
    """Socket client that calls an ASE calculator to get interactions.

    Atomic Simulation Environment:
    https://wiki.fysik.dtu.dk/ase/
    """

    def __init__(self, atoms, address='localhost', port=31415, mode='unix', verbose=False):
        """Store provided data and initialize the base class.

        `atoms` is an ASE `Atoms` object, the rest gets passed to `Client`.
        """

        # store the provided data
        self.atoms = atoms

        # prepare unit conversions
        # (ASE uses Angstrom and eV)
        self.eV = units.unit_to_internal('energy', 'electronvolt', 1.0)
        self.Angstrom = units.unit_to_internal('length', 'angstrom', 1.0)

        # prepare arrays
        self._positions = np.zeros_like(atoms.get_positions())
        self._force = np.zeros_like(atoms.get_positions())
        self._potential = np.zeros(1)

        # call base class constructor
        super(ClientASE, self).__init__(address, port, mode, verbose)

    def _getforce(self):
        """Update stored potential energy and forces using ASE."""

        # for convenience
        atoms = self.atoms

        # update current coordinates and cell
        atoms.set_positions(self._positions / self.Angstrom)
        atoms.set_cell(self._cellh / self.Angstrom)

        # get data out, trigger calculation in the process
        self._force[:] = atoms.get_forces() * self.eV / self.Angstrom
        self._potential[:] = np.array([atoms.get_potential_energy() * self.eV])

        # DEBUG
        #print 'positions for ASE:'
        #print self._positions / self.Angstrom
        #print
        #print 'cell for ASE:'
        #print self._cellh / self.Angstrom
        #print
        #print 'potential energy [Ha]:'
        #print self._potential
        #print
        #print 'forces [atomic units]:'
        #print self._force
        #print
