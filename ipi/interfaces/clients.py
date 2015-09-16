"""Clients that get interactions from specific programs.

Here we implement classes that provide interactions obtained from specific
programs and serve them back to i-PI.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import sys
import socket
import time

import numpy as np

from .sockets import DriverSocket, Message
from ..utils import units


class Client(DriverSocket):
    """Base class for the implementation of a client in Python.

    Handles sending and receiving data from the client code.

    Attributes:
        havedata: Boolean giving whether the client calculated the forces.
    """

    def __init__(self, address="localhost", port=31415, mode="unix", verbose=False, _socket=True):
        """Initialise Client.

        Args:
            - address: A string giving the name of the host network.
            - port: An integer giving the port the socket will be using.
            - mode: A string giving the type of socket used - 'inet' or 'unix'.
            - _socket: If a socket should be opened. Can be False for testing purposes.
        """

        if _socket:
            # open client socket
            if mode == "inet":
                _socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                _socket.connect((address, int(port)))
            elif mode == "unix":
                try:
                    _socket = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
                    _socket.connect("/tmp/ipi_" + address)
                except socket.error:
                    print 'Could not connect to UNIX socket: %s' % ("/tmp/ipi_" + address)
                    sys.exit(1)
            else:
                raise NameError("Interface mode " + mode + " is not implemented (should be unix/inet)")
            super(Client, self).__init__(socket=_socket)
        else:
            super(Client, self).__init__(socket=None)

        self.verbose = verbose

        # allocate data
        self.havedata = False
        self._vir = np.zeros((3,3), np.float64)
        self._cellh = np.zeros((3,3), np.float64)
        self._cellih = np.zeros((3,3), np.float64)
        self._nat = np.int32()
        self._callback = None

    def _getforce(self):
        """Dummy _getforce routine.

        This function must be implemented by subclassing or providing a callback function.
        This function is assumed to calculate the following:
            - self._force: The force of the current positions at self._positions.
            - self._potential: The potential of the current positions at self._positions.
        """
        if self._callback is not None:
            self._force, self._potential = self._callback(self._positions)
        else:
            raise NotImplementedError("_getforce must be implemented by providing a self.callback function or overwritten.")

    def run(self):
        """Serve forces until asked to finish or socket disconnects.

        Serve force and potential, that are calculated in the user provided
        routine _getforce.
        """

        if self.verbose:
            print 'Starting communication loop.'

        try:
            while True:
                msg = self.recv_msg()
                if msg == "":
                    print "Server shut down."
                    break
                elif msg == Message("status"):
                    if self.havedata:
                        self.send_msg("havedata")
                    else:
                        self.send_msg("ready")
                elif msg == Message("posdata"):
                    self._cellh = self.recvall(self._cellh)
                    self._cellih = self.recvall(self._cellih)
                    self._nat = self.recvall(self._nat)
                    self._positions = self.recvall(self._positions)
                    t0 = time.time()
                    self._getforce()
                    if self.verbose:
                        print 'Calculation time: {:7.3f} s'.format(time.time() - t0)
                    self.havedata = True
                elif msg == Message("getforce"):
                    self.sendall(Message("forceready"))
                    self.sendall(self._potential, 8)
                    self.sendall(self._nat, 4)
                    self.sendall(self._force, 8*self._force.size)
                    self.sendall(self._vir, 9*8)
                    self.sendall(np.int32(0), 4)
                    self.havedata = False
                else:
                    print >> sys.stderr, "Client could not understand command:", msg
                    break
        except socket.error as e:
            print 'Error communicating through socket: [{}] {}'.format(e.errno, e.strerror)
        except KeyboardInterrupt:
            print ' Keyboard interrupt.'

        if self.verbose:
            print 'Communication loop finished.'


class ClientASE(Client):
    """Socket client that calls an ASE calculator to get interactions.

    Atomic Simulation Environment:
    https://wiki.fysik.dtu.dk/ase/
    """

    def __init__(self, atoms, verbose=False, address='localhost', port=31415, mode='unix', _socket=True):
        """Store provided data and initialize the base class.

        Arguments:
            - `atoms`: an ASE `Atoms` object
            - `verbose`: enable priting of wall clock time per step
            - the rest gets passed to the `Client` base class
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
        super(ClientASE, self).__init__(address, port, mode, verbose, _socket)

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
