"""Contains the classes that connect the driver to the python code.

ForceField objects are force providers, i.e. they are the abstraction
layer for a driver that gets positions and returns forces (and energy).
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import time
import threading

import numpy as np

from ipi.utils.softexit import softexit
from ipi.utils.messages import verbosity
from ipi.utils.messages import info
from ipi.interfaces.sockets import InterfaceSocket
from ipi.utils.depend import dobject
from ipi.utils.depend import depstrip


__all__ = ['ForceField', 'FFSocket', 'FFLennardJones', 'FFDebye']


class ForceRequest(dict):
    """An extension of the standard Python dict class which only has a == b
    if a is b == True, rather than if the elements of a and b are identical.

    Standard dicts are checked for equality if elements have the same value.
    Here I only care if requests are instances of the very same object.
    This is useful for the `in` operator, which uses equality to test membership.
    """

    def __eq__(self, y):
        """Overwrites the standard equals function."""
        return self is y


class ForceField(dobject):
    """Base forcefield class.

    Gives the standard methods and quantities needed in all the forcefield
    classes.

    Attributes:
        pars: A dictionary of the parameters needed to initialize the forcefield.
            Of the form {'name1': value1, 'name2': value2, ... }.
        name: The name of the forcefield.
        latency: A float giving the number of seconds the socket will wait
            before updating the client list.
        requests: A list of all the jobs to be given to the client codes.
        dopbc: A boolean giving whether or not to apply the periodic boundary
            conditions before sending the positions to the client code.
        _thread: The thread on which the socket polling loop is being run.
        _doloop: A list of booleans. Used to decide when to stop running the
            polling loop.
        _threadlock: Python handle used to lock the thread held in _thread.
    """

    def __init__(self, latency=1.0, name="", pars=None, dopbc=True):
        """Initialises ForceField.

        Args:
            latency: The number of seconds the socket will wait before updating
                the client list.
            name: The name of the forcefield.
            pars: A dictionary used to initialize the forcefield, if required.
                Of the form {'name1': value1, 'name2': value2, ... }.
            dopbc: Decides whether or not to apply the periodic boundary conditions
                before sending the positions to the client code.
        """

        if pars is None:
            self.pars = {}
        else:
            self.pars = pars

        self.name = name
        self.latency = latency
        self.requests = []
        self.dopbc = dopbc
        self._thread = None
        self._doloop = [False]
        self._threadlock = threading.Lock()

    def queue(self, atoms, cell, reqid=-1):
        """Adds a request.

        Note that the pars dictionary need to be sent as a string of a
        standard format so that the initialisation of the driver can be done.

        Args:
            atoms: An Atoms object giving the atom positions.
            cell: A Cell object giving the system box.
            pars: An optional dictionary giving the parameters to be sent to the
                driver for initialisation. Defaults to {}.
            reqid: An optional integer that identifies requests of the same type,
               e.g. the bead index

        Returns:
            A list giving the status of the request of the form {'pos': An array
            giving the atom positions folded back into the unit cell,
            'cell': Cell object giving the system box, 'pars': parameter string,
            'result': holds the result as a list once the computation is done,
            'status': a string labelling the status of the calculation,
            'id': the id of the request, usually the bead number, 'start':
            the starting time for the calculation, used to check for timeouts.}.
        """

        par_str = " "

        if not self.pars is None:
            for k, v in self.pars.items():
                par_str += k + " : " + str(v) + " , "
        else:
            par_str = " "

        pbcpos = depstrip(atoms.q).copy()
        if self.dopbc:
            cell.array_pbc(pbcpos)

        newreq = ForceRequest({
            "id": reqid,
            "pos": pbcpos,
            "cell": (depstrip(cell.h).copy(), depstrip(cell.ih).copy()),
            "pars": par_str,
            "result": None,
            "status": "Queued",
            "start": -1,
            "t_queued": time.time(),
            "t_dispatched": 0,
            "t_finished": 0
            })

        self._threadlock.acquire()
        try:
            self.requests.append(newreq)
        finally:
            self._threadlock.release()

        return newreq

    def poll(self):
        """Polls the forcefield object to check if it has finished."""

        for r in self.requests:
            if r["status"] == "Queued":
                r["t_dispatched"] = time.time()
                r["result"] = [0.0, np.zeros(len(r["pos"]), float), np.zeros((3,3), float), ""]
                r["status"] = "Done"
                r["t_finished"] = time.time()

    def _poll_loop(self):
        """Polling loop.

        Loops over the different requests, checking to see when they have
        finished.
        """

        info(" @ForceField: Starting the polling thread main loop.", verbosity.low)
        while self._doloop[0]:
            time.sleep(self.latency)
            self.poll()

    def release(self, request):
        """Shuts down the client code interface thread.

        Args:
            request: The id of the job to release.
        """

        """Frees up a request."""

        self._threadlock.acquire()
        try:
            if request in self.requests:
                try:
                    self.requests.remove(request)
                except ValueError:
                    print "failed removing request", id(request), ' ',
                    print [id(r) for r in self.requests], "@", threading.currentThread()
                    raise
        finally:
            self._threadlock.release()

    def stop(self):
        """Dummy stop method."""

        self._doloop[0] = False
        for r in self.requests:
            r["status"] = "Exit"

    def run(self):
        """Spawns a new thread.

        Splits the main program into two threads, one that runs the polling loop
        which updates the client list, and one which gets the data.

        Raises:
            NameError: Raised if the polling thread already exists.
        """

        if not self._thread is None:
            raise NameError("Polling thread already started")

        self._doloop[0] = True
        self._thread = threading.Thread(target=self._poll_loop, name="poll_" + self.name)
        self._thread.daemon = True
        self._thread.start()
        softexit.register_function(self.softexit)
        softexit.register_thread(self._thread, self._doloop)

    def softexit(self):
        """ Takes care of cleaning up upon softexit """

        self.stop()


class FFSocket(ForceField):
    """Interface between the PIMD code and a socket for a single replica.

    Deals with an individual replica of the system, obtaining the potential
    force and virial appropriate to this system. Deals with the distribution of
    jobs to the interface.

    Attributes:
        socket: The interface object which contains the socket through which
            communication between the forcefield and the driver is done.
    """

    def __init__(self, latency=1.0, name="", pars=None, dopbc=True, interface=None):
        """Initialises FFSocket.

        Args:
           latency: The number of seconds the socket will wait before updating
              the client list.
           name: The name of the forcefield.
           pars: A dictionary used to initialize the forcefield, if required.
              Of the form {'name1': value1, 'name2': value2, ... }.
           dopbc: Decides whether or not to apply the periodic boundary conditions
              before sending the positions to the client code.
           interface: The object used to create the socket used to interact
              with the client codes.
        """

        # a socket to the communication library is created or linked
        super(FFSocket, self).__init__(latency, name, pars, dopbc)
        if interface is None:
            self.socket = InterfaceSocket()
        else:
            self.socket = interface
        self.socket.requests = self.requests

    def poll(self):
        """Function to check the status of the client calculations."""

        self.socket.poll()

    def run(self):
        """Spawns a new thread."""

        self.socket.open()
        super(FFSocket, self).run()

    def stop(self):
        """Closes the socket and the thread."""

        super(FFSocket, self).stop()
        if self._thread is not None:
            # must wait until loop has ended before closing the socket
            self._thread.join()
        self.socket.close()


class FFLennardJones(ForceField):
    """Basic fully pythonic force provider.

    Computes LJ interactions without minimum image convention, cutoffs or
    neighbour lists. Parallel evaluation with threads.

    Attributes:
        parameters: A dictionary of the parameters used by the driver. Of the
            form {'name': value}.
        requests: During the force calculation step this holds a dictionary
            containing the relevant data for determining the progress of the step.
            Of the form {'atoms': atoms, 'cell': cell, 'pars': parameters,
                         'status': status, 'result': result, 'id': bead id,
                         'start': starting time}.
    """

    def __init__(self, latency=1.0e-3, name="", pars=None, dopbc=False):
        """Initialises FFLennardJones.

        Args:
           pars: Optional dictionary, giving the parameters needed by the driver.
        """

        # check input - PBCs are not implemented here
        if dopbc:
            raise ValueError("Periodic boundary conditions are not supported by FFLennardJones.")

        # a socket to the communication library is created or linked
        super(FFLennardJones, self).__init__(latency, name, pars, dopbc=False)
        self.epsfour = float(self.pars["eps"]) * 4
        self.sixepsfour = 6 * self.epsfour
        self.sigma2 = float(self.pars["sigma"]) * float(self.pars["sigma"])

    def poll(self):
        """Polls the forcefield checking if there are requests that should
        be answered, and if necessary evaluates the associated forces and energy."""

        # We have to be thread-safe, as in multi-system mode this might get
        # called by many threads at once.
        self._threadlock.acquire()
        try:
            for r in self.requests:
                if r["status"] == "Queued":
                    r["status"] = "Running"
                    r["t_dispatched"] = time.time()
                    self.evaluate(r)
        finally:
            self._threadlock.release()

    def evaluate(self, r):
        """Just a silly function evaluating a non-cutoffed, non-pbc and
        non-neighbour list LJ potential."""

        q = r["pos"].reshape((-1, 3))
        nat = len(q)

        v = 0.0
        f = np.zeros(q.shape)
        for i in range(1, nat):
            dij = q[i] - q[:i]
            rij2 = (dij**2).sum(axis=1)

            x6 = (self.sigma2 / rij2)**3
            x12 = x6**2

            v += (x12 - x6).sum()
            dij *= (self.sixepsfour * (2.0*x12 - x6) / rij2)[:,np.newaxis]
            f[i] += dij.sum(axis=0)
            f[:i] -= dij

        v *= self.epsfour

        r["result"] = [v, f.reshape(nat*3), np.zeros((3,3), float), ""]
        r["status"] = "Done"

class FFDebye(ForceField):
   """Debye crystal harmonic reference potential

   Computes a harmonic forcefield. 

   Attributes:
      parameters: A dictionary of the parameters used by the driver. Of the
         form {'name': value}.
      requests: During the force calculation step this holds a dictionary
         containing the relevant data for determining the progress of the step.
         Of the form {'atoms': atoms, 'cell': cell, 'pars': parameters,
                      'status': status, 'result': result, 'id': bead id,
                      'start': starting time}.
   """
   
   def __init__(self, latency = 1.0, name = "",  H=None, xref=None, vref=0.0, pars=None, dopbc = False, threaded=True):
      """Initialises FFDebye.

      Args:
         pars: Optional dictionary, giving the parameters needed by the driver.
      """

      # a socket to the communication library is created or linked
      # NEVER DO PBC -- forces here are computed without.
      super(FFDebye,self).__init__(latency, name, pars, dopbc=False)
            
      if H is None:
          raise ValueError("Must provide the Hessian for the Debye crystal.")
      if xref is None:
          raise ValueError("Must provide a reference configuration for the Debye crystal.")
  		  
      self.H = H
      self.xref = xref
      self.vref = vref

      eigsys=np.linalg.eigh(self.H) 
      info(" @ForceField: Hamiltonian eigenvalues: " + ' '.join(map(str, eigsys[0])), verbosity.medium)                 

   def poll(self):
      """ Polls the forcefield checking if there are requests that should
      be answered, and if necessary evaluates the associated forces and energy. """

      # we have to be thread-safe, as in multi-system mode this might get called by many threads at once
      self._threadlock.acquire()
      try:
         for r in self.requests:
            if r["status"] == "Queued":
               r["status"] = "Running"

               self.evaluate(r)
      finally:
         self._threadlock.release()

   def evaluate(self, r):
      """ A simple evaluator for a harmonic Debye crystal potential. """

      q = r["pos"]
      n3 = len(q)
      if self.H.shape != (n3,n3): 
          raise ValueError("Hessian size mismatch")
      if self.xref.shape != (n3,): 
          raise ValueError("Reference structure size mismatch")
      
      d = q-self.xref
      mf = np.dot(self.H, d)
            
      r["result"] = [ self.vref + 0.5*np.dot(d,mf), -mf, np.zeros((3,3),float), ""]
      r["status"] = "Done"
      r["t_finished"] = time.time()
