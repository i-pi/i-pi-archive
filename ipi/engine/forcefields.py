"""Contains the classes that connect the driver to the python code.

Copyright (C) 2013, Joshua More and Michele Ceriotti

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


Communicates with the driver code, obtaining the force, virial and potential.
Deals with creating the jobs that will be sent to the driver, and
returning the results to the python code.

Classes:
   ForceField: Base forcefield class with the generic methods and attributes.
   FFSocket: Deals with a single replica of the system
   ForceBeads: Deals with the parallelization of the force calculation over
      different beads.
   Forces: Deals with the parallelizatoin of the force calculation over
      different forcefields.
"""

__all__ = ['ForceField', 'FFSocket']

import numpy as np
import sys, os
import time, threading
from ipi.utils.softexit import softexit
from ipi.utils.messages import verbosity, warning, info
from ipi.utils.depend import *
from ipi.utils.nmtransform import nm_rescale
from ipi.interfaces.sockets import InterfaceSocket
from ipi.engine.beads import Beads

# standard dicts are checked for equality if elements have the same value.
# here I only care if requests are instances of the very same object
class ForceRequest(dict):
   def __eq__(self, y):
      return self is y
      
class ForceField(dobject):
   """Base forcefield class.

   Gives the standard methods and quantities needed in all the forcefield
   classes.
   """

   def __init__(self, latency = 1.0, name = "", pars = None, dopbc = True):
      """Initialises ForceField."""

      if pars is None:
         self.pars = {}
      else:
         self.pars = pars

      self.name = name
      self.latency = latency
      self.requests = []
      self.dopbc = dopbc
      self._thread = None
      self._doloop = [ False ]
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
         A list giving the status of the request of the form {'atoms': Atoms
         object giving the atom positions, 'cell': Cell object giving the
         system box, 'pars': parameter string, 'result': holds the result as a
         list once the computation is done, 'status': a string labelling the
         status, 'id': the id of the request, usually the bead number, 'start':
         the starting time for the calculation, used to check for timeouts.}.
      """

      par_str = " "

      if not self.pars is None:
         for k,v in self.pars.items():
            par_str += k + " : " + str(v) + " , "
      else:
         par_str = " "

      pbcpos = depstrip(atoms.q).copy()
      if self.dopbc:
         cell.array_pbc(pbcpos)

      newreq = ForceRequest({ "id": reqid, "pos": pbcpos, "cell": ( depstrip(cell.h).copy(), depstrip(cell.ih).copy() ),
          "pars": par_str,
          "result": None, "status": "Queued",
          "start": -1 })

      self._threadlock.acquire()
      try:
         self.requests.append(newreq)
      finally:
         self._threadlock.release()

      return newreq

   def poll(self):
      for r in self.requests:
         if r["status"] == "Queued":
            r["result"] = [ 0.0, np.zeros(len(r["pos"]),float), np.zeros((3,3),float), ""]
            r["status"] = "Done"

   def _poll_loop(self):
      info(" @ForceField: Starting the polling thread main loop.", verbosity.low)
      while self._doloop[0]:
         time.sleep(self.latency)
         self.poll()

   def release(self, request):

      self._threadlock.acquire()
      try:
         if request in self.requests:   
            try:
               self.requests.remove(request)
            except:
               print "failed removing request", id(request), [id(r) for r in self.requests], "@", threading.currentThread()
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
   """Interface between the PIMD code and the socket for a single replica.

   Deals with an individual replica of the system, obtaining the potential
   force and virial appropriate to this system. Deals with the distribution of
   jobs to the interface.

   Attributes:
      parameters: A dictionary of the parameters used by the driver. Of the
         form {'name': value}.
      socket: The interface object which contains the socket through which
         communication between the forcefield and the driver is done.
      requests: During the force calculation step this holds a dictionary
         containing the relevant data for determining the progress of the step.
         Of the form {'atoms': atoms, 'cell': cell, 'pars': parameters,
                      'status': status, 'result': result, 'id': bead id,
                      'start': starting time}.
   """

   def __init__(self, latency = 1.0, name = "",  pars=None, dopbc = True, interface=None):
      """Initialises FFSocket.

      Args:
         pars: Optional dictionary, giving the parameters needed by the driver.
         interface: Optional Interface object, which contains the socket.
      """

      # a socket to the communication library is created or linked
      super(FFSocket,self).__init__(latency, name, pars, dopbc)
      if interface is None:
         self.socket = InterfaceSocket()
      else:
         self.socket = interface
      self.socket.requests = self.requests

   def poll(self):
      self.socket.poll()

   def run(self):
      self.socket.open()
      super(FFSocket,self).run()

   def stop(self):
      super(FFSocket,self).stop()
      if not self._thread is None:   # must wait until loop has ended before closing the socket
         self._thread.join()
      self.socket.close()

