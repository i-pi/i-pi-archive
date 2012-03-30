"""Deals with the socket communication between the PIMD and driver code.

Deals with creating the socket, transmitting and receiving data, accepting and
removing different driver routines and the parallelization of the force
calculation.

Classes:
   Status: Simple class to keep track of the status, uses bitwise or to give
      combinations of different status options.
   Driver: Class to deal with communication between a client and the driver
      code.
   RestartInterface: Deals with creating the Interface object from a file, and
      writing the checkpoints.
   Interface: Host server class. Deals with distribution of all the jobs 
      between the different client servers.

Functions:
   Message: Sends a header string through the socket.

Exceptions:
   Disconnected: Raised if client has been disconnected.
   InvalidStatus: Raised if client has the wrong status. Shouldn't have to be
      used if the structure of the program is correct.
"""

__all__ = ['Message', 'Disconnected', 'InvalidStatus', 'Status', 'Driver',
           'RestartInterface', 'Interface']

import socket, select, threading, signal, string, os, time
import numpy as np
from utils.restart import Restart, RestartValue

HDRLEN = 12
UPDATEFREQ = 100
TIMEOUT = 1.0 
SERVERTIMEOUT = 2.0*TIMEOUT

def Message(mystr):
   """Returns a header of standard length HDRLEN."""

   return string.ljust(string.upper(mystr), HDRLEN)


class Disconnected(Exception): 
   """Disconnected: Raised if client has been disconnected."""

   pass


class InvalidStatus(Exception):
   """InvalidStatus: Raised if client has the wrong status. 

   Shouldn't have to be used if the structure of the program is correct.
   """

   pass

class Status: 
   """Simple class used to keep track of the status of the client.

   Uses bitwise or to give combinations of different status options.
   i.e. Status.Up | Status.Ready would be understood to mean that the client
   was connected and ready to receive the position and cell data.

   Attributes: 
      Disconnected: Flag for if the client has disconnected.
      Up: Flag for if the client is running.
      Ready: Flag for if the client has ready to receive position and cell data.
      NeedsInit: Flag for if the client is ready to receive forcefield
         parameters.
      HasData: Flag for if the client is ready to send force data.
      Busy: Flag for if the client is busy.
      Timeout: Flag for if the connection has timed out.
   """

   Disconnected = 0
   Up = 1
   Ready = 2
   NeedsInit = 4
   HasData = 8
   Busy = 16
   Timeout = 32


class Driver(socket.socket):
   """Deals with communication between the client and driver code.

   Deals with sending and receiving the data from the driver code. Keeps track
   of the status of the driver. Initialises the driver forcefield, sends the
   position and cell data, and receives the force data.

   Attributes: 
      _buf: A string buffer to hold the reply from the driver.
      busyonstatus: Boolean giving whether the driver is busy.
      status: Keeps track of the status of the driver.
      lastreq: The ID of the last request processed by the client. 
   """

   def __init__(self, socket):
      """Initialises Driver.

      Args:
         socket: A socket through which the communication should be done.
      """

      super(Driver,self).__init__(_sock=socket)
      self._buf = np.zeros(0,np.byte)
      self.busyonstatus = False
      self.status = Status.Up
      self.lastreq = None
      
   def poll(self):
      """Waits for driver status."""

      self.status = self._getstatus()
      
   def _getstatus(self):   
      """Gets driver status.

      Returns:
         An integer labelling the status via bitwise or of the relevant members
         of Status.
      """

      if not self.busyonstatus:
         try:
            self.sendall(Message("status"))      
         except:
            return Status.Disconnected

      readable, writable, errored = select.select([self], [], [], 0.0)
      if not self in readable:
         self.busyonstatus = True
         return Status.Up | Status.Busy
      
      self.busyonstatus = False
      try:
         reply = self.recv(HDRLEN)
      except socket.timeout: 
         print " @SOCKET:   Timeout in status recv!"
         return Status.Up | Status.Busy | Status.Timeout
      except:
         return Status.Disconnected
         
      if not len(reply) == HDRLEN:
         return Status.Disconnected
      elif reply == Message("ready"):
         return Status.Up | Status.Ready
      elif reply == Message("needinit"):
         return Status.Up | Status.NeedsInit
      elif reply == Message("havedata"):
         return Status.Up | Status.HasData
      else:
         print " @SOCKET:    Unrecognized reply: ", reply         
         return Status.Up
   
   def recvall(self, dest):      
      """Gets the potential energy, force and virial from the driver.

      Args:
         dest: Object to be read into.

      Raises:
         Disconnected: Raised if client is disconnected.

      Returns:
         The data read from the socket to be read into dest.
      """

      blen = dest.itemsize*dest.size
      if (blen>len(self._buf)): 
         self._buf.resize(blen)
      bpos = 0
      while bpos < blen:
         timeout = False

#   pre-2.5 version. 
         try:
            bpart = 1
            bpart = self.recv( blen-bpos )
            self._buf[bpos:bpos+len(bpart)]=np.fromstring(bpart, np.byte)
         except socket.timeout:
            print " @SOCKET:   Timeout in status recvall, trying again!"
            timeout = True
            pass
         if (not timeout and bpart == 0):
            raise Disconnected()
         bpos+=len(bpart)

#   post-2.5 version: slightly more compact for modern python versions
#         try:
#            bpart = 1
#            bpart = self.recv_into(self._buf[bpos:], blen-bpos)
#         except socket.timeout:
#            print " @SOCKET:   Timeout in status recvall, trying again!"
#            timeout = True
#            pass
#         if (not timeout and bpart == 0):
#            raise Disconnected()
#         bpos += bpart
#TODO this Disconnected() exception currently just causes the program to hang.
#This should do something more graceful

      if np.isscalar(dest):
         return np.fromstring(self._buf[0:blen], dest.dtype)[0]
      else:
         return np.fromstring(self._buf[0:blen], dest.dtype).reshape(dest.shape)
       
   def initialize(self, pars):
      """Sends the initialisation string to the driver.

      Args:
         pars: The parameter string to be sent to the driver.

      Raises:
         InvalidStatus: Raised if the status is not NeedsInit.
      """

      if self.status & Status.NeedsInit:
         try:
            self.sendall(Message("init"))
            self.sendall(np.int32(len(pars)))
            self.sendall(pars)                        
         except: 
            self.poll()
            return
      else:
         raise InvalidStatus("Status in init was " + self.status)
               
   def sendpos(self, atoms, cell):
      """Sends the position and cell data to the driver.

      Args:
         atoms: An Atoms object giving the atom positions.
         cell: A cell object giving the system box.

      Raises:
         InvalidStatus: Raised if the status is not Ready.
      """

      if (self.status & Status.Ready):
         try:
            self.sendall(Message("posdata"))
            self.sendall(cell.h, 9*8)
            self.sendall(cell.ih, 9*8)      
            self.sendall(np.int32(len(atoms)))
            self.sendall(atoms.q, len(atoms)*3*8)
         except: 
            self.poll()
            return
      else:
         raise InvalidStatus("Status in sendpos was " + self.status)

   def getforce(self):
      """Gets the potential energy, force and virial from the driver.

      Raises:
         InvalidStatus: Raised if the status is not HasData.
         Disconnected: Raised if the driver has disconnected.

      Returns:
         A list of the form [potential, force, virial].
      """

      if (self.status & Status.HasData):
         self.sendall(Message("getforce"));
         reply = ""
         while True:
            try:
               reply = self.recv(HDRLEN)
            except socket.timeout: 
               print " @SOCKET:   Timeout in getforce, trying again!"
               continue
            if reply == Message("forceready"):
               break          
            else:
               print " @SOCKET:   Unexpected getforce reply: ", reply
            if reply == "":
               raise Disconnected()
      else:
         raise InvalidStatus("Status in getforce was " + self.status)      

      mu = np.float64()
      mu = self.recvall(mu)

      mlen = np.int32()
      mlen = self.recvall(mlen)
      mf = np.zeros(3*mlen,np.float64)
      mf = self.recvall(mf)

      mvir = np.zeros((3,3),np.float64)
      mvir = self.recvall(mvir)
      return [mu, mf, mvir]

                  
class RestartInterface(Restart):         
   """Interface restart class.

   Handles generating the apporopriate interface class from the xml
   input file, and generating the xml checkpoin tags and data from an
   instance of the object.

   Attributes:
      address: A string giving the host name.
      port: An integer giving the port used by the socket.
      slots: An integer giving the maximum allowed backlog of queued clients.
      latency: A float giving the number of seconds that the interface waits
         before updating the client list.
      mode: A string giving the type of socket used.
   """

   fields = { "address" : (RestartValue, (str, "localhost")), "port" : (RestartValue, (int,31415)),
            "slots" : (RestartValue, (int, 4) ), "latency" : (RestartValue, (float, 1e-3))  }
   attribs = { "mode": (RestartValue, (str, "unix") ) }

   def store(self, iface):
      """Takes an Interface instance and stores a minimal representation of it.

      Args:
         iface: An interface object.
      """

      self.latency.store(iface.latency)
      self.mode.store(iface.mode)
      self.address.store(iface.address)
      self.port.store(iface.port)
      self.slots.store(iface.slots)
      
   def fetch(self):
      """Creates an Interface object.

      Returns:
         An interface object with the appropriate socket given the attributes
         of the RestartInterface object.
      """

      return Interface(address=self.address.fetch(), port=self.port.fetch(), slots=self.slots.fetch(), mode=self.mode.fetch(), latency=self.latency.fetch())

            
class Interface(object):
   """Host server class.

   Deals with distribution of all the jobs between the different client servers
   and both initially and as clients either finish or are disconnected.
   Deals with cleaning up after all calculations are done. Also deals with the
   threading mechanism, and cleaning up if the interface is killed.

   Attributes:
      address: A string giving the name of the host network.
      port: An integer giving the port the socket will be using.
      slots: An integer giving the maximum allowed backlog of queued clients.
      mode: A string giving the type of socket used.
      latency: A float giving the number of seconds the interface will wait 
         before updating the client list.
      server: The socket used for data transmition.
      clients: A list of the driver clients connected to the server.
      requests: A list of all the jobs required in the current PIMD step.
      jobs: A list of all the jobs currently running.
      _poll_thread: The thread the poll loop is running on.
      _prev_kill: Holds the signals to be sent to clean up the main thread 
         when a kill signal is sent.
      _poll_true: A boolean giving whether the thread is alive.
   """

   def __init__(self, address="localhost", port=31415, slots=4, mode="unix", latency=1e-3):
      """Initialises interface.

      Args:
         address: An optional string giving the name of the host server.
            Defaults to 'localhost'.
         port: An optional integer giving the port number. Defaults to 31415.
         slots: An optional integer giving the maximum allowed backlog of 
            queueing clients. Defaults to 4.
         mode: An optional string giving the type of socket. Defaults to 'unix'.
         latency: An optional float giving the time in seconds the socket will 
            wait before updating the client list. Defaults to 1e-3.

      Raises:
         NameError: Raised if mode is not 'unix' or 'inet'.
      """

      self.address = address
      self.port = port
      self.slots = slots
      self.mode = mode
      self.latency = latency
      
      if self.mode == "unix":
         self.server = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
         self.server.bind("/tmp/wrappi_" + address)
      elif self.mode == "inet":
         self.server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
         self.server.bind((address,port))
      else:
         raise NameError("Interface mode " + self.mode + " is not implemented (should be unix/inet)")

      self.server.listen(slots)
      self.server.settimeout(SERVERTIMEOUT)
      self.clients = []
      self.requests = []
      self.jobs = []      

      self._poll_thread = None
      self._prev_kill = {}
      self._poll_true = False
      
   def queue(self, atoms, cell, pars=None, reqid=0):
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
         status}.
      """

      par_str = " "
      
      if not pars is None: 
         for k,v in pars.items():
            par_str += k + " : " + str(v) + " , "
      else:
         par_str = " "

      newreq = {"atoms": atoms, "cell": cell, "pars": par_str, "result": None, "status": "Queued", "id": reqid }
      self.requests.append(newreq)
      return newreq
      
   def release(self, request):
      """Empties the list of requests once finished.

      Args:
         request: A list of requests that are done.
      """

      if request in self.requests:         
         self.requests.remove(request)
         
   def pool_update(self):
      """Deals with keeping the pool of client drivers up-to-date during a 
      force calculation step.

      Deals with maintaining the client list. Clients that have
      disconnected are removed and their jobs removed from the list of
      running jobs and new clients are connected to the server.
      """
     
      for c in self.clients[:]:         
         if not (c.status & Status.Up):
            try:
               print " @SOCKET:   Client died or got unresponsive. Removing from the list."
               c.shutdown(socket.SHUT_RDWR)
               c.close()
            except:   pass
            self.clients.remove(c)
            for [k,j] in self.jobs[:]:
               if j is c:
                  self.jobs.remove([k,j])
               k["status"] = "Queued"

      keepsearch = True
      while keepsearch:
         readable, writable, errored = select.select([self.server], [], [], 0.0)
         if self.server in readable:
            client, address = self.server.accept()
            client.settimeout(TIMEOUT)
            driver = Driver(client)
            print " @SOCKET:   Client asked for connection from ", address, ". Now hand-shaking."
            driver.poll()
            if (driver.status | Status.Up):
               self.clients.append(driver)
               print " @SOCKET:   Handshaking was successful. Added to the client list."
            else:
               print " @SOCKET:   Handshaking failed. Dropping connection."
               client.shutdown(socket.SHUT_RDWR)
               client.close()
         else:
            keepsearch = False
         
   def pool_distribute(self):
      """Deals with keeping the list of jobs up-to-date during a force 
      calculation step.

      Deals with maintaining the jobs list. Gets data from drivers that have
      finished their calculation and removes that job from the list of running
      jobs, adds jobs to free clients and initialises the forcefields of new
      clients.
      """

      for c in self.clients:
         if not c.status & ( Status.Ready | Status.NeedsInit ):
            c.poll()

      for [r,c] in self.jobs[:]:
         if c.status & Status.HasData:
            try:
               r["result"] = c.getforce()
            except Disconnected:
               c.status = 0
               continue
            c.poll()
            while c.status & Status.Busy:
               c.poll()
            if not (c.status & Status.Up): 
               print " @SOCKET:   Client died a horrible death while getting forces. Will try to cleanup."
               continue
            r["status"] = "Done"
            c.lastreq = r["id"] # saves the ID of the request that the client has just processed
            self.jobs.remove([r,c])
                     
      for r in self.requests:
         if r["status"] == "Queued":
            freec = self.clients[:]
            for [r2, c] in self.jobs:
               freec.remove(c)            
            
            for match_ids in ( True, False):
               matched = False
               for fc in freec:
                  if not (fc.status & Status.Up):
                     self.clients.remove(fc)
                     continue                  
                  if fc.status & Status.HasData:
                     continue                              
                  if not (fc.status & (Status.Ready | Status.NeedsInit | Status.Busy) ): 
                     print " @SOCKET:   (1) Client is in an unexpected status ",fc.status, ". Will try to keep calm and carry on."
                     continue
                  # First, tries to match request ids and lastreq clients.
                  # If it can't match on the first round, gives up and assigns requests
                  # on a first-come first-serve basis
                  if match_ids and not fc.lastreq is r["id"]:
                     continue

                  print " @SOCKET: Assigning request id ", r["id"], " to client with last-id ", fc.lastreq
                  while fc.status & Status.Busy:
                     fc.poll()            
                  if fc.status & Status.NeedsInit:
                     fc.initialize(r["pars"])
                     fc.poll()               
                     while fc.status & Status.Busy: # waits for initialization to finish. hopefully this is fast
                        fc.poll()
                  if fc.status & Status.Ready:
                     fc.sendpos(r["atoms"], r["cell"])
                     r["status"] = "Running"
                     fc.poll()
                     self.jobs.append([r,fc])
                     matched = True
                     break
                  else:
                     print " @SOCKET:   (2) Client is in an unexpected status ",fc.status,". Will try to keep calm and carry on."
               if matched: 
                  break # doesn't do a second round if it managed to 
                        # assign the job

   def _kill_handler(self, signal, frame):
      """Deals with handling a kill call gracefully.

      Prevents any of the threads becoming zombies, by intercepting a
      kill signal using the standard python function signal.signal() and
      then closing the socket and the spawned threads before closing the main
      thread. Called when signals SIG_INT and SIG_TERM are received.

      Args:
         signal: An integer giving the signal number of the signal received 
            from the socket.
         frame: Current stack frame.
      """

      print " @SOCKET:   Kill signal. Trying to make a clean exit."
      self.end_thread()
      try:
         self.__del__()
      except:
         pass      
      if signal in self._prev_kill:
         self._prev_kill[signal](signal, frame)
      
   def _poll_loop(self):   
      """The main thread loop.

      Runs until either the program finishes or a kill call is sent. Updates
      the pool of clients every UPDATEFREQ loops and loops every latency
      seconds until _poll_true becomes false.
      """
      
      print " @SOCKET:   Starting the polling thread main loop."
      poll_iter = 0
      while self._poll_true:
         time.sleep(self.latency)
         # makes sure to remove the last dead client as soon as possible. 
         if poll_iter > UPDATEFREQ or (len(self.clients) > 0 and not(self.clients[0].status & Status.Up)):
            self.pool_update()
            poll_iter = 0
         poll_iter += 1
         self.pool_distribute()
         
         if os.path.exists("EXIT"): # soft-exit
            print " @SOCKET:   Soft exit request. Flushing job queue."
            freec = self.clients[:]     
            # releases all pending requests       
            for r in self.requests:
               r["status"] = "Exit"
            for fc in freec:
               self.clients.remove(fc)
      self._poll_thread = None   
   
   def start_thread(self):
      """Spawns a new thread.

      Splits the main program into two threads, one that runs the polling loop
      which updates the client list, and one which gets the data. Also sets up
      the machinery to deal with a kill call, in the case of a Ctrl-C or
      similar signal the signal is intercepted by the _kill_handler function, 
      which cleans up the spawned thread before closing the main thread.

      Raises:
         NameError: Raised if the polling thread already exists.
      """

      if not self._poll_thread is None:
         raise NameError("Polling thread already started")      
      self._poll_thread = threading.Thread(target=self._poll_loop, name="poll_" + self.address)
      self._poll_thread.daemon = True
      self._prev_kill[signal.SIGINT] = signal.signal(signal.SIGINT, self._kill_handler)
      self._prev_kill[signal.SIGTERM] = signal.signal(signal.SIGTERM, self._kill_handler)    
      self._poll_true = True      
      self._poll_thread.start()
   
   def end_thread(self):
      """Closes the spawned thread.

      Deals with cleaning up the spawned thread cleanly. First sets 
      _poll_true to false to indicate that the poll_loop should be exited, then
      closes the spawned thread and removes it.
      """

      self._poll_true = False
      if not self._poll_thread is None:
         self._poll_thread.join()
      self._poll_thread = None
   
   def __del__(self):
      """Removes the socket.

      Closes the interface and removes the socket. Also unlinks the socket if
      possible, so that the port can be reused immediately.
      """

      print " @SOCKET:   Shutting down the server interface."
      self.server.shutdown(socket.SHUT_RDWR)
      self.server.close()
      if self.mode == "unix":
         os.unlink("/tmp/wrappi_" + self.address)
