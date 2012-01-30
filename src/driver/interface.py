import socket, select, threading, signal, string, os
from utils.restart import Restart, RestartValue
import numpy as np
import time, pdb

HDRLEN=12
UPDATEFREQ=100
TIMEOUT=1.0 # try to trigger more timeout events
SERVERTIMEOUT=2.0*TIMEOUT

def Message(mystr): return string.ljust(string.upper(mystr),HDRLEN)

class Disconnected(Exception): pass
class InvalidStatus(Exception): pass
class Status: 
   Disconnected=0;   Up = 1;  Ready=2;   NeedsInit=4;  HasData=8;  Busy=16;   Timeout=32
   
class Driver(socket.socket):
   def __init__(self, socket):
      super(Driver,self).__init__(_sock=socket)
      self._buf=np.zeros(0,np.byte)
      self.busyonstatus=False
      self.status=Status.Up
      
   def poll(self):
      self.status=self._getstatus()
      
   def _getstatus(self):   
      if not self.busyonstatus:
         try:
            self.sendall(Message("status"))      
         except:
            return Status.Disconnected

      readable, writable, errored = select.select([self], [], [], 0.0)
      if not self in readable: 
         self.busyonstatus=True
         return Status.Up | Status.Busy
      
      self.busyonstatus=False
      try:         
         reply=self.recv(HDRLEN)
      except socket.timeout: 
         print "timeout in status recv"
         return Status.Up | Status.Busy | Status.Timeout
      except:
         return Status.Disconnected
         
      #print "status:", reply

      if not len(reply) == HDRLEN:
         return Status.Disconnected
      elif reply==Message("ready"):
         return Status.Up | Status.Ready
      elif reply==Message("needinit"):
         return Status.Up | Status.NeedsInit
      elif reply==Message("havedata"):
         return Status.Up | Status.HasData
      else:
         print "Unrecognized reply", reply         
         #return Status.Up
         return Status.Disconnected
   
   def recvall(self, dest):      
      blen=dest.itemsize*dest.size
      if (blen>len(self._buf)) : self._buf.resize(blen)  # keeps a permanent buffer, which is expanded if necessary
      bpos=0
      bpart=0
      while bpos<blen:
         timeout=False
         try:
            bpart = 1
            bpart=self.recv_into(self._buf[bpos:], blen-bpos )
         except socket.timeout:
            print "timeout in recvall, trying again"; timeout=True; pass
         if (not timeout and bpart == 0 ): raise Disconnected()
         bpos+=bpart
#TODO this Disconnected() exception currently just causes the program to hang.
#This should do something more graceful

#   pre-2.5 version. needed to run on the good ole' neptune
#         try:
#            bpart = 0
#            bpart=self.recv( blen-bpos )
#            self._buf[bpos:bpos+len(bpart)]=np.fromstring(bpart, np.byte)
#         except socket.timeout:
#            print "timeout in recvall, trying again"; pass
#         if (len(bpart) == 0 ): raise Disconnected()
#         bpos+=len(bpart)

      if np.isscalar(dest):       
         return np.fromstring(self._buf[0:blen], dest.dtype)[0]
      else:
         return np.fromstring(self._buf[0:blen], dest.dtype).reshape(dest.shape)
      
       
   def initialize(self, pars):
      if self.status & Status.NeedsInit:
         try:
            self.sendall(Message("init"));
            self.sendall(np.int32(len(pars)))
            self.sendall(pars)                        
         except: 
            self.poll(); return
      else: raise InvalidStatus("Status in init was "+self.status)
               
   def sendpos(self, atoms, cell):
      if (self.status & Status.Ready):
         try:
            self.sendall(Message("posdata"));
            self.sendall(cell.h, 9*8)
            self.sendall(cell.ih,9*8)      
            self.sendall(np.int32(len(atoms)))
            self.sendall(atoms.q,len(atoms)*3*8)
         except: 
            self.poll(); return
      else: raise InvalidStatus("Status in sendpos was "+self.status)

   def getforce(self):
      if (self.status & Status.HasData):
         self.sendall(Message("getforce"));
         reply=""
         while True:
            try:
               reply=self.recv(HDRLEN)
            except socket.timeout: 
               print "timeout in recvforce, try again!!"
               continue
            if reply==Message("forceready"): break          
            else: print "oh-oh. got ", reply, "in getforce"
            #pdb.set_trace()
            if reply=="": raise Disconnected()
      else: raise InvalidStatus("Status in getforce was "+self.status)      
      mu=np.float64(); mu=self.recvall(mu)         
      mlen=np.int32(); mlen=self.recvall(mlen)

      mf=np.zeros(3*mlen,np.float64)
      mf=self.recvall(mf)

      mvir=np.zeros((3,3),np.float64)
      mvir=self.recvall(mvir)
      return  [ mu, mf, mvir ]
                  
class RestartInterface(Restart):         
   fields={ "address" : (RestartValue, (str, "localhost")), "port" : (RestartValue, (int,31415)),
            "slots" : (RestartValue, (int, 4) ), "latency" : (RestartValue, (float, 1e-3))  }
   attribs={ "mode": (RestartValue, (str, "unix") ) }

   def store(self, iface):
      self.latency.store(iface.latency)
      self.mode.store(iface.mode)
      self.address.store(iface.address)
      self.port.store(iface.port)
      self.slots.store(iface.slots)
      
   def fetch(self):
      return Interface(address=self.address.fetch(), port=self.port.fetch(), slots=self.slots.fetch(), mode=self.mode.fetch(), latency=self.latency.fetch())

            
class Interface(object):

   def __init__(self, address="localhost", port=31415, slots=4, mode="unix", latency=1e-3):
      self.address = address; self.port = port; self.slots = slots; self.mode=mode; self.latency=latency
      
      if self.mode == "unix":
         self.server = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
         self.server.bind("/tmp/wrappi_"+address)
      elif self.mode=="inet":
         self.server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
         self.server.bind((address,port))
      else:
         raise NameError("Interface mode "+self.mode+" is not implemented (shall be unix/inet)")

      self.server.listen(slots)
      self.server.settimeout(SERVERTIMEOUT)
      self.clients=[]
      self.requests=[]
      self.jobs=[]      

      self._poll_thread=None
      self._prev_kill={}
      self._poll_true=False
      self.time_update=0.0
      self.time_distribute=0.0      
      
   def queue(self, atoms,cell, pars={}):
      if not pars == {}: 
         par_str = str(pars["eps"]) + " " + str(pars["sigma"]) + " " + str(pars["cutoff"]) + " " + str(pars["nearest_neighbour"])
      else : par_str=" "
      newreq={"atoms":atoms,"cell":cell, "pars":par_str, "result":None, "status":"Queued"}
      self.requests.append(newreq);
      return newreq
      
   def release(self, request):
      if request in self.requests:         
         self.requests.remove(request)
         
   def pool_update(self):
      start=time.time()
      # first, remove clients which had disconnected
      for c in self.clients[:]:         
         if not (c.status & Status.Up):
            try:
               print "removing dead client"
               c.shutdown(socket.SHUT_RDWR); c.close()
            except: "Exception removing client"; pass
            self.clients.remove(c)
            # also make sure that there are no jobs involving client c
            for [k,j] in self.jobs[:]:
               if j is c: self.jobs.remove([k,j])
               k["status"]="Queued"

      keepsearch=True
      while keepsearch:
         readable, writable, errored = select.select([self.server], [], [], 0.0)         
         if self.server in readable:
            client, address = self.server.accept()
            client.settimeout(TIMEOUT)
            driver=Driver(client)
            print "client ", address, "accepted, now handshaking"
            driver.poll()
            if (driver.status | Status.Up):
               self.clients.append(driver)
               print "handshake completed"
            else:
               print "handshake failed"
               client.shutdown(socket.SHUT_RDWR); client.close()
         else: keepsearch=False
         
      self.time_update+=time.time()-start
         
   def pool_distribute(self):
      start=time.time()
      # first, pools which were busy
      for c in self.clients:
         if not c.status & ( Status.Ready | Status.NeedsInit ): c.poll()

      # then checks for finished jobs      
      for [r,c] in self.jobs[:]:
         if c.status & Status.HasData:
            try:
               r["result"]=c.getforce()
            except Disconnected:
               c.status=0
               continue
            c.poll(); 
            while c.status & Status.Busy: c.poll()
            if not (c.status & Status.Up): 
               print "Client died badly while getting forces!"
               continue
            r["status"]="Done"
            self.jobs.remove([r,c])
                     
      # then starts new queued requests, if possible
      for r in self.requests:
         if r["status"] == "Queued":
            # gets a free client (if any)
            freec=self.clients[:]
            for [r2, c] in self.jobs:
               freec.remove(c)            
            
            for fc in freec:
               if not (fc.status & Status.Up): # client died waiting for a request, must remove and clean up
                  self.clients.remove(fc)
                  continue                  
               if fc.status & Status.HasData : continue                              
               if not (fc.status & (Status.Ready | Status.NeedsInit | Status.Busy) ): 
                  print "something weird is going on with a client, status is:",fc.status," - trying again later"
                  continue
               
#               fc.poll()
               while fc.status & Status.Busy: fc.poll()
                              
               if fc.status & Status.NeedsInit:
                  fc.initialize(r["pars"])
                  fc.poll()               
                  while fc.status & Status.Busy: fc.poll()
               if fc.status & Status.Ready:
                  fc.sendpos(r["atoms"], r["cell"])
                  r["status"]="Running"
                  fc.poll()
                  self.jobs.append([r,fc])
                  break # we are done for this request
               else: print "something very weird is going on with a client: status is:",fc.status," - trying again later"
      self.time_distribute+=time.time()-start

   # threading and signaling machinery
   # handles sigint gracefully (terminates interface, etc)
   def _kill_handler(self, signal, frame):
      print "kill called"
      self.end_thread()
      try:    self.__del__()
      except: pass      
      if signal in self._prev_kill: self._prev_kill[signal](signal, frame)
      
   def _poll_loop(self):   
      poll_iter=0
      while self._poll_true:
         time.sleep(self.latency)
         if poll_iter>UPDATEFREQ:
            self.pool_update(); poll_iter=0
         poll_iter+=1
         self.pool_distribute()
      self._poll_thread=None   
   
   def start_thread(self):
      if not self._poll_thread is None: raise NameError("Polling thread already started")      
      self._poll_thread=threading.Thread(target=self._poll_loop, name="poll_"+self.address)
      self._poll_thread.daemon=True
      self._poll_true=True
      self._prev_kill[signal.SIGINT]=signal.signal(signal.SIGINT, self._kill_handler)
      self._prev_kill[signal.SIGTERM]=signal.signal(signal.SIGTERM, self._kill_handler)    
      self._poll_thread.start()
   
   def end_thread(self):
      self._poll_true=False
      if not self._poll_thread is None: self._poll_thread.join()
      self._poll_thread=None
   
   def __del__(self):
      print "shutting down interface"
      self.server.shutdown(socket.SHUT_RDWR); self.server.close()
      if self.mode=="unix": os.unlink("/tmp/wrappi_"+self.address)                 
