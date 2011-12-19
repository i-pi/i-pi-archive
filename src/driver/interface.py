import socket, select, string
from utils.restart import Restart, RestartValue
import numpy as np
import time, pdb

HDRLEN=12
TIMEOUT=10.0
SERVERTIMEOUT=2.0*TIMEOUT

def Message(mystr): return string.ljust(string.upper(mystr),HDRLEN)

class Disconnected(Exception): pass
class InvalidStatus(Exception): pass
class Status: 
   Disconnected=0;   Up = 1;  Ready=2;   NeedsInit=4;  HasData=8;  Busy=16;
   
class Driver(socket.socket):
   def __init__(self, socket):
      super(Driver,self).__init__(_sock=socket)
      self._buf=np.zeros(0,np.byte)
      
   def status(self):
      try:
         self.sendall(Message("status"))      
      except:
         return Status.Disconnected

      try:
         reply=self.recv(HDRLEN)
      except socket.timeout:
         return Status.Up | Status.Busy
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
         return Status.Up
   
   def recvall(self, dest):      
      blen=dest.itemsize*dest.size
      if (blen>len(self._buf)) : self._buf.resize(blen)  # keeps a permanent buffer, which is expanded if necessary
      bpos=0
      while bpos<blen:
         try:
            bpart=self.recv_into(self._buf[bpos:], blen-bpos )
         except socket.timeout:
            print "timeout, v2 trying again"; pass
         if (bpart == 0 ): raise Disconnected()
         bpos+=bpart

      if np.isscalar(dest):       
         return np.fromstring(self._buf[0:blen], dest.dtype)[0]
      else:
         return np.fromstring(self._buf[0:blen], dest.dtype).reshape(dest.shape)
      
       
   def initialize(self, pars):
      self.sendall(Message("init"));
      self.sendall(np.int32(len(pars)))
      self.sendall(pars)                        
               
   def sendpos(self, atoms, cell):
      if (self.status() & Status.Ready):
         self.sendall(Message("posdata"));
         self.sendall(cell.h, 9*8)
         self.sendall(cell.ih,9*8)      
         self.sendall(np.int32(len(atoms)))
         self.sendall(atoms.q,len(atoms)*3*8)
      else: print "status was",self.status(); raise InvalidStatus()

   def getforce(self):
      if (self.status() & Status.HasData):
         self.sendall(Message("getforce"));
         while True:
            try:
               reply=self.recv(HDRLEN)
            except socket.timeout: pass
            if reply==Message("forceready"): break
            if reply=="": raise Disconnected()
            
         mu=np.float64(); mu=self.recvall(mu)         
         mlen=np.int32(); mlen=self.recvall(mlen)

         mf=np.zeros(3*mlen,np.float64)
         mf=self.recvall(mf)

         mvir=np.zeros((3,3),np.float64)
         mvir=self.recvall(mvir)
         return  [ mu, mf, mvir ]
         
      else: raise InvalidStatus()
         
         
class RestartInterface(Restart):         
   fields={ "address" : (RestartValue, (str, "localhost")), "port" : (RestartValue, (int,31415)),
            "slots" : (RestartValue, (int, 2) ) }
   def store(self, iface):
      self.address.store(iface.address)
      self.port.store(iface.port)
      self.slots.store(iface.slots)
      
   def fetch(self):
      return Interface(address=self.address.fetch(), port=self.port.fetch(), slots=self.slots.fetch())
            
class Interface(object):
   def __init__(self, address="localhost", port=3141, slots=1):
      self.server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
      self.address = address; self.port = port; self.slots = slots
      self.server.bind((address,port))
      self.server.listen(slots)
      self.server.settimeout(SERVERTIMEOUT)
      self.clients=[]
      self.requests=[]
      self.jobs=[]      
      
   def __del__(self):
      print "shutting down interface"
      self.server.shutdown(socket.SHUT_RDWR); self.server.close()
      
   def queue(self, atoms,cell, pars=""):
      newreq={"atoms":atoms,"cell":cell, "pars":pars, "result":None, "status":"Queued"}
      self.requests.append(newreq);
      return newreq
      
   def release(self, request):
      if request in self.requests:         
         self.requests.remove(request)
         
   def pool_update(self):
      # first, remove clients which had disconnected
      for c in self.clients[:]:
         if not (c.status() & Status.Up):
            try:
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
            if (driver.status() | Status.Up):
               self.clients.append(driver)
               print "handshake completed"
            else:
               print "handshake failed"
               client.shutdown(socket.SHUT_RDWR); client.close()
         else: keepsearch=False
         
   def pool_distribute(self):
      # first, checks if jobs have been completed
      for [r,c] in self.jobs[:]:
         if c.status() & Status.HasData:
            r["result"]=c.getforce()
            r["status"]="Done"
            self.jobs.remove([r,c])
                     
      for r in self.requests:
         if r["status"] == "Queued":
            # gets a free client (if any)
            freec=self.clients[:]
            for [r2, c] in self.jobs:
               freec.remove(c)
               
            if len(freec) > 0 : 
               fc=freec[0]
               st=fc.status()
               if not (st & (Status.Ready | Status.NeedsInit) ):
                  print "something weird is going on with a client"
               else:
                  if st & Status.NeedsInit:
                     fc.initialize(r["pars"])                              
                  fc.sendpos(r["atoms"], r["cell"])
                  r["status"]="Running"
                  self.jobs.append([r,fc])
                  
