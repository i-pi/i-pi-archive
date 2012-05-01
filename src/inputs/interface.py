"""Deals with creating the interface between the wrapper and the socket.

Classes:
   InputInterface: Deals with creating the Interface object from a file, and
      writing the checkpoints.
"""

__all__ = [ 'InputInterface' ]

import socket, select, threading, signal, string, os, time
import numpy as np
from utils.inputvalue import Input, InputValue
from driver.interface import *

class InputInterface(Input):         
   """Interface input class.

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
      timeout: A float giving a number of seconds after which a calculation core
         is considered dead. Defaults to zero, i.e. no timeout.
   """

   fields = {"address": (InputValue, {"dtype"   : str,
                                      "default" : "localhost", 
                                      "help"    : "This gives the server address that the socket will run on." } ), 
             "port":    (InputValue, {"dtype"   : int, 
                                      "default" : 31415, 
                                      "help"    : "This gives the port number that defines the socket."} ),
             "slots":   (InputValue, {"dtype"   : int, 
                                      "default" : 4, 
                                      "help"    : "This gives the number of client codes that can queue at any one time."} ), 
             "latency": (InputValue, {"dtype"   : float, 
                                      "default" : 1e-3,
                                      "help"    : "This gives the number of seconds between each check for new clients."} ),
             "timeout": (InputValue, {"dtype"   : float, 
                                      "default" : 0.0, 
                                      "help"    : "This gives the number of seconds before assuming a calculation has died. If 0 there is no timeout." } )}
   attribs = { "mode": (InputValue, {"dtype"    : str,
                                     "options"  : [ "unix", "inet" ],
                                     "default"  : "inet", 
                                     "help"     : "Specifies whether the driver interface will listen onto a internet socket [inet] or onto a unix socket [unix]." } )}

   default_help = "Specifies the parameters for the socket interface."
   default_label = "INTERFACE"

   def store(self, iface):
      """Takes an Interface instance and stores a minimal representation of it.

      Args:
         iface: An interface object.
      """

      super(InputInterface,self).store(iface)
      self.latency.store(iface.latency)
      self.mode.store(iface.mode)
      self.address.store(iface.address)
      self.port.store(iface.port)
      self.slots.store(iface.slots)
      self.timeout.store(iface.timeout)
      
   def fetch(self):
      """Creates an Interface object.

      Returns:
         An interface object with the appropriate socket given the attributes
         of the InputInterface object.
      """

      super(InputInterface,self).fetch()
      return Interface(address=self.address.fetch(), port=self.port.fetch(), 
            slots=self.slots.fetch(), mode=self.mode.fetch(), 
            latency=self.latency.fetch(), timeout=self.timeout.fetch())

   def check(self):
      """Function that deals with optional arguments."""

      super(InputInterface,self).check()
      if self.port.fetch() < 1 or self.port.fetch() > 65535:
         raise ValueError("Port number out of acceptable range.")
      elif self.port.fetch() < 1025: 
         print "Warning, low port number being used, this may interrupt important system processes."
      
      if self.slots.fetch() < 1 or self.slots.fetch() > 5:
         raise ValueError("Slot number out of acceptable range.")
      if self.latency.fetch() < 0:
         raise ValueError("Negative latency parameter specified.")
      if self.timeout.fetch() < 0.0:
         raise ValueError("Negative timeout parameter specified.")
