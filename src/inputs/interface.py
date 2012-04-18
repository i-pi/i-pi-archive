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

__all__ = [ 'RestartInterface' ]

import socket, select, threading, signal, string, os, time
import numpy as np
from utils.restart import Restart, RestartValue
from utils.inputvalue import Input, InputValue
from driver.interface import *



class RestartInterface(Input):         
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
      timeout: A float giving a number of seconds after which a calculation core
         is considered dead. Defaults to zero, i.e. no timeout.
   """

   fields = { "address" : (InputValue, { "dtype" : str, "default" : "localhost" } ), 
              "port" :    (InputValue, { "dtype" : int, "default" :  31415      } ),
              "slots" :   (InputValue, { "dtype" : int, "default" :  4          } ), 
              "latency" : (InputValue, { "dtype" : float, "default" : 1e-3      } ), 
              "timeout":  (InputValue, { "dtype" : float, "default" : 0.0       } ) }
   attribs = { "mode": (InputValue, { "dtype" : str,
                                      "options" : [ "unix", "inet" ],
                                      "default" : "inet", 
                                      "help"    : "Specifies whether the driver interface will listen onto a internet socket [inet] or onto a unix socket[unix]" } ) 
             }

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
      self.timeout.store(iface.timeout)
      
   def fetch(self):
      """Creates an Interface object.

      Returns:
         An interface object with the appropriate socket given the attributes
         of the RestartInterface object.
      """

      return Interface(address=self.address.fetch(), port=self.port.fetch(), 
            slots=self.slots.fetch(), mode=self.mode.fetch(), 
            latency=self.latency.fetch(), timeout=self.timeout.fetch())
