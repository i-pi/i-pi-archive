"""Creates objects that deal with the evaluation of interactions."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


from copy import copy

from ipi.engine.forcefields import ForceField, FFSocket, FFLennardJones
from ipi.interfaces.sockets import InterfaceSocket
from ipi.utils.inputvalue import *


__all__ = ["InputFFSocket", 'InputFFLennardJones']


class InputForceField(Input):
   """ForceField input class.

   Handles generating one instance of a particular forcefield class from the xml
   input file, and generating the xml checkpoint tags and data from an
   instance of the object.

   Attributes:
      name: The name by which the forcefield will be identified in the System forces section.
      pbc: A boolean describing whether periodic boundary conditions will
         be applied to the atom positions before they are sent to the driver
         code.

   Fields:
      latency: The number of seconds to sleep between looping over the requests.
      parameters: A dictionary containing the forcefield parameters.
   """

   attribs = { "name" : ( InputAttribute, { "dtype"   : str,
                                         "help"    : "Mandatory. The name by which the forcefield will be identified in the System forces section." } ),
               "pbc":  ( InputAttribute, { "dtype"   : bool,
                                         "default" : True,
                                         "help"    : "Applies periodic boundary conditions to the atoms coordinates before passing them on to the driver code." })

            }
   fields = {
            "latency" : ( InputValue, { "dtype"   : float,
                                         "default" : 0.01,
                                         "help"    : "The number of seconds the polling thread will wait between exhamining the list of requests." } ),
            "parameters" : (InputValue, { "dtype" : dict,
                                     "default" : {},
                                     "help" : "The parameters of the force field"} )
   }

   default_help = "Base forcefield class that deals with the assigning of force calculation jobs and collecting the data."
   default_label = "FORCEFIELD"

   def store(self, ff):
      """Takes a ForceField instance and stores a minimal representation of it.

      Args:
         forceb: A ForceField object.
      """

      Input.store(self,ff)
      self.name.store(ff.name)
      self.latency.store(ff.latency)
      self.parameters.store(ff.pars)
      self.pbc.store(ff.dopbc)

   def fetch(self):
      """Creates a ForceField object.

      Returns:
         A ForceField object.
      """

      super(InputForceField,self).fetch()

      return ForceField(pars = self.parameters.fetch(), name = self.name.fetch(), latency = self.latency.fetch(), dopbc = self.pbc.fetch())


class InputFFSocket(InputForceField):
   """Creates a ForceField object with a socket interface.

   Handles generating one instance of a socket interface forcefield class.

   Attributes:
      mode: Describes whether the socket will be a unix or an internet socket.

   Fields:
      address: The server socket binding address.
      port: The port number for the socket.
      slots: The number of clients that can queue for connections at any one
         time.
      timeout: The number of seconds that the socket will wait before assuming
         that the client code has died. If 0 there is no timeout.
   """

   fields = {"address": (InputValue, {"dtype"   : str,
                                      "default" : "localhost",
                                      "help"    : "This gives the server address that the socket will run on." } ),
             "port":    (InputValue, {"dtype"   : int,
                                      "default" : 65535,
                                      "help"    : "This gives the port number that defines the socket."} ),
             "slots":   (InputValue, {"dtype"   : int,
                                      "default" : 4,
                                      "help"    : "This gives the number of client codes that can queue at any one time."} ),
             "timeout": (InputValue, {"dtype"   : float,
                                      "default" : 0.0,
                                      "help"    : "This gives the number of seconds before assuming a calculation has died. If 0 there is no timeout." } )}
   attribs = {
               "mode": (InputAttribute, {"dtype"    : str,
                                     "options"  : [ "unix", "inet" ],
                                     "default"  : "inet",
                                     "help"     : "Specifies whether the driver interface will listen onto a internet socket [inet] or onto a unix socket [unix]." } ),

              }

   attribs.update(InputForceField.attribs)
   fields.update(InputForceField.fields)

   default_help = "Deals with the assigning of force calculation jobs to different driver codes, and collecting the data, using a socket for the data communication."
   default_label = "FFSOCKET"

   def store(self, ff):
      """Takes a ForceField instance and stores a minimal representation of it.

      Args:
         ff: A ForceField object with a FFSocket forcemodel object.
      """

      if (not type(ff) is FFSocket):
         raise TypeError("The type " + type(ff).__name__ + " is not a valid socket forcefield")


      super(InputFFSocket,self).store(ff)

      self.address.store(ff.socket.address)
      self.port.store(ff.socket.port)
      self.timeout.store(ff.socket.timeout)
      self.slots.store(ff.socket.slots)
      self.mode.store(ff.socket.mode)

   def fetch(self):
      """Creates a ForceSocket object.

      Returns:
         A ForceSocket object with the correct socket parameters.
      """

      return FFSocket(pars = self.parameters.fetch(), name = self.name.fetch(), latency = self.latency.fetch(), dopbc = self.pbc.fetch(),
              interface=InterfaceSocket(address=self.address.fetch(), port=self.port.fetch(),
            slots=self.slots.fetch(), mode=self.mode.fetch(), timeout=self.timeout.fetch() ) )


   def check(self):
      """Deals with optional parameters."""

      super(InputFFSocket,self).check()
      if self.port.fetch() < 1 or self.port.fetch() > 65535:
         raise ValueError("Port number " + str(self.port.fetch()) + " out of acceptable range.")
      elif self.port.fetch() < 1025:
         warning("Low port number being used, this may interrupt important system processes.", verbosity.low)

      if self.slots.fetch() < 1 or self.slots.fetch() > 5:
         raise ValueError("Slot number " + str(self.slots.fetch()) + " out of acceptable range.")
      if self.latency.fetch() < 0:
         raise ValueError("Negative latency parameter specified.")
      if self.timeout.fetch() < 0.0:
         raise ValueError("Negative timeout parameter specified.")


class InputFFLennardJones(InputForceField):

   attribs = {}
   attribs.update(InputForceField.attribs)

   default_help = """Simple, internal LJ evaluator without cutoff, neighbour lists or minimal image convention.
                   Expects standard LJ parameters, e.g. { eps: 0.1, sigma: 1.0 }. """
   default_label = "FFLJ"

   def store(self, ff):
      super(InputFFLennardJones,self).store(ff)

   def fetch(self):
      super(InputFFLennardJones,self).fetch()

      return FFLennardJones(pars = self.parameters.fetch(), name = self.name.fetch(),
               latency = self.latency.fetch(), dopbc = self.pbc.fetch())
