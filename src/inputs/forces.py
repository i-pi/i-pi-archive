"""Deals with creating the forcefield class.

Classes:
   InputForce: Deals with creating the ForceField object from a file, and
      writing the checkpoints.
"""

__all__ = ['InputForces']

from engine.forces import *
from inputs.interface import InputInterface
from utils.inputvalue import *
from copy import copy

class InputFFSocket(InputInterface):
   """FFSocket input class.

   Handles generating the appropriate forcefield class from the xml
   input file, and generating the xml checkpoint tags and data from an
   instance of the object.

   Attributes:
      type: A string indicating the type being used. 'socket' is currently
         the only available option.
      interface: A restart interface instance.
      parameters: A dictionary of the parameters used by the driver. Of the
         form {"name": value}.
      nreduced: An integer giving the number of beads to use if a ring polymer
         contraction scheme is being used.
   """

   attribs = copy(InputInterface.attribs)
   attribs["nbeads"] = ( InputValue, { "dtype"   : int,
                                         "default" : 0,
                                         "help"    : "If the forcefield is to be evaluated on a contracted ring polymer, this gives the number of beads that are used. If not specified, the forcefield will be evaluated on the full ring polymer." } )
   attribs["weight"] = ( InputValue, { "dtype"   : float,
                                         "default" : 1.0,
                                         "help"    : "This force term will be added to give the total force using this weight." } )

   default_help = "Deals with the assigning of jobs to different driver codes, and collecting the data."
   default_label = "FORCES"

   def store(self, force):
      """Takes a ForceField instance and stores a minimal representation of it.

      Args:
         force: A forcefield object.
      """

      if (not type(force) is FFSocket):
         raise TypeError("The type " + type(force).__name__ + " is not a valid socket forcefield")

      super(InputFFSocket,self).store(force.socket)
      self.nbeads.store(force.nbeads)
      self.weight.store(force.weight)

   def fetch(self):
      """Creates a forcefield object.

      Returns:
         A forcefield object of the appropriate type and with the appropriate
         interface given the attributes of the InputForce object.
      """

      interface=super(InputFFSocket,self).fetch()
      force = FFSocket( interface=interface,nbeads=self.nbeads.fetch(),weight=self.weight.fetch() )

      print "reading force", self.nbeads.fetch(), " ", self.weight.fetch()
      return force


class InputForces(Input):

   dynamic = {  "socket" : (InputFFSocket, { "help" : "Each of the <properties> tags specify how to create a file in which one or more properties are written, one line per frame. " } )
            }

   def fetch(self):
      """ Returs a list of the output objects included in this dynamic container. """

      super(InputForces, self).fetch()
      flist = [ (n, f.fetch()) for (n, f) in self.extra ]

      return flist

   def store(self, flist):
      """ Stores a list of the output objects, creating a sequence of dynamic containers. """

      super(InputForces, self).store()
      self.extra = []

      print "storing forces", flist
      for el in flist:
         if el[0]=="socket":
            iff = InputFFSocket()
            iff.store(el[1])
            self.extra.append(("socket", iff))
