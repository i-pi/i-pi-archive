"""Deals with creating the forcefield class.

Classes:
   InputForce: Deals with creating the ForceField object from a file, and
      writing the checkpoints.
"""

__all__ = ['InputForces', 'InputForceBeads', "InputFBSocket"]

from engine.forces import *
from inputs.interface import InputInterface
from utils.inputvalue import *
from copy import copy

class InputForceBeads(Input):
   """ForceBeads input class.

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

   attribs = { "nbeads" : ( InputValue, { "dtype"   : int,
                                         "default" : 0,
                                         "help"    : "If the forcefield is to be evaluated on a contracted ring polymer, this gives the number of beads that are used. If not specified, the forcefield will be evaluated on the full ring polymer." } ),
               "weight" : ( InputValue, { "dtype"   : float,
                                         "default" : 1.0,
                                         "help"    : "This force term will be added to give the total force using this weight." } )
            }

   default_help = "Deals with the assigning of jobs to different driver codes, and collecting the data."
   default_label = "FORCEBEADS"

   def store(self, forceb):
      """Takes a ForceField instance and stores a minimal representation of it.

      Args:
         force: A forcefield object.
      """

      #~ if (not type(force) is FFSocket):
         #~ raise TypeError("The type " + type(force).__name__ + " is not a valid socket forcefield")

      Input.store(self,forceb)
      self.nbeads.store(forceb.nbeads)
      self.weight.store(forceb.weight)

   def fetch(self):
      """Creates a forcefield object.

      Returns:
         A forcefield object of the appropriate type and with the appropriate
         interface given the attributes of the InputForce object.
      """

      super(InputForceBeads,self).fetch()

      return ForceBeads(model=ForceField(), nbeads=self.nbeads.fetch(), weight=self.weight.fetch())


class InputFBSocket(InputForceBeads, InputInterface):

   attribs = copy(InputInterface.attribs)
   attribs.update(InputForceBeads.attribs)


   def store(self, forceb):
      """Takes a ForceField instance and stores a minimal representation of it.

      Args:
         force: A forcefield object.
      """

      if (not type(forceb.f_model) is FFSocket):
         raise TypeError("The type " + type(forceb.f_model).__name__ + " is not a valid socket forcefield")

      InputForceBeads.store(self,forceb)
      InputInterface.store(self,forceb.f_model.socket)

   def fetch(self):

      return ForceBeads(model=FFSocket( interface=InputInterface.fetch(self) ),nbeads=self.nbeads.fetch(),weight=self.weight.fetch() )



class InputForces(Input):

   dynamic = {  "socket" : (InputFBSocket, { "help" : "Each of the <properties> tags specify how to create a file in which one or more properties are written, one line per frame. " } )
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

      for el in flist:
         if el[0]=="socket":
            iff = InputFBSocket()
            iff.store(el[1])
            self.extra.append(("socket", iff))
