"""Deals with creating the forcefield class.

Classes:
   InputForces: Deals with creating all the forcefield objects.
   InputForceBeads: Base class to deal with one particular forcefield object.
   InputFBSocket: Deals with creating a forcefield using sockets.
"""

__all__ = ['InputForces', 'InputForceBeads', "InputFBSocket"]

from engine.forces import *
from inputs.interface import InputInterface
from utils.inputvalue import *
from copy import copy

class InputForceBeads(Input):
   """ForceBeads input class.

   Handles generating one instance of a particular forcefield class from the xml
   input file, and generating the xml checkpoint tags and data from an
   instance of the object.

   Attributes:
      nbeads: The number of beads that the forcefield will be evaluated on.
      weight: A scaling factor for the contribution from this forcefield.
   """

   attribs = { "nbeads" : ( InputValue, { "dtype"   : int,
                                         "default" : 0,
                                         "help"    : "If the forcefield is to be evaluated on a contracted ring polymer, this gives the number of beads that are used. If not specified, the forcefield will be evaluated on the full ring polymer." } ),
               "weight" : ( InputValue, { "dtype"   : float,
                                         "default" : 1.0,
                                         "help"    : "This force term will be added to give the total force using this weight." } )
            }

   default_help = "Deals with the assigning of force calculation jobs and collecting the data."
   default_label = "FORCEBEADS"

   def store(self, forceb):
      """Takes a ForceBeads instance and stores a minimal representation of it.

      Args:
         forceb: A ForceBeads object.
      """

      #~ if (not type(force) is FFSocket):
         #~ raise TypeError("The type " + type(force).__name__ + " is not a valid socket forcefield")

      Input.store(self,forceb)
      self.nbeads.store(forceb.nbeads)
      self.weight.store(forceb.weight)

   def fetch(self):
      """Creates a ForceBeads object.

      Returns:
         A ForceBeads object.
      """

      super(InputForceBeads,self).fetch()

      return ForceBeads(model=ForceField(), nbeads=self.nbeads.fetch(), weight=self.weight.fetch())


class InputFBSocket(InputForceBeads, InputInterface):
   """Creates a ForceBeads object with a socket interface.

   Handles generating one instance of a socket interface forcefield class.
   Shares its attributes between InputForceBeads, which deals with creating the
   forcefield, and InputInterface, which deals with creating the socket
   interface.
   """

   attribs = copy(InputInterface.attribs)
   attribs.update(InputForceBeads.attribs)

   default_help = "Deals with the assigning of force calculation jobs to different driver codes, and collecting the data, using a socket for the data communication."
   default_label = "FORCESOCKET"

   def store(self, forceb):
      """Takes a ForceField instance and stores a minimal representation of it.

      Args:
         forceb: A ForceBeads object with a FFSocket forcemodel object.
      """

      if (not type(forceb.f_model) is FFSocket):
         raise TypeError("The type " + type(forceb.f_model).__name__ + " is not a valid socket forcefield")

      InputForceBeads.store(self,forceb)
      InputInterface.store(self,forceb.f_model.socket)

   def fetch(self):
      """Creates a ForceBeads object.

      Returns:
         A ForceBeads object with the correct socket parameters.
      """

      return ForceBeads(model=FFSocket( interface=InputInterface.fetch(self) ),nbeads=self.nbeads.fetch(),weight=self.weight.fetch() )


class InputForces(Input):
   """Deals with creating all the forcefield objects required in the 
   simulation.

   Attributes:
      extra: A list of all the forcefield objects read in dynamically from
         the xml input file.
   """

   dynamic = {  "socket" : (InputFBSocket, { "help" : InputFBSocket.default_help } )
            }

   default_help = "Deals with creating all the necessary forcefield objects."
   default_label = "MULTIFORCE"

   def fetch(self):
      """Returns a list of the output objects included in this dynamic 
      container.

      Returns:
         A list of tuples, with each tuple being of the form ('type', 'object'),
         where 'type' is the type of forcefield, and 'object' is a
      """

      super(InputForces, self).fetch()
      flist = [ (n, f.fetch()) for (n, f) in self.extra ]

      return flist

   def store(self, flist):
      """Stores a list of the output objects, creating a sequence of 
      dynamic containers.

      Args:
         A list of tuples, with each tuple being of the form ('type', 'object'),
         where 'type' is the type of forcefield, and 'object' is a
         forcefield object of that type.
      """

      super(InputForces, self).store()
      self.extra = []

      for el in flist:
         if el[0]=="socket":
            iff = InputFBSocket()
            iff.store(el[1])
            self.extra.append(("socket", iff))
