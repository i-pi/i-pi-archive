"""Deals with creating the forcefield class.

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


Classes:
   InputForces: Deals with creating all the forcefield objects.
   InputForceBeads: Base class to deal with one particular forcefield object.
   InputFBSocket: Deals with creating a forcefield using sockets.
"""

__all__ = ['InputForces', 'InputForceComponent']

from copy import copy
from ipi.engine.forces import *
from ipi.inputs.interface import InputInterfaceSocket
from ipi.utils.inputvalue import *

class InputForceComponent(InputValue):
   """ForceComponent input class.

   Uses the forcefield object whose name is specified as the value of the 
   field (matching one of the forcefields defined in the simulation tag)
   to compute one component of the force acting on the ring polymer.
      

   Attributes:
      nbeads: The number of beads that the forcefield will be evaluated on.
      weight: A scaling factor for the contribution from this forcefield.
   """

   attribs = { "nbeads" : ( InputAttribute, { "dtype"   : int,
                                         "default" : 0,
                                         "help"    : "If the forcefield is to be evaluated on a contracted ring polymer, this gives the number of beads that are used. If not specified, the forcefield will be evaluated on the full ring polymer." } ),
               "weight" : ( InputAttribute, { "dtype"   : float,
                                         "default" : 1.0,
                                         "help"    : "A scaling factor for this forcefield, to be applied before adding the force calculated by this forcefield to the total force." } ),
               "name" : ( InputAttribute, { "dtype" : str,
                                          "default" : "",
                                          "help" : "An optional name to refer to this force component." } )
            }

   default_help = "Base class that deals with the assigning of force calculation jobs and collecting the data."
   default_label = "FORCECOMPONENT"
   
   def __init__(self, help=None, dimension=None, units=None, default=None, dtype=None):
      """Initializes InputCell.

      Just calls the parent initialization function with appropriate arguments.
      """

      super(InputForceComponent,self).__init__(dtype=str, dimension=dimension, default=default, help=help)
      
   def store(self, forceb):
      """Takes a ForceBeads instance and stores a minimal representation of it.

      Args:
         forceb: A ForceBeads object.
      """

      super(InputForceComponent,self).store(forceb.ffield)
      self.nbeads.store(forceb.nbeads)
      self.weight.store(forceb.weight)
      self.name.store(forceb.name)

   def fetch(self):
      """Creates a ForceBeads object.

      Returns:
         A ForceBeads object.
      """

      val=super(InputForceComponent,self).fetch()      
      return ForceComponent(ffield=val, nbeads=self.nbeads.fetch(), weight=self.weight.fetch(), name=self.name.fetch())

   def check(self):
      """Checks for optional parameters."""

      super(InputForceComponent,self).check()
      if self.nbeads.fetch() < 0:
         raise ValueError("The forces must be evaluated over a positive number of beads.")


class InputForces(Input):
   """Deals with creating all the forcefield objects required in the
   simulation.

   Dynamic fields:
      socket: Socket object to create the server socket.
   """

   #At the moment only socket driver codes implemented, other types
   #could be used in principle
   dynamic = {  "force" : (InputForceComponent, { "help" : InputForceComponent.default_help } )
            }

   default_help = "Deals with creating all the necessary forcefield objects."
   default_label = "FORCES"

   def fetch(self):
      """Returns a list of the output objects included in this dynamic
      container.

      Returns:
         A list of tuples, with each tuple being of the form ('type', 'object'),
         where 'type' is the type of forcefield, and 'object' is a
      """

      super(InputForces, self).fetch()   
      flist = [ f.fetch() for (n, f) in self.extra ]      
      return flist

   def store(self, flist):
      """Stores a list of the output objects, creating a sequence of
      dynamic containers.

      Args:
         flist: A list of tuples, with each tuple being of the form
         ('type', 'object') where 'type' is the type of forcefield
         and 'object' is a forcefield object of that type.
      """

      super(InputForces, self).store()
      self.extra = []

      for el in flist:
         iff = InputForceComponent()
         iff.store(el)
         self.extra.append(("force", iff))
