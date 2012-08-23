"""Deals with creating the multiforce class.

Classes:
   InputMulti: Deals with creating the MultiForce object from a file, and
      writing the checkpoints.
"""

__all__ = ['InputMulti']

import numpy as np
from engine.forces import *
from inputs.forces import *
from utils.inputvalue import *

class InputMulti(Input):
   """Multiforce input class.

   Handles generating the appropriate multiforce class from the xml
   input file, and generating the xml checkpoint tags and data from an 
   instance of the object.

   Attributes:
      force: Forcefield objects which can calculate the potential and forces
         on contracted ring polymers. Multiple force objects may be specified.
   """

   dynamic =  { "force"  : ( InputForce, { "help"    : "Force field that calculates the force acting on a contracted ring polymer. More than one force field can be specified using the same input fields. If the 'nreduced' attribute is not specified then the forces will be evaluated on the full ring polymer.",
                                           "default" : input_default(factory=FFSocket) } )}

   default_help = "Deals with ring polymer contraction, and assigning the different jobs to the different driver codes."
   default_label = "MULTIFORCE"
   
   def store(self, mforce):
      """Takes a multiforce instance and stores a minimal representation of it.

      Args:
         force: A multiforce object.
      """

      super(InputMulti,self).store()
      self.extra = []

      for ff in mforce.forcelist:
         iff = InputForce()
         iff.store(ff)
         self.extra.append(("force", iff))

   def fetch(self):
      """Creates a multiforce object.

      Returns:
         A multiforce object of the appropriate type and with the appropriate
         interface given the attributes of the InputMulti object.
      """

      super(InputMulti,self).fetch()
      forcelist = [ff.fetch() for (f, ff) in self.extra]

      return MultiForce(forcelist)

   def check(self):
      """Function that deals with optional arguments.

      Makes sure that the number of forcefields specified in the reduced_beads
      tag matches with the number of forcefields given.
      """

      super(InputMulti,self).check()
      #TODO make sure the reduced ring polymer sizes are odd and less than 
      #nbeads. This will need to be done after we know what nbeads is, so
      #probably in ensemble.
