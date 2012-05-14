"""Deals with creating the multiforce class.

Classes:
   InputMulti: Deals with creating the MultiForce object from a file, and
      writing the checkpoints.
"""

__all__ = ['InputMulti']

import numpy as np
from engine.forces import *
from utils.inputvalue import *

class InputMulti(Input):
   """Multiforce input class.

   Handles generating the appropriate multiforce class from the xml
   input file, and generating the xml checkpoint tags and data from an 
   instance of the object.

   Attributes:
   """

   fields =  { "main_force" : ( InputForce, { "help" : "Main force field which calculates the forces on the full ring polymer." } ),
               "force2"  : ( InputForce, { "help"    : "Force field that calculates the force acting on a contracted ring polymer.",
                                           "default" : ForceBeads() } ),
               "force3"  : ( InputForce, { "help"    : "Force field that calculates the force acting on a contracted ring polymer.",
                                           "default" : ForceBeads() } ),
               "reduced_beads" : ( InputArray, { "help"    : "An array giving the number of beads for the contracted ring polymers.",
                                                 "dtype"   : int,
                                                 "default" : np.zeros(0, int)})}

   default_help = "Deals with ring polymer contraction, and assigning the different jobs to the different driver codes."
   default_label = "MULTIFORCE"
   
   def store(self, force):
      """Takes a multiforce instance and stores a minimal representation of it.

      Args:
         force: A multiforce object.
      """

      super(InputMulti,self).store()
      self.main_force.store(force._forces[0])
      if len(force._forces) > 1:
         self.reduced_beads.store(force.nreduced)
         self.force2.store(force._forces[1])
      if len(force._forces) > 2:
         self.force3.store(force._forces[2])

   def fetch(self):
      """Creates a multiforce object.

      Returns:
         A multiforce object of the appropriate type and with the appropriate
         interface given the attributes of the InputMulti object.
      """

      super(InputForce,self).fetch()
      forcelist = [self.main_force.fetch(), self.force2.fetch(), 
                  self.force3.fetch()]
      force = MultiForce(nreduced=self.reduced_beads.fetch(), forces=forcelist)

      return force

   def check(self):
      """Function that deals with optional arguments.

      Makes sure that the number of forcefields specified in the reduced_beads
      tag matches with the number of forcefields given.
      """

      super(InputMulti,self).check()
      nforces = 0
      if self.force2._explicit:
         nforces += 1
      if self.force3._explicit:
         nforces += 1
      if len(self.reduced_beads.fetch()) != nforces:
         raise ValueError("The number of terms given in reduced_beads does not correspond to the number of forcefields specified.")
      #TODO make sure nreduced are odd and less than nbeads
