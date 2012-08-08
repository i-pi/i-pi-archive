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
   """

   attribs={"nfields"  : (InputValue, {"dtype"   : int,
                                       "help"    : "The total number of forcefields that will be used. If not specified, then it will assume that the last number 'x' for which 'forcex' is specified is the correct value. Note that this will mean that any error in the number or name of the forcefields will not be checked.",
                                       "default" : 1}) }
   fields =  { "force1"  : ( InputForce, { "help"    : "Force field that calculates the force acting on a contracted ring polymer. If more than one force field is needed, you can specify 'force2', 'force3', etc. and the same input fields will work. If the 'type' attribute is 'socket', then the forces will be evaluated on the full ring polymer.",
                                           "default" : FFSocket() } )}

   default_help = "Deals with ring polymer contraction, and assigning the different jobs to the different driver codes."
   default_label = "MULTIFORCE"
   
   def parse(self, xml=None, text=""): 
      """Parses an xml file.

      Overwrites the standard parse function so that a variable number of 
      forcefields can be used, all with the same input interface.

      Args:
         xml: An xml_node object containing the all the data for the parent
            tag.
         text: The data held between the start and end tags.

      Raises:
         NameError: Raised if one of the tags in the xml input file is 
            incorrect.
         ValueError: Raised if the user does not specify a required field.
      """

      self._explicit = True
      if not xml is None:
         for a, v in xml.attribs.iteritems() :
            if a in self.attribs: 
               self.__dict__[a].parse(text=v)
            elif a == "_text":
               pass
            else:
               raise NameError("Attribute name '" + a + "' is not a recognized property of '" + xml.name + "' objects")
         
         self.adapt()      
         
         for f, v in xml.fields.iteritems():
            if f in self.fields:
               self.__dict__[f].parse(xml=v)
            elif f == "_text":
               pass
            elif f[0:5] == "force": #Creates new tags as necessary
               try:
                  force_no = int(f[5:])
               except ValueError:
                  raise ValueError("Non-integer specifier for one of the force tags in multiforce.")
               force1_v = self.fields["force1"]
               self.__dict__[f] = force1_v[0](**force1_v[1]) 
               self.__dict__[f].parse(text=v)
            else:
               raise NameError("Tag name '" + f + "' is not a recognized property of '" + xml.name + "' objects")
               
      for a in self.attribs:
         va = self.__dict__[a]
         if not (va._explicit or va._optional):
            raise ValueError("Attribute name '" + a + "' is mandatory and was not found in the input for the property " + xml.name)
      for f in self.fields:
         vf = self.__dict__[f]
         if not (vf._explicit or vf._optional):
            raise ValueError("Field name '" + f + "' is mandatory and was not found in the input for the property " + xml.name)

   def store(self, force):
      """Takes a multiforce instance and stores a minimal representation of it.

      Args:
         force: A multiforce object.
      """

      super(InputMulti,self).store()
      self.nfields.store(len(force._forces))
      for ff in range(len(force._forces)):
         force_name = "force" + str(ff+1)
         self.__dict__[force_name].store(force._forces[ff])

   def write(self, name="", indent=""): 
      """Writes data in xml file format.

      Overwrites the standard write function so that a variable number of 
      forcefields can be used, all with the same input interface.

      Args:
         name: An optional string giving the tag name. Defaults to "".
         indent: An optional string giving the string to be added to the start
            of the line, so usually a number of tabs. Defaults to "".

      Returns:
         A string giving all the data contained in the fields and attribs
         dictionaries, in the appropriate xml format.
      """

      rstr = indent + "<" + name;
      for a in self.attribs:      
         rstr += " " + a + "='" + str(self.__dict__[a].fetch()) + "'"
      rstr += ">\n"
      for f in self.fields:
         rstr += self.__dict__[f].write(f, "   " + indent)
      for ff in range(self.nfields.fetch()-1):
         force_name = "force" + str(ff+2)
         rstr += self.__dict__[force_name].write(force_name, "   " + indent)
      rstr += indent + "</" + name + ">\n"
      return rstr

   def fetch(self):
      """Creates a multiforce object.

      Returns:
         A multiforce object of the appropriate type and with the appropriate
         interface given the attributes of the InputMulti object.
      """

      super(InputMulti,self).fetch()
      forcelist = []
      force_no = 1
      while (True)
         try:
            ff_name = "force" + str(force_no)
            forcelist.append(self.__dict__[ff_name].fetch())
            force_no += 1
         except KeyError:
            if not self.nfields._explicit:
               self.nfields.store(len(forcelist))
               print "'nfields' tag in 'multiforce' not given, and has been set to ", len(forcelist)
            elif len(forcelist) != self.nfields.fetch():
               raise ValueError("The number of forcefields does not match the 'nfields' attribute")
            elif len(forcelist) == 0:
               raise ValueError("'force1' tag not specified in multiforce.")
            break
      force = MultiForce(forces=forcelist)

      return force

   def check(self):
      """Function that deals with optional arguments.

      Makes sure that the number of forcefields specified in the reduced_beads
      tag matches with the number of forcefields given.
      """

      super(InputMulti,self).check()
      #TODO make sure the reduced ring polymer sizes are odd and less than 
      #nbeads. This will need to be done after we know what nbeads is, so
      #probably in ensemble.
