"""Contains the classes that are used to write to and read from restart files.

The classes defined in this module define the base functions which parse the
data in the restart files. Each restart object defined has a fields and an
attributes dictionary, which are filled with the tags and attributes that 
are allowed to be present, along with their default values and data type.

These are then filled with the data from the xml file when the program
is initialised, and are filled by the values calculated in the program which 
are then output to the checkpoint file when a restart file is required.

Also deals with checking for user input errors, of the form of misspelt tags,
bad data types, and failure to input required fields.

Classes:
   Restart: Base restart class with the generic methods and attributes.
   RestartValue: Restart class for scalar objects.
   RestartArray: Restart class for arrays.
"""

__all__ = ['Restart', 'RestartValue', 'RestartArray']

import numpy as np
from  io.io_xml import *

class Restart(object):
   """Base class for restart handling.

   Has the generic methods for dealing with the restart file. Parses the input
   data, outputs the output data, and deals with storing and returning the 
   data obtained during the simulation for the restart files.

   Attributes:
      fields: A dictionary holding the possible tags contained within the 
         tags for this restart object, which are then turned into the objects
         held by the object given by this restart object. The dictionary is
         of the form:
         {"tag": (restart_type, (data_type, default_value)), ... }.
      attribs: A dictionary holding the attribute data for the tag for this
         restart object. The dictionary is of the form:
         {"attrib_name": (restart_type, (data_type, default_value)), ... }.
   """
 
   fields = {}
   attribs = {}

   def __init__(self):
      for f, v in self.fields.iteritems():
         self.__dict__[f] = v[0](*v[1])
      for a, v in self.attribs.iteritems():
         self.__dict__[a] = v[0](*v[1])     
      
   def store(self, value): 
      pass
            
   def fetch(self): 
      if not self.check():
         raise AssertionError("Consistency check failed in Restart.fetch")
            
   def check(self):
      return True
   
   def write(self, name="", indent=""): 
      rstr = indent + "<" + name;
      for a in self.attribs:      
         rstr += " " + a + "='" + str(self.__dict__[a].fetch()) + "'"
      rstr += ">\n"
      for f in self.fields:
         rstr += self.__dict__[f].write(f, "   " + indent)
      rstr += indent + "</" + name + ">\n"
      return rstr
   
   def parse(self, xml=None, text=""): 
      if not xml is None:
         for a, v in xml.attribs.iteritems() :
            if a in self.attribs: 
               self.__dict__[a].parse(text=v)
            else:
               pass #TODO raise an exception
         for f, v in xml.fields.iteritems():
            if f in self.fields:
               self.__dict__[f].parse(xml=v)
            else:
               pass #TODO raise an exception            

   
class RestartValue(Restart):
   def __init__(self, dtype=None, value=None, default=None):      
      super(RestartValue,self).__init__()
      if not dtype is None:
         self.type = dtype
      elif not self.value is None:
         self.type = type(value)
      else:
         raise TypeError("You must provide either value or dtype")

      self.default = default
      if not value is None:
         self.store(value)
      elif not default is None:
         self.store(default)
      
   def store(self, value):
      self.value = self.type(value)
      
   def fetch(self): 
      return self.value
      
   def write(self, name="", indent=""): 
      return indent + "<" + name + ">" + write_type(self.type, self.value) + "</" + name + ">\n"
   
   def parse(self, xml=None, text=""):
      if xml is None:
         self.value = read_type(self.type, text)
      else:   
         self.value = read_type(self.type, xml.fields["_text"])

         
ELPERLINE = 5
class RestartArray(Restart):
   attribs = { "shape" : (RestartValue,(tuple, ())) }
   def __init__(self, dtype=None, value=None, default=None):
      super(RestartArray,self).__init__()
      
      if not dtype is None:
         self.type = dtype
      elif not self.value is None:
         self.type = value.dtype
      else:
         raise TypeError("You must provide either value or dtype")

      self.default = default
      if not value is None:
         self.store(value)
      elif not default is None:
         self.store(default)
   
   def store(self, value):
      self.shape.store(value.shape)
      self.value = np.array(value, dtype=self.type).flatten().copy()
      
   def fetch(self): 
      if self.shape.fetch() == (0,):
         return np.resize(self.value,0).copy()
      else:
         return self.value.reshape(self.shape.fetch()).copy()

   def write(self, name="", indent=""): 
      rstr = indent + "<" + name + " shape='" + write_tuple(self.shape.fetch()) + "'> ";
      if (len(self.value)>ELPERLINE):
         rstr += "\n" + indent + " [ "
      else:
         rstr += " [ "
      for i,v in enumerate(self.value):          
         if (len(self.value) > ELPERLINE and i > 0 and i%ELPERLINE == 0):
            rstr += "\n" + indent + "   "
         rstr += write_type(self.type, v) + ", "
      rstr=rstr.rstrip(", ")
      if (len(self.value)>ELPERLINE):
         rstr += " ]\n" + indent
      else:
         rstr += " ] "
      rstr += "</" + name + ">\n"
      return rstr
      
   def parse(self, xml=None, text=""):
      if xml is None:
         pass #TODO should print an error!
      else:   
         self.value = read_array(self.type,xml.fields["_text"]) 
         if "shape" in xml.attribs: 
            self.shape.store(read_type(tuple,xml.attribs["shape"]))
         else:
            self.shape.store(self.value.shape)
      
