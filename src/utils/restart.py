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
      """Initialises Restart.

      Automatically adds all the fields and attribs names to the objects
      dictionary, along with a default initial value, if specified.
      """

      for f, v in self.fields.iteritems():
         self.__dict__[f] = v[0](*v[1])
      for a, v in self.attribs.iteritems():
         self.__dict__[a] = v[0](*v[1])     
      
   def store(self, value):
      """Dummy function for storing data."""

      pass
            
   def fetch(self):
      """Dummy function to retrieve data."""

      if not self.check():
         raise AssertionError("Consistency check failed in Restart.fetch")
            
   def check(self):
      """Dummy function to check for input errors."""

      return True
   
   def write(self, name="", indent=""): 
      """Writes data to the xml file.

      Writes the tag, attributes, data and closing tag appropriate to the
      particular fields and attribs data. Writes in a recursive manner, so
      that objects contained in the fields dictionary have their write function
      called, so that their tags are written between the start and end tags 
      of this object, as is required for the xml format.
      
      This also adds an indent to the lower levels of the xml heirarchy, 
      so that it is easy to see which tags contain other tags.

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
      rstr += indent + "</" + name + ">\n"
      return rstr
   
   def parse(self, xml=None, text=""): 
      """Parses an xml file.

      Uses the xml_node class defined in io_xml to read all the information 
      contained within the root tags, and uses it to give values for the attribs
      and fields data recursively. It does this by giving all the data between 
      the appropriate field tag to the appropriate field restart object as a 
      string, and the appropriate attribute data to the appropriate attribs
      restart object as a string. These data are then parsed by these objects
      until all the information is read, or an input error is found. 

      Args:
         xml: An xml_node object containing the all the data for the parent
            tag.
         text: The data held between the start and end tags.

      Raises:
         To be implemented.
      """

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
   """Scalar class for restart handling.

   Has the methods for dealing with simple data tags of the form:
   <tag_name> data </tag_name>, where data is just a value. Takes the data and
   converts it to the required data_type, so that it can be used in the 
   simulation.

   Attributes:
      type: Data type of the data.
      value: Value of data. Also specifies data type if type is None.
      default: Default value if the tag is not specified.
   """
 
   def __init__(self, dtype=None, value=None, default=None):      
      """Initialises RestartValue.

      Args:
         dtype: An optional data type. Defaults to None.
         value: An optional value for the data. Defaults to None. Also
            specifies the data type if type is None.
         default: Optional default value if the tag is not specified. Defaults
            to None.
      """

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
      """Converts the data to the appropriate data type and stores it.

      Args:
         value: The raw data to be stored.
      """

      self.value = self.type(value)
      
   def fetch(self): 
      """Returns the stored data."""

      return self.value
      
   def write(self, name="", indent=""): 
      """Writes data to the xml file.

      Writes the data in the appropriate format between appropriate tags.

      Args:
         name: An optional string giving the tag name. Defaults to "".
         indent: An optional string giving the string to be added to the start
            of the line, so usually a number of tabs. Defaults to "".

      Returns:
         A string giving the stored value in the appropriate xml format.
      """

      return indent + "<" + name + ">" + write_type(self.type, self.value) + "</" + name + ">\n"
   
   def parse(self, xml=None, text=""):
      """Reads the data for a single value from an xml file.

      Args:
         xml: An xml_node object containing the all the data for the parent
            tag.
         text: The data held between the start and end tags.
      """

      if xml is None:
         self.value = read_type(self.type, text)
      else:   
         self.value = read_type(self.type, xml.fields["_text"])

         
ELPERLINE = 5
class RestartArray(Restart):
   """Array class for restart handling.

   Has the methods for dealing with simple data tags of the form:
   <tag_name shape="(shape)"> data </tag_name>, where data is just an array and    of the form [data[0], data[1], ... , data[length]]. 
   
   Takes the data and converts it to the required data_type, 
   so that it can be used in the simulation. Also holds the shape of the array,
   so that we can use a simple 1D list of data to specify a multi-dimensional
   array.

   Attributes:
      type: Data type of the data.
      value: Value of data. Also specifies data type if type is None.
      default: Default value if the tag is not specified.
   """
 
   attribs = { "shape" : (RestartValue,(tuple, ())) }
   def __init__(self, dtype=None, value=None, default=None):
      """Initialises RestartArray.

      Args:
         dtype: An optional data type. Defaults to None.
         value: An optional value for the data. Defaults to None. Also
            specifies the data type if type is None.
         default: Optional default value if the tag is not specified. Defaults
            to None.
      """

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
      """Converts the data to the appropriate data type and shape and stores it.

      Args:
         value: The raw data to be stored.
      """

      self.shape.store(value.shape)
      self.value = np.array(value, dtype=self.type).flatten().copy()
      
   def fetch(self): 
      """Fetches the stored data.

      If the shape has not been specified, it assumes a 1D array of the 
      same length as the data list.
      """

      if self.shape.fetch() == (0,):
         return np.resize(self.value,0).copy()
      else:
         return self.value.reshape(self.shape.fetch()).copy()

   def write(self, name="", indent=""): 
      """Writes data to the xml file.

      Writes the data in the appropriate format between appropriate tags. Note
      that only ELPERLINE values are printed on each line if there are more than
      this in the array. If the values are floats, or another data type with
      a fixed width of data output, then they are aligned in columns.

      Args:
         name: An optional string giving the tag name. Defaults to "".
         indent: An optional string giving the string to be added to the start
            of the line, so usually a number of tabs. Defaults to "".

      Returns:
         A string giving the stored value in the appropriate xml format.
      """

      rstr = indent + "<" + name + " shape='" + write_tuple(self.shape.fetch()) + "'> "
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
      """Reads the data for a single array from an xml file.

      Args:
         xml: An xml_node object containing the all the data for the parent
            tag.
         text: The data held between the start and end tags.
      """

      if xml is None:
         pass #TODO should print an error!
      else:   
         self.value = read_array(self.type,xml.fields["_text"]) 
         if "shape" in xml.attribs: 
            self.shape.store(read_type(tuple,xml.attribs["shape"]))
         else:
            self.shape.store(self.value.shape)
      
