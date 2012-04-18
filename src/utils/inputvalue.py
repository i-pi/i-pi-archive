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
   Input: Base input class with the generic methods and attributes.
   InputValue: Input class for scalar objects.
   InputArray: Input class for arrays.
"""

__all__ = ['Input', 'InputValue', 'InputArray']

import numpy as np
from  io.io_xml import *
from units import UnitMap

class Input(object):
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
   default_help = "Generic input value"
   default_dimension = "undefined"
   default_units = ""
   default_value = None

   def __init__(self, help=None, dimension=None, units = None, default=None):
      """Initialises Input.

      Automatically adds all the fields and attribs names to the objects
      dictionary as keys, then initialises all the appropriate restart objects
      as the corresponding values.
      """

      if help is None:
         self._help = self.default_help  
      else:
         self._help = help
      if dimension is None:
         self._dimension = self.default_dimension
      else:
         self._dimension = dimension
      if units is None:
         self.units = self.default_units
      else:
         self.units = units
      if default is None:     
         self._default = self.default_value
         self._optional = False
      else:                   
         self._default = default
         self._optional = True
         
      self._explicit = False

      for f, v in self.fields.iteritems():
         self.__dict__[f] = v[0](**v[1])
      for a, v in self.attribs.iteritems():
         self.__dict__[a] = v[0](**v[1])
      
   def adapt(self):
      """ Dummy function being called after the parsing of attributes
          and before the parsing of fields. """
      
      pass
      
   def store(self, value):
      """Dummy function for storing data."""

      self._explicit = True
      pass
            
   def fetch(self):
      """Dummy function to retrieve data."""

      self.check()
      pass
            
   def check(self):
      """Dummy function to check for input errors."""

      if not (self.units in UnitMap[self._dimension]):
         raise TypeError("Unit type " + self.units + " is not compatible with dimension " + self._dimension)
      if not (self._explicit or self._optional):
         raise ValueError("Uninitialized Input value of type " + type(self).__name__)
   
   def write(self, name="", indent=""): 
      """Writes data in xml file format.

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
         NameError: Raised if one of the tags in the xml input file is 
            incorrect.
         AssertionError: Raised if a consistancy check is failed in reading
            the data.
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
      
   def help(self, level=0, stop_level=None, format=None):

      if (not stop_level is None and level > stop_level):
         return ""

      rstr = ""
      if level == 0:
         rstr += r"\documentclass[12pt,fleqn]{report}"
         rstr += "\n\\begin{document}\n"
      
      rstr += self._help + "\n"
      
      if self._dimension != "undefined": 
         rstr += r"{\\ \bf DIMENSION: }" + self._dimension + "\\\\\n"

      if self._default != None: 
         rstr += r"{\\ \bf DEFAULT: }" + str(self._default) + "\\\\\n"

      if hasattr(self, "type"): 
         rstr += r"{\\ \bf DATA TYPE: }" + self.type.__name__ + "\\\\\n"
      
      if len(self.attribs) != 0 and level != stop_level:
         rstr += "\paragraph{Attributes}\n\\begin{itemize}\n"
         for a in self.attribs:
            rstr += "\item {\\bf " + a + "}:\n " + self.__dict__[a].help(level+1, stop_level)
         rstr += "\end{itemize}\n\n"
            
      if len(self.fields) != 0 and level != stop_level:
         rstr += "\paragraph{Fields}\n\\begin{itemize}\n"
         for f in self.fields:
            rstr += "\item {\\bf " + f + "}:\n " + self.__dict__[f].help(level+1, stop_level)
         rstr += "\end{itemize}\n\n"

      if level == 0:
         rstr += r"\end{document}"
       
      return rstr
       
class InputValue(Input):
   """Scalar class for restart handling.

   Has the methods for dealing with simple data tags of the form:
   <tag_name> data </tag_name>, where data is just a value. Takes the data and
   converts it to the required data_type, so that it can be used in the 
   simulation.

   Attributes:
      type: Data type of the data.
      value: Value of data. Also specifies data type if type is None.
      default: Default value if the tag is not specified.
      units: a Units object or None, used for automatic input conversion
   """
 
   def __init__(self,  help=None, dimension=None, units = None, default=None, dtype=None, value=None, options=None):
      """Initialises InputValue.

      Args:
         dtype: An optional data type. Defaults to None.
         value: An optional value for the data. Defaults to None. Also
            specifies the data type if type is None.
         default: Optional default value if the tag is not specified. Defaults
            to None.
         units: A string defining the _kind_ of units assigned to the variable
            e.g. "energy", "time"...
         attribs: A list of the attribs found in the XML file.
      """

      super(InputValue,self).__init__(help, dimension, units, default)

      if not dtype is None:
         self.type = dtype
      elif not self.value is None:
         self.type = type(value)
      else:
         raise TypeError("You must provide either value or dtype")

      if options is not None:
         self._valid = options
         if not default is None and not self._default in self._valid:
            raise ValueError("Default value not in option list " + str(self._valid))
      else:
         self._valid = None
      
      if not value is None:
         self.store(value)
      elif not self._default is None:
         self.store(self._default)

      self.attribs={}
      
   def store(self, value):
      """Converts the data to the appropriate data type and stores it.

      Args:
         value: The raw data to be stored.
      """     
      
      super(InputValue,self).store(value)
      self.value = self.type(value)
      if self._dimension != "undefined":
         self.value /= UnitMap[self._dimension][self.units]
      
   def fetch(self): 
      """Returns the stored data."""

      super(InputValue,self).fetch()
      if self._dimension != "undefined":
         return self.value*UnitMap[self._dimension][self.units]
      else:
         return self.value

   def check(self):
      
      super(InputValue,self).check()
      if not (self._valid is None or self.value in self._valid):         
         raise ValueError(str(self.value) + " is not a valid option ("+str(self._valid)+")")

      
   def write(self, name="", indent=""): 
      """Writes data in xml file format.

      Writes the data in the appropriate format between appropriate tags.

      Args:
         name: An optional string giving the tag name. Defaults to "".
         indent: An optional string giving the string to be added to the start
            of the line, so usually a number of tabs. Defaults to "".

      Returns:
         A string giving the stored value in the appropriate xml format.
      """

      rstr = indent + "<" + name
      if self.units == "":
         rstr += ">"
      else:
         rstr += " units='" + self.units + "'>"
      rstr += write_type(self.type, self.value) + "</" + name + ">\n"
      return rstr
   
   def parse(self, xml=None, text=""):
      """Reads the data for a single value from an xml file.

      Args:
         xml: An xml_node object containing the all the data for the parent
            tag.
         text: The data held between the start and end tags.
      """

      self._explicit = True
      if xml is None:
         self.value = read_type(self.type, text)
      else:   
         self.value = read_type(self.type, xml.fields["_text"])      
         if "units" in xml.attribs:
            self.units = xml.attribs["units"]
      if hasattr(self, "opt_list"):
         self.index = self.opt_list.name2index(self.value)           


ELPERLINE = 5
class InputArray(Input):
   """Array class for restart handling.

   Has the methods for dealing with simple data tags of the form:
   <tag_name shape="(shape)"> data </tag_name>, where data is just an array and
   of the form [data[0], data[1], ... , data[length]]. 
   
   Takes the data and converts it to the required data_type, 
   so that it can be used in the simulation. Also holds the shape of the array,
   so that we can use a simple 1D list of data to specify a multi-dimensional
   array.

   Attributes:
      type: Data type of the data.
      value: Value of data. Also specifies data type if type is None.
      default: Default value if the tag is not specified.
   """
 
   attribs = { "shape" : (InputValue, 
                 {"dtype": tuple, 
                  "help": "The shape of the array",
                  "default": (0,)}) }

   def __init__(self, help=None, dimension=None, units=None, default=None, dtype=None, value=None):
      """Initialises InputArray.

      Args:
         dtype: An optional data type. Defaults to None.
         value: An optional value for the data. Defaults to None. Also
            specifies the data type if type is None.
         default: Optional default value if the tag is not specified. Defaults
            to None.
      """

      super(InputArray,self).__init__(help, dimension, units, default)
      
      if not dtype is None:
         self.type = dtype
      elif not self.value is None:
         self.type = type(value)
      else:
         raise TypeError("You must provide either value or dtype")

      if not value is None:
         self.store(value)
      elif not self._default is None:
         self.store(self._default)
         
   def store(self, value):
      """Converts the data to the appropriate data type and shape and stores it.

      Args:
         value: The raw data to be stored.
      """

      super(InputArray,self).store(value)
      self.value = np.array(value, dtype=self.type).flatten().copy()
      if self._dimension != "undefined":
         self.value /= UnitMap[self._dimension][self.units]
      if self.shape.fetch() == (0,):
         self.shape.store((len(self.value),))
      
   def fetch(self): 
      """Returns the stored data."""

      super(InputArray,self).fetch()

      if self.shape.fetch() == (0,):
         value = np.resize(self.value,0).copy()
      else:
         value = self.value.reshape(self.shape.fetch()).copy()

      if self._dimension != "undefined":
         return value*UnitMap[self._dimension][self.units]
      else:
         return value

   def write(self, name="", indent=""): 
      """Writes data in xml file format.

      Writes the data in the appropriate format between appropriate tags. Note
      that only ELPERLINE values are printed on each line if there are more
      than this in the array. If the values are floats, or another data type
      with a fixed width of data output, then they are aligned in columns.

      Args:
         name: An optional string giving the tag name. Defaults to "".
         indent: An optional string giving the string to be added to the start
            of the line, so usually a number of tabs. Defaults to "".

      Returns:
         A string giving the stored value in the appropriate xml format.
      """

      rstr = indent + "<" + name + " shape='" + write_tuple(self.shape.fetch())
      if self.units == "":
         rstr += "'>"
      else:
         rstr += "' units='" + self.units + "'>"

      if (len(self.value) > ELPERLINE):
         rstr += "\n" + indent + " [ "
      else:
         rstr += " [ "

      for i, v in enumerate(self.value):
         if (len(self.value) > ELPERLINE and i > 0 and i%ELPERLINE == 0):
            rstr += "\n" + indent + "   "
         rstr += write_type(self.type, v) + ", "

      rstr = rstr.rstrip(", ")
      if (len(self.value) > ELPERLINE):
         rstr += " ]\n" + indent
      else:
         rstr += " ] "

      rstr += "</" + name + ">\n"
      return rstr

   def parse(self, xml=None, text=""):
      """Reads the data for a single value from an xml file.

      Args:
         xml: An xml_node object containing the all the data for the parent
            tag.
         text: The data held between the start and end tags.
      """

      self._explicit = True
      if xml is None:
         pass 
      else:   
         self.value = read_array(self.type, xml.fields["_text"])      
         if "units" in xml.attribs:
            self.units = xml.attribs["units"]
         if "shape" in xml.attribs:
            self.shape.store(read_type(tuple, xml.attribs["shape"]))
         else:
            self.shape.store(self.value.shape)



class InputTest(Input):
   
   fields = {
        "timestep" : ( InputValue, { "help"  : "time step of the simulation",
                                     "dtype" : float,
                                     "dimension" : "time",
                                     "default" : 1.0     } ), 
        "label" :    ( InputValue, { "help" : "name of something", 
                                     "dtype" : str,
                                     "default" : "name2",
                                     "options" : ["name1", "name2", "name3"]}),
        "mandatory" : ( InputValue, {"help" : "an integer which must be there", 
                                     "dtype": int } ), 
        "array" : ( InputArray, { "help" : "lets try an array variable", 
                                  "dtype": float,
                                  "dimension": "mass",
                                  "default": np.array([1,2,3,4,5])}),
        "mandarray" : ( InputArray, { "help" : "more array variables", 
                                      "dtype": int})
        }

   default_help = "A test for the new input class"        

class InputTest2(Input):
   fields = { "nest" : (InputTest, {"help" : "This should test nesting"}),
              "other": (InputValue,{"help" : "Another field",
                                    "dtype": bool,
                                    "default": True})}

   default_help = "A test for a nested new input class"
     
     
def main(file_name):
   # TEMPORARY WHILE WE TEST THE NEW INPUT    
   myinput = InputTest()
   myinput2 = InputTest2()


   ifile = open("test.xml","r")
   ifile2 = open("test2.xml","r")
   xmlrestart = xml_parse_file(ifile) # Parses the file.
   xmlrestart2 = xml_parse_file(ifile2)
      
   myinput.parse(xmlrestart.fields["test"]) 
   myinput2.parse(xmlrestart2.fields["test"]) 

   print "timestep: ", myinput.timestep.fetch()
   print "label: ", myinput.label.fetch()
   print "mandatory: ", myinput.mandatory.fetch()
   print "array: ", myinput.array.fetch()
   print "opt_list: ", myinput.label._valid

   print
   print myinput.write("test")
   print
   print
   print "other: ", myinput2.other.fetch()
   print "label: ", myinput2.nest.label.fetch()
   print
   print myinput2.write("test")

#helpfile = open("inputtest/helpfile.tex","w")
#helpfile.write(myinput2.help(stop_level = 2))
#helpfile.write(myinput2.help())

