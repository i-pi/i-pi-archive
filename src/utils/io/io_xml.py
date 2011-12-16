import numpy as np
import sys, string
sys.path.append('./')
import utils.restart
import xml.sax.handler, xml.sax, pprint

def read_dict(data):
   """Takes a line with an map of the form:
      {kw[0]: value[0], kw[1]: value[1], ...}, and interprets it"""

   try:
      begin = data.index("{")
      end = data.index("}")
   except:
      print "Error in map syntax"
      exit()
   
   elements = data.count(",") + 1
   if data.count(":") != elements:
      print "Error in map syntax"
      exit()

   length = len(data)
   comma_list = [i for i in range(length) if data[i] == ","]
   colon_list = [i for i in range(length) if data[i] == ":"]

   try:
      if elements == 1:
         output = {}
         kw = data[begin+1:colon_list[0]]
         kw = string.strip(kw)
         value = data[colon_list[0]+1:end]
         output[kw] = string.strip(value)
      else:
         output = {}
         kw = data[begin+1:colon_list[0]]
         kw = string.strip(kw)
         value = data[colon_list[0]+1:comma_list[0]]
         output[kw] = string.strip(value)

         kw = data[comma_list[elements-2]+1:colon_list[elements-1]]
         kw = string.strip(kw)
         value = data[colon_list[elements-1]+1:end]
         output[kw] = string.strip(value)

         for i in range(1, elements-1):
            kw = data[comma_list[i-1]+1:colon_list[i]]
            kw = string.strip(kw)
            value = data[colon_list[i]+1:comma_list[i]]
            output[kw] = string.strip(value)
      return output
   except ValueError:
      print "Error in writing to string in map"
      exit()

def read_array(data, dtype):
   """Takes a line with an array of the form: 
      [array[0], array[1], array[2],...], and interprets it"""

   try:
      begin = data.index("[")
      end = data.index("]")
   except ValueError:
      print "Error in array syntax"
      exit()

   elements = data.count(",") + 1
   length = len(data)
   comma_list = [i for i in range(length) if data[i] == ","]
   for i in range(length):
      if data[i] == "D":
         data = data[0:i] + "E" + data[i+1:length]
  
   try:
      if elements == 1:
         output = np.zeros(elements, dtype)
         output[0] = dtype(data[begin+1:end])
      else:
         output = np.zeros(elements, dtype)
         output[0] = dtype(data[begin+1:comma_list[0]])
         output[elements-1] = dtype(data[comma_list[elements-2]+1:end])
         for i in range(1,elements-1):
            output[i] = dtype(data[comma_list[i-1]+1:comma_list[i]])
      return output
   except ValueError:
      print "Tried to write a value to array not of form ", dtype
      exit()

def read_tuple(data):
   """Takes a line with a tuple of the form: 
      (tuple[0], tuple[1], tuple[2],...), and interprets it"""

   try:
      begin = data.index("(")
      end = data.index(")")
   except ValueError:
      print "Error in tuple syntax"
      exit()

   elements = data.count(",") + 1
   length = len(data)
   comma_list = [i for i in range(length) if data[i] == ","]
   for i in range(length):
      if data[i] == "D":
         data = data[0:i] + "E" + data[i+1:length]
  
   try:
      if elements == 1:
         output = np.zeros(elements)
         output[0] = float(data[begin+1:end])
      else:
         output = np.zeros(elements)
         output[0] = float(data[begin+1:comma_list[0]])
         output[elements-1] = float(data[comma_list[elements-2]+1:end])
         for i in range(1,elements-1):
            output[i] = float(data[comma_list[i-1]+1:comma_list[i]])
      return output
   except ValueError:
      print "Tried to write NaN to tuple"
      exit()

def read_float(data):
   """Takes a formatted line with a float and interprets it"""

   output = 0.0
   length = len(data)
   for i in range(length):
      if data[i] == "D":
         data = data[0:i] + "E" + data[i+1:length]
   try:
      output = float(data)
      return output
   except ValueError:
      print "Tried to write NaN to float"
      exit()

def read_int(data):
   """Takes a formatted line with an integer and interprets it"""

   output = 0
   try:
      output = int(data)
      return output
   except ValueError:
      print "Tried to write NaN to int"
      exit()

def read_bool(data):
   """Takes a formatted line with an integer and outputs a boolean"""

   output = False
   try:
      output = bool(int(data))
      return output
   except ValueError:
      print "Tried to write NaN to int in bool"
      exit()

class xml_object:
   def __init__(self, name, default = None, container = None, parent = None, child_list = [], func = None, args = None):
      self.inside = False
      self.name = name
      self.value = None
      self.attributes = None 
      self.default = default
      self.container = container
      self.parent = parent
      self.child_list = child_list
      self.func = func
      self.args = args
      if self.container is None and self.func is None:
         print "Error in xml_object ", self.name, " initialisation, no way of computing value"
      if self.container is not None and self.func is not None:
         print "Error in xml_object ", self.name, " initialisation, container made for primitive value"

   def start_tag(self, attributes):
      if self.parent:
         if not self.parent.inside:
            print "Tag ", self.name, " is not within ", selfparent.name
            exit()
      self.attributes = attributes
      self.inside = True

   def fill(self, data):
      if self.container is None:
         if self.args:
            self.value = self.func(data, self.args)
         else:
            self.value = self.func(data)

   def end_tag(self):
      if self.container is not None:
         for child in self.child_list:
            if child.value is None:
               if child.default is not None:
                  self.container.__dict__[child.name] = child.default
                  print "Using default value for ", child.name
               else:
                  print "Value for required field ", child.name, " not specified"
                  exit()
            else:
               self.container.__dict__[child.name] = child.value
               
         self.value = self.container
      self.inside = False

class empty_container:
   pass

class Init_read(xml.sax.handler.ContentHandler):
   
   def __init__(self, simulation_template):
      self.func_dict = {np.ndarray: read_array, dict: read_dict, float: read_float, int: read_int, bool: read_bool, str: string.strip, tuple: read_tuple}

      self.object_dict = {}
      #self.object_dict["simulation"] = xml_object(name = "simulation", container = empty_container())
      self.object_dict["force"] = xml_object(name = "force", container = empty_container())

      self.init_template(simulation_template, self.object_dict["force"])

   def init_template(self, fields, parent):
      parent.child_list = []
      for name in fields:
         field = fields[name]
         if field[0] is utils.restart.RestartValue:
            for val_type in self.func_dict:
               if field[1][0] is val_type:
                  self.object_dict[name] = xml_object(name = name, default = field[1][1], func = self.func_dict[val_type], parent = parent)
                  parent.child_list.append(self.object_dict[name])
         elif field[0] == utils.restart.RestartArray:
            self.object_dict[name] = xml_object(name = name, default = field[1][1], func = read_array, args = field[1][0], parent = parent)
            parent.child_list.append(self.object_dict[name])
         else:
            try:
               self.object_dict[name] = xml_object(name = name, container = empty_container(), default = field[1][1], parent = parent)
               parent.child_list.append(self.object_dict[name]) 
               self.init_template(field[1][0].fields, parent = self.object_dict[name])
            except AttributeError:
               print "Template has unrecognized type ", field[1][0]
               exit()

   def startElement(self, name, attributes):
      try:
         self.object_dict[name].start_tag(attributes) 
      except KeyError:
         print "Unrecognized starting tag ", name
         exit()
      
   def characters(self, data):
      for name in self.object_dict:
         if self.object_dict[name].inside is True:
            self.object_dict[name].fill(data)

   def endElement(self, name):
      try:
         self.object_dict[name].end_tag()
      except KeyError:
         print "Unrecognized ending tag ", name

def init_from_xml(input_file, template):
   parser = xml.sax.make_parser()
   handler = Init_read(template)
   parser.setContentHandler(handler)
   parser.parse(input_file)

   return handler
