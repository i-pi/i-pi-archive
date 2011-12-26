from xml.sax import parseString
from xml.sax.handler import ContentHandler 
import numpy as np
import string

class xml_node(object):
   def __init__(self, attribs={}, name="", fields={}):
      self.attribs=attribs
      self.name=name         
      self.fields=fields

class xml_handler(ContentHandler):
   def __init__(self):
      self.root=xml_node(name="root", fields={})
      self.open=[self.root]
      self.level=0
      self.buffer=[""]
      
   def startElement(self, name, attrs): 
      newnode=xml_node(attribs=dict((k,attrs[k]) for k in attrs.keys()), name=name, fields={})
      self.open.append(newnode)
      self.open[self.level].fields[name]=newnode
      self.buffer.append("")
      self.level+=1      

   def characters(self, data):
      self.buffer[self.level]+=data
      
   def endElement(self, name):
      self.open[self.level].fields["_text"]=self.buffer[self.level]
      self.buffer.pop(self.level)
      self.open.pop(self.level)
      self.level-=1

def parse_xml(buf):
   myhandle=xml_handler()
   parseString(buf, myhandle)
   return myhandle.root

def read_type(type, data):
   if not type in readtype_funcs: raise TypeError("Conversion not available for given type")
   return type(readtype_funcs[type](data))

def read_float(data):
   """Takes a formatted line with a float and interprets it"""
   return float(data.replace("D","E"))

def read_int(data):
   """Takes a formatted line with an integer and interprets it"""
   
   return int(data)

def read_bool(data):
   """Takes a formatted line with an integer and outputs a boolean"""
   
   if   data.upper() == "TRUE":  return True
   elif data.upper() == "FALSE": return False
   else: raise ValueError(data+" does not represent a bool value")

def read_list(data, delims="[]", split=",", strip=" \n\t'"):
   """Takes a line with an array of the form: 
      [array[0], array[1], array[2],...], and interprets it"""

   try:
      begin = data.index(delims[0])
      end = data.index(delims[1])
   except ValueError: raise ValueError("Error in list syntax: could not locate delimiters")

   rlist= data[begin+1:end].split(split)
   for i in range(len(rlist)):
      rlist[i]=rlist[i].strip(strip)

   return rlist

def read_array(dtype, data):
   """ Reads a np.ndarray of type dtype """
   rlist=read_list(data)
   for i in range(len(rlist)): rlist[i]=read_type(dtype,rlist[i])
   
   return np.array(rlist, dtype)

def read_tuple(data):
   """Takes a line with a tuple of the form: 
      (tuple[0], tuple[1], tuple[2],...), and interprets it"""

   rlist=read_list(data, delims="()")
   return tuple([int(i) for i in rlist])

def read_dict(data):
   """Takes a line with an map of the form:
      {kw[0]: value[0], kw[1]: value[1], ...}, and interprets it"""

   rlist=read_list(data,delims="{}")
   def mystrip(data): return data.strip(" \n\t'")
   rdict={}
   for s in rlist:
      rtuple=map(mystrip,s.split(":"))      
      if not len(rtuple)==2: raise ValueError("Format for a key:value format is wrong for item"+s)
      rdict[rtuple[0]]=rtuple[1]
      
   return rdict   
      
readtype_funcs = {np.ndarray: read_array, dict: read_dict, float: read_float, int: read_int, bool: read_bool, str: string.strip, tuple: read_tuple}

def write_type(type, data):
   if not type in writetype_funcs: raise TypeError("Conversion not available for given type")
   return writetype_funcs[type](data)

def write_list(data, delims="[]"):
   """Takes a line with an array of the form: 
      [array[0], array[1], array[2],...], and interprets it"""

   rstr=delims[0];
   
   for v in data:  rstr+=str(v)+", "
   
   rstr=rstr.rstrip(", ")
   rstr+=delims[1]
   return rstr

def write_tuple(data): return write_list(data, delims="()")

def write_float(data): return "%16.8e" % (data)

def write_bool(data):  return "%5.5s" % (str(data))

def write_dict(data, delims="{}"):
   
   rstr=delims[0]
   for v in data:
      rstr += str(v)+": "+str(data[v])+", "
   rstr = rstr.strip(", ")
   rstr+=delims[1]
   return rstr

writetype_funcs = {float: write_float, dict: write_dict, int: str, bool: write_bool, str: string.strip, tuple: write_tuple, np.uint : str}
