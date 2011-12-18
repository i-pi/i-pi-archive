import numpy as np
import types

import pdb


from  io.io_xml import *

class Restart(object):   
   fields={}
   attribs={}
   def __init__(self):
      for f, v in self.fields.iteritems():    self.__dict__[f]=v[0](*v[1])
      for a, v in self.attribs.iteritems():   self.__dict__[a]=v[0](*v[1])     
      
   def store(self, value): pass
   
   def fetch(self): pass
      
   def check(self): return True
   
   def write(self, name="", indent=""): 
      rstr=indent+"<"+name;
      for a in self.attribs:      
         rstr+=" "+a+"='"+str(self.__dict__[a].fetch())+"'"
      rstr+=">\n"
      for f in self.fields:
         rstr+=self.__dict__[f].write(f,"   "+indent)
      rstr+=indent+"</"+name+">\n"
      return rstr
   
   def parse(self, xml=None, text=""): 
      if not xml is None:
         for a, v in xml.attribs.iteritems() :
            if a in self.attribs: 
               self.__dict__[a].parse(text=v)
            else : pass #raise an exception
         for f, v in xml.fields.iteritems() :
            if f in self.fields:
               self.__dict__[f].parse(xml=v)
            else : pass #raise an exception            
   
class RestartValue(Restart):
   def __init__(self, dtype=None, value=None, default=None):      
      super(RestartValue,self).__init__()
      if not dtype is None:   self.type=dtype
      elif not self.value is None: self.type=type(value)
      else: raise TypeError("You must provide either value or dtype")

      self.default=default
      if not value is None:  self.store(value)
      elif not default is None: self.store(default)
      
   def store(self, value):
      self.value=self.type(value)
      
   def fetch(self): 
      return self.value
      
   def write(self, name="", indent=""): 
      return indent+"<"+name+">"+str(self.value)+"</"+name+">\n"
   
   def parse(self, xml=None, text=""):
      if xml is None:
         self.value = read_type(self.type, text)
      else:   
         self.value = read_type(self.type, xml.fields["_text"])
         
ELPERLINE=5
class RestartArray(Restart):
   attribs={ "shape" : (RestartValue,(tuple, ())) }
   def __init__(self, dtype=None, value=None, default=None):
      super(RestartArray,self).__init__()
      
      if not dtype is None:   self.type=dtype
      elif not self.value is None: self.type=value.dtype
      else: raise TypeError("You must provide either value or dtype")

      self.default=default
      if not value is None:  self.store(value)
      elif not default is None: self.store(default)
   
   def store(self, value):
      self.shape.store(value.shape)
      self.value=np.array(value, dtype=self.type).flatten().copy()
      
   def fetch(self): 
      return self.value.reshape(self.shape.fetch()).copy()

   def write(self, name="", indent=""): 
      rstr=indent+"<"+name+" shape='"+write_tuple(self.shape.fetch())+"'> ";
      if (len(self.value)>ELPERLINE): rstr+="\n"+indent+" [ "
      else: rstr+=" [ "
      for i,v in enumerate(self.value):          
         if (len(self.value)>ELPERLINE and i>0 and i%ELPERLINE==0): rstr+="\n"+indent
         rstr+=str(v)+", "         
      rstr=rstr.rstrip(", ")
      if (len(self.value)>ELPERLINE): rstr+=" ]\n"+indent
      else: rstr+=" ] "
      rstr+="</"+name+">\n"
      return rstr
      
   def parse(self, xml=None, text=""):
      if xml is None: pass #should print an error!
      else:   
         self.shape.store(read_type(tuple,xml.attribs["shape"]))
         self.value=read_array(self.type,xml.fields["_text"]) 
      

      
