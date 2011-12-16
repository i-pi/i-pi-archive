import numpy as np
import types

import pdb

class Restart(object):   
   fields={}
   attribs={}
   def __init__(self):
      print "generating from", self.fields
      for f, v in self.fields.iteritems():    self.__dict__[f]=v[0](*v[1])
      for a, v in self.attribs.iteritems():   self.__dict__[a]=v[0](*v[1])     
      
   def store(self, value): pass
   
   def fetch(self): pass
      
   def check(self): return True
   
   def write(self, name="", indent=""): 
      rstr=indent+"<"+name+" ATTRIBS GO HERE>\n"
      for f in self.fields:
         rstr+=self.__dict__[f].write(f,"   "+indent)
      rstr+="</"+name+">\n"
      return rstr
   
   def parse(self, text): pass
   # <name attr1=val1 attr2=val2>ENCLOSED TEXT</name>
#      parsed=xml_parse(text) ==> { "name" : ( {"attr1":"val1", "attr2"="val2"} , "ENCLOSED TEXT" ) }
#      for p, vp in parsed.iteritems():
#         for a, va in vp[0]:
#            if not a in self.__dict__[p].
#            self.__dict__[p].__dict__[a].parse(va)
#         self.__dict__[p].parse(vp[1])
         
      
      
   
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
   
         
class RestartArray(Restart):
   attribs={ "size" : (RestartValue,(int,0)),  "shape" : (RestartValue,(tuple, ())) }
   def __init__(self, dtype=None, value=None, default=None):
      super(RestartArray,self).__init__()
      
      if not dtype is None:   self.type=dtype
      elif not self.value is None: self.type=value.dtype
      else: raise TypeError("You must provide either value or dtype")

      self.default=default
      if not value is None:  self.store(value)
      elif not default is None: self.store(default)
   
   def store(self, value):
      self.size.store(value.size)
      self.shape.store(value.shape)
      self.value=np.array(value, dtype=self.type).flatten().copy()
      
   def fetch(self): 
      return self.value.reshape(self.shape.fetch()).copy()

   def write(self, name="", indent=""): 
      rstr=indent+"<"+name+" size="+str(self.size.fetch()) +" shape="+str(self.shape.fetch())+">"+str(self.value)+"</"+name+">\n "
      return rstr
