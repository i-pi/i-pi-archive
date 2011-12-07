import numpy
# if A depends upon B, then A.dep_up-->B and B.dep_dw-->A
class depend_proxy(object):
   """Prototype class for dependency handling"""
   def __init__(self, dependants=[], dependencies=[]):
      self._dependants=[]
      self._tainted=False
      self._value=None
      for item in dependencies:
         item.add_dependant(self)
      for item in dependants:
         self.add_dependant(item)
   
   def link_value(self, value):
      self._value=value
   
   def add_dependant(self,newdep):
      """Makes newdep dependent on self"""
      self._dependants.append(newdep)
      newdep.taint(taintme=True)      

   def add_dependency(self,newdep):
      """Makes self dependent on newdep"""
      newdep._dependants.append(self)      
      self.taint(taintme=True)

   def taint(self,taintme=True):
      """Recursively sets tainted flag on dependent objects."""
      self._tainted = True     #this is to prevent circular dependencies to hang forever
      for item in self._dependants: 
         if (not item.tainted()):
            item.taint()
      self._tainted = taintme
      
   def tainted(self):
      """Returns tainted flag"""
      return self._tainted
      
   def val_update(self, manual=True): pass
   def val_set(self, manual=True): pass
         

class depend_func(depend_proxy):
   """Proxy which defines a function to compute the value of the property. 
      Setting the property manually is forbidden."""
   def __init__(self, func, dependants=[], dependencies=[]):
      self._func=func
      super(depend_func, self).__init__(dependants=dependants, dependencies=dependencies)
         
   def val_update(self, manual=True):
      self._value.set(self._func(), manual=manual)

   def val_set(self, manual=True):
      if (manual):
         print "Cannot set manually the value of the automatically-computed property <",self._value.name,">"
         raise NameError(self._value.name)

class synchronizer(object):
   def __init__(self, deps=dict()):
      self._synced=deps
      self._manual=None

class depend_sync(depend_proxy):
   """Proxy which allows to keep different representations of the same 
      quantity synchronized. Only one can be set manually."""
   def __init__(self, func, synchro, dependants=[], dependencies=[]):
      self._func=func
      self.synchro=synchro
      self._tainted=False
      super(depend_sync, self).__init__(dependants=dependants, dependencies=dependencies)

   def link_value(self, value):
      self._value=value
      self.synchro._synced[value.name]=self
      self.synchro._manual=value.name
      
   def taint(self,taintme=True):
      """Recursively sets tainted flag on dependent objects."""
      super(depend_sync,self).taint(taintme)
      # Also taints object within the sync group, making sure that the one which is manually set is not tainted
      self._tainted=True
      for v in self.synchro._synced.values():
         if (not v.tainted()) and (not v is self) and (not v._value.name == self.synchro._manual): v.taint(taintme=True)
      self._tainted=taintme
               
   def val_update(self, manual=True):
      self._tainted=False
      if (not self._value.name == self.synchro._manual):
         self._value.set(self._func[self.synchro._manual](), manual=manual)

   def val_set(self, manual=True):         
      if (manual):
         for v in self.synchro._synced.values():
            v.taint(taintme=True)
         self.synchro._manual=self._value.name
         self._tainted=False
      
class depend_base(object):
   def __init__(self, deps=None, name=None, tainted=True):
      self.name=name
      if deps==None:
         self.deps=depend_proxy()
      else:
         self.deps=deps
      if (self.deps._value is None): self.deps.link_value(self)
      self.deps.taint(taintme=tainted)  # ALWAYS taint on init      

class depend_value(depend_base):
   def __init__(self, value=None, deps=None, name=None, tainted=True):
      super(depend_value,self).__init__(deps, name, tainted)
      self.deps._value=self
      self._value=value

   def get(self):
      if self.deps.tainted():  
         self.deps.val_update(manual=False)
         self.deps.taint(taintme=False)         

      return self._value
      
   def __get__(self, instance, owner):
      return self.get() 

   def set(self, value, manual=True):
      self._value=value
      self.deps.taint(taintme=False)
      self.deps.val_set(manual=manual)
      
   def __set__(self, instance, value): 
      self.set(value)   

   
class depend_array(numpy.ndarray, depend_base):
   def __new__(cls, value, deps=None, name=None, tainted=True):
      # Input array is an already formed ndarray instance
      # We first cast to be our class type
      obj = numpy.asarray(value).view(cls)
      # add the new attribute to the created instance
      return obj
      
   def __init__(self, value, deps=None, name=None, tainted=True):
      super(depend_array,self).__init__(deps, name, tainted)
      if self.deps._value is None: 
         self.deps._value = self
   
   def __array_finalize__(self, obj): pass

   def __getitem__(self,index):
      if self.deps.tainted():  
         self.deps.val_update(manual=False)         
      if (not numpy.isscalar(index) or self.ndim > 1 ):
         return depend_array(super(depend_array,self).__getitem__(index), deps=self.deps, name=self.name, tainted=self.deps._tainted)
      else:
         return super(depend_array,self).__getitem__(index)   

   def __getslice__(self,i,j):
      return self.__getitem__(slice(i,j))

   def get(self):
      return self[:]
      
   def __get__(self, instance, owner):
      return self.get() 

   def __setitem__(self,index,value,manual=True):
      self.deps.taint(taintme=False)
      self.deps.val_set(manual=manual)
      super(depend_array,self).__setitem__(index,value)   

   def __setslice__(self,i,j,value):
      return self.__setitem__(slice(i,j),value)

   def set(self, value, manual=True):
      self.__setitem__(slice(None,None),value=value,manual=manual)
      
   def __set__(self, instance, value): 
      self.set(value)

def dget(obj,member):
   return obj.__dict__[member]

def depget(obj,member):
   return obj.__dict__[member].deps

class dobject(object): 
   def __getattribute__(self, name):
      value = object.__getattribute__(self, name)
      if hasattr(value, '__get__'):
         value = value.__get__(self, self.__class__)
      return value

   def __setattr__(self, name, value):
      try:
         obj = object.__getattribute__(self, name)
      except AttributeError:
         pass
      else:
         if hasattr(obj, '__set__'):
            return obj.__set__(self, value)
      return object.__setattr__(self, name, value)
