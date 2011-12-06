import numpy
# if A depends upon B, then A.dep_up-->B and B.dep_dw-->A
class depend_proxy(object):
   """Prototype class for dependency handling"""
   def __init__(self, value=None, name=None, dependants=[], dependencies=[]):
      self._value=value
      self._dependants=[]
      self._name=name
      self._tainted=False
      for item in dependencies:
         item.add_dependant(self)
      for item in dependants:
         self.add_dependant(item)
   
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
   def __init__(self, func, value=None, name=None, dependants=[], dependencies=[]):
      self._func=func
      super(depend_func, self).__init__(value=value,name=name, dependants=dependants, dependencies=dependencies)
         
   def val_update(self, manual=True):
      self._value.set(self._func(), manual=manual)

   def val_set(self, manual=True):
      if (manual):
         print "Cannot set manually the value of the automatically-computed property <",self.name,">"
         raise NameError(self.name)

class synchronizer(object):
   def __init__(self, deps=dict()):
      self._synced=deps
      self._manual=None

class depend_sync(depend_proxy):
   """Proxy which allows to keep different representations of the same 
      quantity synchronized. Only one can be set manually."""
   def __init__(self, func, synchro, value=None, name=None, dependants=[], dependencies=[]):
      self._func=func
      self.synchro=synchro
      self.synchro._synced[name]=self
      self.synchro._manual=name
      self._tainted=False
      super(depend_sync, self).__init__(value=value,name=name, dependants=dependants, dependencies=dependencies)

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
   def __init__(self, deps=None, name=None):
      if deps==None:
         self.deps=depend_proxy()
      else:
         self.deps=deps

class depend_value(depend_base):
   def __init__(self, value, deps=None, name=None):
      super(depend_value,self).__init__(deps, name)
      self.deps._value=self
      self._value=value
      self.deps.taint(taintme=False)

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
   def __new__(cls, input_array, deps=None, name=None, parent=None):
      # Input array is an already formed ndarray instance
      # We first cast to be our class type
      obj = numpy.asarray(input_array).view(cls)
      # add the new attribute to the created instance
      return obj
      
   def __init__(self, input_array, deps=None, name=None, parent=None):
      super(depend_array,self).__init__(deps, name)
      self.parent = parent
      if parent is not None:
         self.deps = parent.deps
      else:
         self.deps._value=self
   
   def __array_finalize__(self, obj): pass

   def __getitem__(self,index):
#      if self.parent is not None:
#         if self.parent.deps.tainted():
#            self.deps.val_update(manual=False)
#      else:
      if self.deps.tainted():  
         self.deps.val_update(manual=False)
         
      if (not numpy.isscalar(index)):      
         #if self.parent is not None:
         #   parent = self.parent
         #else:
         #   parent = self
         a = depend_array(super(depend_array,self).__getitem__(index), parent=self)#parent)#, deps=self.deps   )
         return a
      else:
         return super(depend_array,self).__getitem__(index)   

   def __getslice__(self,i,j):
      return self.__getitem__(slice(i,j))

   def get(self):
      return self[:]
      
   def __get__(self, instance, owner):
      return self.get() 

   def __setitem__(self,index,value,manual=True):
#      if self.parent is not None:
#         self.parent.deps.taint(taintme=False)
#         self.parent.deps.val_set(manual=manual)
#      else:
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
