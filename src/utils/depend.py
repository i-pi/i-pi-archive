"""Contains the classes that are used to define the dependency network.

The classes defined in this module overload the standard __get__ and __set__
routines of the numpy ndarray class and standard library object class so that 
they automatically keep track of whether anything they depend on has been
altered, and so only recalculate their value when necessary. 

Basic quantities that depend on nothing else can be manually altered in the 
usual way, all other quantities are updated automatically and cannot be changed 
directly. 

The exceptions to this are synchronized properties, which are in effect 
multiple basic quantities all related to each other, for example the bead and
normal mode representations of the positions and momenta. In this case any of
the representations can be set manually, and all the other representations
must keep in step.

For a more detailed discussion, see the reference manual.

Classes:
   depend_base: Base depend class with the generic methods and attributes.
   depend_value: Depend class for scalar objects.
   depend_array: Depend class for arrays.
   synchronizer: Class that holds the different objects that are related to each
      other and the functions to change between the representations,
      and keeps track of which property has been set manually.
   dobject: An extension of the standard library object that overloads __get__
      and __set__, so that we can use the standard syntax for setting and 
      getting the depend object, i.e. foo = value, not foo.set(value).

Functions:
   dget: Gets the dependencies of a depend object.
   dset: Sets the dependencies of a depend object.
   depstrip: Used on an 
   depcopy:
   deppipe:   
"""

__all__ = ['depend_base', 'depend_value', 'depend_array', 'synchronizer', 
           'dobject', 'dget', 'dset', 'depstrip', 'depcopy', 'deppipe']

import numpy as np

class synchronizer(object):
   def __init__(self, deps=None):
      if deps is None:
         self.synced=dict()
      else:
         self.synced=deps
         
      self.manual=None


#TODO put some error checks in the init to make sure that the object is initialized from consistent synchro and func states
class depend_base(object):
   """Prototype class for dependency handling"""
   def __init__(self, name="", synchro=None, func=None, dependants=None, dependencies=None, tainted=None):
      self._dependants=[]
      if tainted is None:
         tainted=np.array([True],bool)
      if dependants is None:
         dependants=[]
      if dependencies is None:
         dependencies=[]
      self._tainted=tainted
      self._func=func
      self._name=name
      self._synchro=synchro
      if not self._synchro is None and not self._name in self._synchro.synced:
         self._synchro.synced[self._name]=self
         self._synchro.manual=self._name
         
      for item in dependencies:
         item.add_dependant(self, tainted)
      
      self._dependants=dependants
      if (tainted): 
         self.taint(taintme=tainted)
   
   def remove_dependant(self, rmdep):
      irm=-1
      for idep in range(len(self._dependants)): 
         if self._dependants[idep] is rmdep: 
            irm=idep
      if irm >=0: 
         self._dependants.pop(irm)
         
   def add_dependant(self, newdep, tainted=True):
      """Makes newdep dependent on self"""
      self._dependants.append(newdep)
      if tainted: 
         newdep.taint(taintme=True)      

   def add_dependency(self, newdep, tainted=True):
      """Makes self dependent on newdep"""
      newdep._dependants.append(self)      
      if tainted: 
         self.taint(taintme=True)

   def taint(self,taintme=True):
      """Recursively sets tainted flag on dependent objects."""
      self._tainted[:] = True
      for item in self._dependants: 
         if (not item.tainted()):
            item.taint()
      if not self._synchro is None:
         for v in self._synchro.synced.values():
            if (not v.tainted()) and (not v is self):
               v.taint(taintme=True)
         self._tainted[:]=(taintme and (not self._name == self._synchro.manual))         
      else: self._tainted[:] = taintme
      
   def tainted(self):
      """Returns tainted flag"""
      return self._tainted[0]
      
   def update_auto(self):
      if not self._synchro is None:
         if (not self._name == self._synchro.manual):
            self.set(self._func[self._synchro.manual](), manual=False)
         else:
            print "####"+self._name+" probably shouldn't be tainted (synchro)!"
      elif not self._func is None: 
         self.set(self._func(), manual=False)
      else: 
         print "####"+self._name+" probably shouldn't be tainted (value)!"

   def update_man(self): 
      if not self._synchro is None:     
         self._synchro.manual=self._name
         for v in self._synchro.synced.values():
            v.taint(taintme=True)
         self._tainted[:]=False      
      elif not self._func is None:
         raise NameError("Cannot set manually the value of the automatically-computed property <"+self._name+">")
      else:
         self.taint(taintme=False)     
            
   def set(self, value, manual=False): pass           
   def get(self): pass      


class depend_value(depend_base):
   def __init__(self, value=None, name="", synchro=None, func=None, dependants=None, dependencies=[], tainted=None):
      self._value=value
      super(depend_value,self).__init__(name, synchro, func, dependants, dependencies, tainted)

   def get(self):
      if self.tainted():  
         self.update_auto()
         self.taint(taintme=False)
         
      return self._value
      
   def __get__(self, instance, owner):
      return self.get() 

   def set(self, value, manual=True):
      self._value=value
      self.taint(taintme=False)
      if (manual):
         self.update_man()
      
   def __set__(self, instance, value): 
      self.set(value)   


class depend_array(np.ndarray, depend_base):
   def __new__(cls, value, name="", synchro=None, func=None, dependants=None, dependencies=None, tainted=None, storage=None):
      obj = np.asarray(value).view(cls)
      return obj

   def __init__(self, value, name="", synchro=None, func=None, dependants=None, dependencies=None, tainted=None, storage=None):
      super(depend_array,self).__init__(name, synchro, func, dependants, dependencies, tainted)
      if storage is None:
         self._storage=value
      else:
         self._storage=storage  #keeps track of where the original data is pointing, as automatic updates should access there
            
   def __array_finalize__(self, obj):  
      depend_base.__init__(self)  #explicitly initialize in case we got here
      self._storage=self
   
   # whenever possible in compound operations just return a regular ndarray
   __array_priority__=-1.0  
      
   def reshape(self, newshape):
      return depend_array(self.base.reshape(newshape), name=self._name, synchro=self._synchro, func=self._func, dependants=self._dependants, tainted=self._tainted, storage=self._storage)  

   def flatten(self):
      return self.reshape(self.size)
   
   @staticmethod
   def __scalarindex(index, depth=1):
      """Checks if an index points at a scalar value.
      
      Returns a logical stating whether a __get__ instruction based
      on index would return a scalar.      

      Arguments:
         index : the index to be checked
         depth : the rank of the array which is being accessed      
      """
      
      if (np.isscalar(index) and depth <= 1):
         return True
      elif (isinstance(index, tuple) and len(index)==depth):
         for i in index:
            if not np.isscalar(i): return False
         return True
      return False
      
   def __getitem__(self,index):
      if self.tainted():  
         self.update_auto()
         self.taint(taintme=False)
           
      if (self.__scalarindex(index, self.ndim)):
         return depstrip(self)[index]
      else:      
         return depend_array(self.base[index], name=self._name, synchro=self._synchro, func=self._func, dependants=self._dependants, tainted=self._tainted, storage=self._storage)

   def __getslice__(self,i,j):
      return self.__getitem__(slice(i,j,None))

   def get(self):
      return self.__getitem__(slice(None,None,None))
            
   def __get__(self, instance, owner):
      return self.__getitem__(slice(None,None,None))

   def __setitem__(self,index,value,manual=True):      
      self.taint(taintme=False)      
      if manual:
         self.base[index]=value
         self.update_man()
      else:
         self._storage[index]=value
      
   def __setslice__(self,i,j,value):
      return self.__setitem__(slice(i,j),value)

   def set(self, value, manual=True):
      self.__setitem__(slice(None,None),value=value,manual=manual)    

   def __set__(self, instance, value): 
      self.__setitem__(slice(None,None),value=value)

def dget(obj,member):
   return obj.__dict__[member]
   
def dset(obj,member,value,name=None):
   obj.__dict__[member]=value
   if not name is None: 
      obj.__dict__[member]._name=name
   
def depstrip(deparray):
   return deparray.view(np.ndarray)

def deppipe(objfrom,memberfrom,objto,memberto):
   dfrom=dget(objfrom,memberfrom)
   dto=dget(objto,memberto) 
   dto._func=lambda : dfrom.get()
   dto.add_dependency(dfrom)
   
def depcopy(objfrom,memberfrom,objto,memberto):
   dfrom=dget(objfrom,memberfrom)
   dto=dget(objto,memberto)
   dto._dependants=dfrom._dependants
   dto._synchro=dfrom._synchro   
   dto._tainted=dfrom._tainted
   dto._func=dfrom._func
   if hasattr(dfrom,"_storage"):
      dto._storage=dfrom._storage
   

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
