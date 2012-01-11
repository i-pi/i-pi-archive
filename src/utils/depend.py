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
      other and keeps track of which property has been set manually.
   dobject: An extension of the standard library object that overloads __get__
      and __set__, so that we can use the standard syntax for setting and 
      getting the depend object, i.e. foo = value, not foo.set(value).

Functions:
   dget: Gets the dependencies of a depend object.
   dset: Sets the dependencies of a depend object.
   depstrip: Used on a depend_array object, to access its value without
      needing the depend machinery, and so much more quickly. Must not be used
      if the value of the array is to be changed.
   depcopy: Copies the dependencies from one object to another
   deppipe: Used to make two objects be synchronized to the same value.
"""

__all__ = ['depend_base', 'depend_value', 'depend_array', 'synchronizer', 
           'dobject', 'dget', 'dset', 'depstrip', 'depcopy', 'deppipe']

import numpy as np

class synchronizer(object):
   """Class to implement synched objects.

   Holds the objects used to keep two or more objects in step with each other. 
   This is shared between all the synched objects.

   Attributes:
      synched: A dictionary containing all the synched objects, of the form
         {"name": depend object}.
      manual: A string containing the name of the object being manually changed.
   """

   def __init__(self, deps=None):
      """Initialises synchronizer.

      Args:
         deps: Optional dictionary giving the synched objects of the form
            {"name": depend object}.
      """

      if deps is None:
         self.synced=dict()
      else:
         self.synced=deps
         
      self.manual=None


#TODO put some error checks in the init to make sure that the object is initialized from consistent synchro and func states
class depend_base(object):
   """Base class for dependency handling.

   Builds the majority of the machinery required for the different depend
   objects. Contains functions to add and remove dependencies, the tainting 
   mechanism by which information about which objects have been updated is
   passed around the dependency network, and the manual and automatic update
   functions to check that depend objects with functions are not manually
   updated and that synchronized objects are kept in step with the one manually
   changed.

   Attributes:
      _tainted: An array containing one boolean, which is True if one of the
         dependencies has been changed since the last time the value was 
         cached.
      _func: A function name giving the method of calculating the value, 
         if required. None otherwise.
      _name: The name of the depend base object.
      _synchro: A synchronizer object to deal with synched objects, if
         required. None otherwise.
      _dependants: A list containing all objects dependent on the self.
   """ 

   def __init__(self, name, synchro=None, func=None, dependants=None, dependencies=None, tainted=None):
      """Initialises depend base.

      An unusual initialisation routine, as it has to be able to deal with the
      depend array mechanism for returning slices as new depend arrays.

      This is the reason for the penultimate if statement; it automatically
      taints objects created from scratch but does nothing to slices which are
      not tainted.

      Also, the last if statement makes sure that if a synchronized property is
      sliced, this initialization routine does not automatically set it to the
      manually updated property.

      Args:
         name: A string giving the name of self.
         tainted: An optional array containing one boolean which is True if one
         of the dependencies has been changed.
         func: An optional argument that can be specified either by a function 
            name, or for synchronized values a dictionary of the form 
            {"name": function name}; where "name" is one of the other
            synched objects and function name is the name of a function to
            get the object "name" from self.
         synchro: An optional synchronizer object.
         dependants: An optional list containing objects that depend on self.
         dependencies: An optional list containing objects that self 
            depends upon.
   """

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
         
      for item in dependencies:
         item.add_dependant(self, tainted)
      
      self._dependants=dependants
      if (tainted): 
         self.taint(taintme=tainted)

      if not self._synchro is None and not self._name in self._synchro.synced:
         self._synchro.synced[self._name]=self
         self._synchro.manual=self._name
   
   def remove_dependant(self, rmdep):
      """Removes a dependency.

      Args:
         rmdep: The depend object to be removed from the dependency list.
      """

      irm=-1
      for idep in range(len(self._dependants)): 
         if self._dependants[idep] is rmdep: 
            irm=idep
      if irm >=0: 
         self._dependants.pop(irm)
         
   def add_dependant(self, newdep, tainted=True):
      """Adds a dependant property.

      Args:
         newdep: The depend object to be added to the dependency list.
         tainted: A boolean that decides whether newdep should be tainted.
            True by default.
      """

      self._dependants.append(newdep)
      if tainted: 
         newdep.taint(taintme=True)      

   def add_dependency(self, newdep, tainted=True):
      """Adds a dependency.

      Args:
         newdep: The depend object self now depends upon.
         tainted: A boolean that decides whether self should 
            be tainted. True by default.
      """

      newdep._dependants.append(self)      
      if tainted: 
         self.taint(taintme=True)

   def taint(self,taintme=True):
      """Recursively sets tainted flag on dependent objects.

      The main function dealing with the dependencies. Taints all objects 
      further down the dependency tree until either all objects have been
      tainted, or it reaches only objects that have already been tainted. Note
      that in the case of a dependency loop the initial setting of _tainted to
      True prevents an infinite loop occuring.

      Also, in the case of a synchro object, the manually set quantity is not
      tainted, as it is assumed that synchro objects only depend on each other.

      Args:
         taintme: A boolean giving whether self should be tainted at the end.
            True by default.
      """

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
      """Returns tainted flag."""

      return self._tainted[0]
      
   def update_auto(self):
      """Automatic update routine.

      Updates the value when get has been called and self has been tainted.
      """

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
      """Manual update routine.

      Updates the value when the value has been manually set. Also raises an
      exception if a calculated quantity has been manually set. Also starts the
      tainting routine.

      Raises:
         NameError: If a calculated quantity has been manually set.
      """

      if not self._synchro is None:     
         self._synchro.manual=self._name
         for v in self._synchro.synced.values():
            v.taint(taintme=True)
         self._tainted[:]=False      
      elif not self._func is None:
         raise NameError("Cannot set manually the value of the automatically-computed property <"+self._name+">")
      else:
         self.taint(taintme=False)     
            
   def set(self, value, manual=False):
      """Dummy setting routine."""

      pass           

   def get(self):      
      """Dummy getting routine."""

      pass


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

      Used so that looking up one item in an array returns a scalar, whereas
      looking up a slice of the array returns a new array with the same
      dependencies as the original, so that changing the slice also changes
      the global array.

      Arguments:
         index : the index to be checked
         depth : the rank of the array which is being accessed      
      
      Returns:
         A logical stating whether a __get__ instruction based
         on index would return a scalar.      
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
