import numpy as np
import pdb, time

class synchronizer(object):
   def __init__(self, deps=None):
      if deps is None:
         self.synced=dict()
      else:
         self.synced=deps
         
      self.manual=None

# if A depends upon B, then A.dep_up-->B and B.dep_dw-->A
class depend_base(object):
   """Prototype class for dependency handling"""
   def __init__(self, name="", synchro=None, func=None, dependants=None, dependencies=[], tainted=None):
      self._dependants=[]
      # initializes defaults
      if tainted is None:     tainted=np.array([True],bool)
      if dependants is None:  dependants=[]
      self._tainted=tainted
      self._func=func
      self._name=name
      self._synchro=synchro
      if not self._synchro is None and not self._name in self._synchro.synced:
         self._synchro.synced[self._name]=self
         self._synchro.manual=self._name
         
      # don't want to 
      for item in dependencies:  item.add_dependant(self, tainted)
      
      self._dependants=dependants
      if (tainted): self.taint(taintme=tainted)
   
   def add_dependant(self,newdep, tainted=True):
      """Makes newdep dependent on self"""
      self._dependants.append(newdep)
      if tainted: newdep.taint(taintme=True)      

   def add_dependency(self,newdep, tainted=True):
      """Makes self dependent on newdep"""
      newdep._dependants.append(self)      
      if tainted: self.taint(taintme=True)

   def taint(self,taintme=True):
      """Recursively sets tainted flag on dependent objects."""
      #if not self._linkto is None and self._linkto.name == "vir": pdb.set_trace()      
      self._tainted[:] = True     #this is to prevent circular dependencies to hang forever
      for item in self._dependants: 
         if (not item.tainted()):  item.taint()
      if not self._synchro is None:
         for v in self._synchro.synced.values():
            if (not v.tainted()) and (not v is self): v.taint(taintme=True)
         self._tainted[:]=(taintme and (not self._name == self._synchro.manual))         
      else: self._tainted[:] = taintme
      
   def tainted(self):
      """Returns tainted flag"""
      return self._tainted[0]
      
   #TODO put some error checks in the init to make sure that the object is initialized from consistent synchro and func states
   def update_auto(self):
      if not self._synchro is None:
         if (not self._name == self._synchro.manual):
            self.set(self._func[self._synchro.manual](), manual=False)
         else: print "####"+self._name+" shouldn't probably be tainted (synchro)!"
      elif not self._func is None: 
         self.set(self._func(), manual=False)
      else: 
         print "####"+self._name+" shouldn't probably be tainted (value)!"
         pass

   def update_man(self): 
      if not self._synchro is None:     
         self._synchro.manual=self._name
         for v in self._synchro.synced.values():
            v.taint(taintme=True)
         self._tainted[:]=False      
      elif not self._func is None:
         raise NameError("Cannot set manually the value of the automatically-computed property <"+self._name+">")
      else: self.taint(taintme=False)     
            
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
      if (manual) : self.update_man()
      
   def __set__(self, instance, value): 
      self.set(value)   

class depend_array(np.ndarray, depend_base):
   def __new__(cls, value, name="", synchro=None, func=None, dependants=None, dependencies=[], tainted=None, storage=None):
#      print "__new__"
      # Input array is an already formed ndarray instance
      # We first cast to be our class type
      obj = np.asarray(value).view(cls)
      return obj

   def __init__(self, value, name="", synchro=None, func=None, dependants=None, dependencies=[], tainted=None, storage=None):
#      print "__init__"
      super(depend_array,self).__init__(name, synchro, func, dependants, dependencies, tainted)
      self._storage=storage
      if storage is None: self._storage=value  #keeps track of where the original data is pointing, as automatic updates should access there
      else: self._storage=storage
            
#      print "init", name
#      super(depend_array,self).__init__(deps=deps, name=name, tainted=tainted)
#      if self.deps._linkto is None: self.deps._linkto = self
   
   def __array_finalize__(self, obj):  
#      print "__finalize__"   
#      print "finalize", type(self), type(obj),  hasattr(self,"name"),  hasattr(obj,"name")
      # makes sure that --if we really mean to return a deparray-- some basic dep things are provided
      depend_base.__init__(self)  #explicitly initialize in case we got here
      self._storage=self
#      self.deps=depend_proxy()
   
   # whenever possible in compound operations just return a regular ndarray
   __array_priority__=-1.0  
#   def __array_wrap__(self, out_arr, context=None):
#      print "array_wrap", self.name, out_arr.name #, context
#      return super(depend_array,self).__array_wrap__(self, out_arr, context)      
#      return super(depend_array,self).__array_wrap__(self, out_arr, context).view(np.ndarray)

#   def __array_prepare__(self, out_arr, context=None):        
#      print "array_prepare"
#      return super(depend_array,self).__array_prepare__(self, out_arr, context).view(np.ndarray)
      
   def reshape(self, newshape):
      return depend_array(self.base.reshape(newshape), name=self._name, synchro=self._synchro, func=self._func, dependants=self._dependants, dependencies=[], tainted=self._tainted, storage=self._storage)  

   def flatten(self):
      return self.reshape(self.size)

   def __getitem__(self,index):
#      print "getitem", self.name, self.deps.tainted(), index
      
      if self.tainted():  
         self.update_auto()
         self.taint(taintme=False)
              
      if (not np.isscalar(index) or self.ndim > 1 ):
#         return depend_array(super(depend_array,self).__getitem__(index), deps=self.deps, name=self.name, tainted=self.deps._tainted)  
#         return depend_array(self.view(np.ndarray)[index], name=self._name, synchro=self._synchro, func=self._func, dependants=self._dependants, dependencies=[], tainted=self._tainted)  
         return depend_array(self.base[index], name=self._name, synchro=self._synchro, func=self._func, dependants=self._dependants, dependencies=[], tainted=self._tainted, storage=self._storage)  
      else:
         return self.view(np.ndarray)[index]

   def __getslice__(self,i,j):
      return self.__getitem__(slice(i,j,None))

   def get(self):
      return self.__getitem__(slice(None,None,None))
            
   def __get__(self, instance, owner):
      return self.__getitem__(slice(None,None,None))

   def __setitem__(self,index,value,manual=True):      
      #print "setitem", manual, self.name, self.deps._tainted    
      self.taint(taintme=False)      
      #super(depend_array,self).__setitem__(index,value)   # directly write to the base array
      if manual:
         #self.view(np.ndarray)[index]=value
         self.base[index]=value
      else:
         self._storage[index]=value
      if (manual) : self.update_man()
      
   def __setslice__(self,i,j,value):
      return self.__setitem__(slice(i,j),value)

   def set(self, value, manual=True):
      self.__setitem__(slice(None,None),value=value,manual=manual)    

   def __set__(self, instance, value): 
      self.__setitem__(slice(None,None),value=value)

#   def __str__(self):
#      return str(self.get())

def dget(obj,member):
   return obj.__dict__[member]
def dset(obj,member,value):
   obj.__dict__[member]=value
   
def depstrip(deparray):
   return deparray.view(np.ndarray)

def depcopy(objfrom,memberfrom,objto,memberto):
   dfrom=dget(objfrom,memberfrom); dto=dget(objto,memberto);
   dto._dependants=dfrom._dependants
   dto._synchro=dfrom._synchro   
   dto._tainted=dfrom._tainted
   dto._func=dfrom._func
   if hasattr(dfrom,"_storage"): dto._storage=dfrom._storage
   
#def depcopy(objfrom,memberfrom,objto,memberto):
#   dfrom=dget(objfrom,memberfrom); dto=dget(objto,memberto);
#   if not type(dfrom) is type(dto) : raise TypeError("Cannot copy between depend storage of different type")
#   depfrom=depget(objfrom,memberfrom); depto=depget(objto,memberto);
#   if not type(depfrom) is type(depto) : raise TypeError("Cannot copy between depend proxies of different type")   
#   
#   for d in depto._dependants: depfrom.add_dependant(d)
#   if depfrom._linkto is None: depfrom._linkto=depto._linkto
#   dset(objto,memberto,dfrom)


#time.doprint=False
class dobject(object): 
   def __getattribute__(self, name):
      #if time.doprint: print "getattr", name
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
