import numpy
#class dobject(object): 
#    def __getattribute__(self, name):
#        value = object.__getattribute__(self, name)
#        if hasattr(value, '__get__'):
#            value = value.__get__(self, self.__class__)
#        return value

#    def __setattr__(self, name, value):
#        try:
#            obj = object.__getattribute__(self, name)
#        except AttributeError:
#            pass
#        else:
#            if hasattr(obj, '__set__'):
#                return obj.__set__(self, value)
#        return object.__setattr__(self, name, value)

## if A depends upon B, then A.dep_up-->B and B.dep_dw-->A
#class depend_proto(object):
#   def __init__(self,name=None,deplist=[]):
#      print "init proto"
#      self._dep_up=[]
#      self._dep_dw=[]
#      self._name=name
#      self._tainted=True
#      for item in deplist:
#         item.add_dependant(self)
#   
#   def add_dependant(self,newdep):
#      """Makes newdep dependent on self"""
#      self._dep_dw.append(newdep)
#      newdep._dep_up.append(self)      
#      newdep.taint(taintme=True)      

#   def add_dependency(self,newdep):
#      """Makes self dependent on newdep"""
#      self._dep_up.append(newdep)
#      newdep._dep_dw.append(self)      
#      newdep.taint(taintme=True)

#   def taint(self,taintme=True):
#      """Recursively sets tainted flag on dependent objects."""
#      self._tainted = True     #this is to prevent circular dependencies to hang forever
#      for item in self._dep_dw: 
#         if (not item._tainted):
#            item.taint()
#      self._tainted = taintme
#      
#   def tainted(self):
#      """Returns tainted flag"""
#      return self._tainted

#class depend_value(depend_proto):
#   def __init__(self,value=None,name=None,deplist=[]):
#      super(depend_value,self).__init__(name=name,deplist=deplist)
#      self._value=value
#      self.taint(taintme=False)
#      
#   def __get__(self, instance, owner):
#      return self.get() 

#   def __set__(self, instance, value): 
#      self.set(value)
#   def get(self):
#      return self._value
#      
#   def set(self, value):
#      self._value=value 
#      self.taint(taintme=False)
#      
#class depend_calc(depend_proto):
#   def __init__(self,func,name=None,deplist=[]):
#      super(depend_calc,self).__init__(name=name,deplist=deplist)
#      self._func=func
#      self._value=None
#      self.taint(taintme=True)
#      
#   def __get__(self, instance, owner):
#      return self.get()
#      
#   def get(self):
#      if (self._tainted):
#         self._value=self._func()
#         self.taint(taintme=False)
#      return self._value

#class depend_sync(depend_proto):       
#   def __init__(self,synclist=[],name=None,deplist=[]):
#      super(depend_sync,self).__init__(name=name,deplist=deplist)
#      super(depend_sync,self).__setattr__("_syncs",  dict())
#      print self._syncs
#      for s in synclist: 
#         self._syncs[s]=depend_value(name=s)
#      self.taint(taintme=False)

#   def __getattribute__(self, name):
#      if not "_syncs" in super(depend_sync,self).__getattribute__("__dict__"): 
#         return super(depend_sync,self).__getattribute__(name)
#      
#      gsync=depend_proto.__getattribute__(self,"_syncs")
#      if name in gsync:
#         print "getting ",name
#         if gsync[name].tainted(): 
#            print "must compute ",name
#            pass
#            # do something to get from the manually set item
#         return gsync[name].get()
#      else:
#         return super(depend_sync,self).__getattribute__(name)

#   def __setattr__(self, name, value):
#      if not "_syncs" in super(depend_sync,self).__getattribute__("__dict__"): 
#         return super(depend_sync,self).__setattr__(name, value)
#         
#      print "setting ", name
#      if name in self._syncs:
#         self._manual=name
#         for s in self._syncs.values(): 
#            s.taint(taintme=True)
#         print "setting ",name
#         self._syncs[name].set(value)
#         self.taint(taintme=False);
#         return self._syncs[name]
#      else:
#         return super(depend_sync,self).__setattr__(name, value)



# if A depends upon B, then A.dep_up-->B and B.dep_dw-->A
class depend_proxy(object):
   def __init__(self, value=None, name=None, dependants=[], dependencies=[]):
      print "init proxy"
      self._value=value
      self._dep_up=[]
      self._dep_dw=[]
      self._name=name
      self._tainted=False
      for item in dependencies:
         item.add_dependant(self)
      for item in dependants:
         self.add_dependant(item)
   
   def add_dependant(self,newdep):
      """Makes newdep dependent on self"""
      self._dep_dw.append(newdep)
      newdep._dep_up.append(self)      
      newdep.taint(taintme=True)      

   def add_dependency(self,newdep):
      """Makes self dependent on newdep"""
      self._dep_up.append(newdep)
      newdep._dep_dw.append(self)      
      self.taint(taintme=True)

   def taint(self,taintme=True):
      """Recursively sets tainted flag on dependent objects."""
      self._tainted = True     #this is to prevent circular dependencies to hang forever
      for item in self._dep_dw: 
         if (not item._tainted):
            item.taint()
      self._tainted = taintme
      
   def tainted(self):
      """Returns tainted flag"""
      return self._tainted
      
   def val_update(self): pass
      

class depend_func(depend_proxy):
   def __init__(self, func, value=None, name=None, dependants=[], dependencies=[]):
      self._func=func
      super(depend_func, self).__init__(value=value,name=name, dependants=dependants, dependencies=dependencies)
         
   def val_update(self):
      self._value.set(self._func())

#class depend_sync(depend_proxy):
#   def __init__(self, funcs, values, name=None, dependants=[], dependencies=[]):
#      self._func=funcs
#      super(depend_func, self).__init__(value=values,name=name, dependants=dependants, dependencies=dependencies)
#         
#   def val_update(self):
#      self._value.set(self._func())


class depend_base(object):
   def __init__(self, deps=None, name=None):
      self.name=name
      if deps==None:
         self.deps=depend_proxy()
      else:
         self.deps=deps

class depend_value(depend_base):
   def __init__(self, value, deps=None, name=None):
      print "init value"
      super(depend_value,self).__init__(deps, name)
      self.deps._value=self
      self._value=value
      self.deps.taint(taintme=False)

   def get(self):
      print "getting value"
      if self.deps.tainted():  
         self.deps.val_update()
         self.deps.taint(taintme=False)         
      
      return self._value
      
   def __get__(self, instance, owner):
      return self.get() 

   def set(self, value):
      self._value=value
      self.deps.taint(taintme=False)
      
   def __set__(self, instance, value): 
      self.set(value)   

   
class depend_array(numpy.ndarray, depend_base):
   def __new__(cls, input_array, deps=None, name=None):
      print "new array"
      # Input array is an already formed ndarray instance
      # We first cast to be our class type
      obj = numpy.asarray(input_array).view(cls)
      # add the new attribute to the created instance
      return obj
      
   def __init__(self, input_array, deps=None, name=None):
      print "init array"
      super(depend_array,self).__init__(deps, name)
      self.deps._value=self
   
   def __array_finalize__(self, obj): pass

   def __getitem__(self,index):
      print "getting ", index
      if self.deps.tainted():  
         self.deps.val_update()
         
      if (not numpy.isscalar(index)):      
         return depend_array(super(depend_array,self).__getitem__(index), deps=self.deps   )
      else:
         return super(depend_array,self).__getitem__(index)   

   def __getslice__(self,i,j):
      return self.__getitem__(slice(i,j))

   def get(self):
      return self[:]
      
   def __get__(self, instance, owner):
      return self.get() 

   def __setitem__(self,index,value):
      self.deps.taint(taintme=False)
      super(depend_array,self).__setitem__(index,value)   

   def __setslice__(self,i,j,value):
      return self.__setitem__(slice(i,j),value)

   def set(self, value):
      self[:]=value
      
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
