import numpy

class depend(object):
   """ A descriptor class for quantities that may depend
   from other quantities or be dependencies for other quantities.
   Usage is:
   x=depend(func=<function>,deplist=<[list of dependencies]>,
            value=<initial value>,name=<name>)

   Get and set must be called explicitly, and it is possible to 
   taint directly the object.
   Whenever one of the dependencies get tainted, then <function>
   will be called on __get__, and a new value generated. Otherwise,
   the cached value is returned. 
   Synced objects are similarly cached, and marked as tainted when
   one explicitly sets their value, but it is assumed that call to 
   function will leave the synced objects in a consistent state. 
   Also, it is ensured that synchronization is a symmetric property.
   """
   
   def add_dependant(self,newdep):
      self.__deps.append(newdep)
      newdep.taint(taintme=True)      

   # dependency groups represent kind of one<--> many dependencies, i.e.
   # A depends from [B1,B2,...] in such a way that tainting any Bi taints A, but 
   # it is known that won't taints the other B's. On the other hand, if A is tainted
   # elsehow, all the B's may be tainted and should be marked as such.
   # An example is A = global momentum vector, Bi = slice corresponding to atom i
   def add_depgrp(self,newgrp):
      self.__depgrp.append(newgrp)

   #recursive tainting
   def taint(self,taintme=True, tainter=None):
      """Recursively sets tainted flag on dependent objects."""
      self.__tainted = True     #this is to prevent circular dependencies to make infinite recursion
      # recurse but stop when an object which is already tainted is encountered
      for item in self.__deps: 
         if (not item._depend__tainted): item.taint(tainter=self)
      # dependency groups: represent a many to one dependecy
      for grp in self.__depgrp: 
         if (not tainter in grp):
            for item in grp:
               if (not item._depend__tainted): item.taint(tainter=self)
         
      self.__tainted = taintme
      
   def tainted(self): return self.__tainted
   
   

   def __init__(self,func=None,deplist=[],value=None,name=None,):
      #print "initializing decorator object",
      self.__deps=[]
      self.__depgrp=[]
#      self.__sync=[]
      self.__func=func
      self.__name=name
      if (name is None and not func is None) : self.__name=self.__func.func_name
      self.__value=value
      if (value is None) : self.__tainted=True; #self.__unsync=True;    
      else:  self.__tainted=False; # self.__unsync=False;     
      # adds dependencies
      for item in deplist:    item.add_dependant(self)
#      for item in synclist:   item.add_synchro(self)
      
   def get(self): 
      #print "  inside decorator getter for", self.__name
      #if ((self.__tainted or self.__unsync) and not self.__func is None):
      if (self.__tainted and not self.__func is None):  
         self.__value=self.__func()
         self.taint(taintme=False, tainter=self); #self.__tainted=False; self.__unsync=False; 
      self.__tainted=False
      return self.__value
      
   def set(self,value): 
      self.taint(taintme=False, tainter=self);  #taints dependencies but not self
      self.__value=value     


## An object class which allows instantiable descriptors and 
## which provides a method to bypass the descriptor machinery
#class dobject(object):
#   def depbind(self, what, to):
#      self.getdesc(to).add_dependant(self.getdesc(what))
#      
#   def syncbind(self, o1, o2):
#      self.getdesc(o1).add_dependant(self.getdesc(o2))
#      self.getdesc(o2).add_dependant(self.getdesc(o1))
#   
#   def getdesc(self, name):
#      if name in self.__class__.__dict__ : return self.__class__.__dict__[name]
#      else: return self.__dict__[name]

#   def __getattribute__(self, name):
#      value = object.__getattribute__(self, name)
#      if hasattr(value, '__get__'): 
#         value = value.__get__(self, self.__class__)
#      return value

#   def __setattr__(self, name, value):
#      try:
#         obj = object.__getattribute__(self, name) 
#      except AttributeError:  
#         pass
#      else:
#         if hasattr(obj, '__set__'):
#            return obj.__set__(self, value)
#      return object.__setattr__(self, name, value) 

#class dep_proxy(numpy.ndarray): 
#   def __new__(cls, input_array, dep, info=None):
#      #print "__new__ called!", input_array[0]
#      # Input array is an already formed ndarray instance
#      # We first cast to be our class type
#      obj = numpy.asarray(input_array).view(cls)
#      # add the new attribute to the created instance
#      obj.info = info
#      # Finally, we must return the newly created object:
#      return obj
#      
#   def __array_finalize__(self, obj):
#      # see InfoArray.__array_finalize__ for comments
#      if obj is None: return
#      self.info = getattr(obj, 'info', None)
#      
#   def __init__(self,input_array, dep):
#      #print "init called", input_array[0], self.__getitem__(0)
#      self.__depobj=dep
#            
#   def __getitem__(self, index):
#      #print "           proxy getter ",index, super(dep_proxy,self).__getitem__(index)
#      return super(dep_proxy,self).__getitem__(index)
#      
#   def __setitem__(self, index, val):
#      #print "           proxy setter ",index,super(dep_proxy,self).__getitem__(index)
#      super(dep_proxy,self).__setitem__(index,val)
#      #self.__dep.unsync(unsyncme=False)
#      self.__depobj.taint(taintme=False, tainter=self.__depobj)


#class depend(object):
#   """ A descriptor class for quantities that may depend
#   from other quantities or be dependencies for other quantities.
#   Usage is:
#   x=depend(func=<function>,
#            deplist=<[list of dependencies]>,synclist=<[list of synced]>,
#            value=<initial value>,name=<name>)

#   Whenever one of the dependencies get tainted, then <function>
#   will be called on __get__, and a new value generated. Otherwise,
#   the cached value is returned. 
#   Synced objects are similarly cached, and marked as tainted when
#   one explicitly sets their value, but it is assumed that call to 
#   function will leave the synced objects in a consistent state. 
#   Also, it is ensured that synchronization is a symmetric property.
#   """
#   
#   def add_dependant(self,newdep):
#      self.__deps.append(newdep)

#   # dependency groups represent kind of one<--> many dependencies, i.e.
#   # A depends from [B1,B2,...] in such a way that tainting any Bi taints A, but 
#   # it is known that won't taints the other B's. On the other hand, if A is tainted
#   # elsehow, all the B's may be tainted and should be marked as such.
#   # An example is A = global momentum vector, Bi = slice corresponding to atom i
#   def add_depgrp(self,newgrp):
#      self.__depgrp.append(newgrp)

#   #recursive tainting
#   def taint(self,taintme=True, tainter=None):
#      """Recursively sets tainted flag on dependent objects."""
#      self.__tainted = True     #this is to prevent circular dependencies to make infinite recursion
#      # recurse but stop when an object which is already tainted is encountered
#      for item in self.__deps: 
#         if (not item._depend__tainted): item.taint(tainter=self)
#      # dependency groups: represent a many to one dependecy
#      for grp in self.__depgrp: 
#         if (not tainter in grp):
#            for item in grp:
#               if (not item._depend__tainted): item.taint(tainter=self)
#         
#      self.__tainted = taintme
#      
##   def unsync(self,unsyncme=True):
##      """Signals unsync with the set of synced objects."""
##      self.__unsync = unsyncme
##      for item in self.__sync: item._depend__unsync=True
#      
#   def __init__(self,func=None,deplist=[],value=None,name=None,):
#      #print "initializing decorator object",
#      self.__deps=[]
#      self.__depgrp=[]
##      self.__sync=[]
#      self.__func=func
#      self.__name=name
#      if (name is None and not func is None) : self.__name=self.__func.func_name
#      self.__value=value
#      if (value is None) : self.__tainted=True; #self.__unsync=True;    
#      else:  self.__tainted=False; # self.__unsync=False;     
#      # adds dependencies
#      for item in deplist:    item.add_dependant(self)
##      for item in synclist:   item.add_synchro(self)
#      
#   def __get__(self,obj,cls): 
#      #print "  inside decorator getter for", self.__name
#      #if ((self.__tainted or self.__unsync) and not self.__func is None):
#      if (self.__tainted and not self.__func is None):  
#         self.__value=self.__func()
#         self.taint(taintme=False, tainter=self); #self.__tainted=False; self.__unsync=False; 
#      self.__tainted=False
#      return self.__value or a proxy object if it is a ndarray instance
#      if (hasattr(self.__value,'__iter__')): 
#         return dep_proxy(self.__value,dep=self)
#      else: 
#         return self.__value
#      
#   def __set__(self,obj,value): 
#      #print "  inside decorator setter for", self.__name
#      self.taint(taintme=False, tainter=self);  #taints dependencies but not self
##      self.unsync(unsyncme=False); 
#      self.__value=value     
