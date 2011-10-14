import numpy
class dep_proxy(numpy.ndarray): 
   def __new__(cls, input_array, instance, dep, info=None):
      #print "__new__ called!", input_array[0]
      # Input array is an already formed ndarray instance
      # We first cast to be our class type
      obj = numpy.asarray(input_array).view(cls)
      # add the new attribute to the created instance
      obj.info = info
      # Finally, we must return the newly created object:
      return obj
      
   def __array_finalize__(self, obj):
      # see InfoArray.__array_finalize__ for comments
      if obj is None: return
      self.info = getattr(obj, 'info', None)
      
   def __init__(self,input_array,instance, dep):
      #print "init called", input_array[0], self.__getitem__(0)
      self.__dep=dep
      self.__inst=instance
            
   def __getitem__(self, index):
      #print "           proxy getter ",index, super(dep_proxy,self).__getitem__(index)
      return super(dep_proxy,self).__getitem__(index)
      
   def __setitem__(self, index, val):
      #print "           proxy setter ",index,super(dep_proxy,self).__getitem__(index)
      super(dep_proxy,self).__setitem__(index,val)
      self.__dep.unsync(self.__inst,unsyncme=False)
      self.__dep.taint(self.__inst,taintme=False)
      
   def add_depend(self, dep):
      print "dep added through proxy", dep.__class__
      self.__dep.add_depend(dep)

#TODO: (maybe)
# * Make list of dependencies unique
# * Circular dependencies detection and avoidance
# * Type check to make sure that dependencies are depend objects
class depend(object):
   """ A decorator class for quantities that may depend
   from other quantities or be dependencies for other quantities.
   Usage is:
   @depend(value=<initial value>,deplist=<[list of dependencies]>,synclist=<[list of synced]>)
   def function(self):
   Whenever one of the dependencies get tainted, then function
   will be called on __get__, and a new value generated. Otherwise,
   the cached value is returned. 
   Synced objects are similarly cached, and marked as tainted when
   one explicitly sets their value, but it is assumed that call to 
   function will leave the synced objects in a consistent state. 
   Also, it is ensured that synchronization is a symmetric property.
   """
   
   def add_depend(self, dep):
      print "dep added", dep.__class__, self.__deps
      self.__deps.append(dep)

   def add_synchro(self,dep):
      self.__sync.append(dep)
      # makes sure that syncing is symmetric
      if (not(self in dep._depend__sync)): dep.add_synchro(self) 
            
   def taint(self, instance, taintme=True):
      """Recursively sets tainted flag on dependent objects."""   
      setattr(instance,self.__name+"_tainted",taintme)
      for item in self.__deps: 
         item.taint(instance)  
      
   def unsync(self,instance,unsyncme=True):
      """Signals unsync with the set of synced objects."""
      setattr(instance,self.__name+"_unsync",unsyncme)
      for item in self.__sync: item._depend__unsync=True
               
   def __call__(self, func):
      print "with function", func.func_name
      self.__func = func
      return self
                     
   def __init__(self, deps=[], sync=[]):
      self.__deps=[]; self.__sync=[]
      for item in deps: item.add_depend(self)
      for item in sync: item.add_synchro(self)

   def __call__(self, func):
      self.__name=func.func_name
      self.__func=func   
      return self

   def __set__(self, instance, val):
      print "setter called", self.__name
      setattr(instance,self.__name+"_val",val)
      self.taint(instance,taintme=False)
      self.unsync(instance,unsyncme=False)

   def __get__(self, instance, cls):    
      print "getter called", self.__name
      if ( ( not hasattr(instance, self.__name+"_val") 
            or getattr(instance,self.__name+"_tainted")
            or getattr(instance,self.__name+"_unsync") )
            and not (self.__func is None)):
        setattr(instance,self.__name+"_val",self.__func(instance))
        self.taint(instance,taintme=False);   self.unsync(instance,unsyncme=False)
        
      if (hasattr(getattr(instance, self.__name+"_val"),'__iter__')): 
         return dep_proxy(getattr(instance, self.__name+"_val"),instance=instance,dep=self)
      else: 
         return getattr(instance, self.__name+"_val")        
#         
#   def __getattribute__(self, name):
#      print "accessing %r.%s" % (self, name)
#      return object.__getattribute__(self, name)
#   def __getattr__(self, name):
#      print "access %r.%s" % (self, name)
#      return object.__getattribute__(self, name)         
         
      
      
#!THIS DOESN'T REALLY WORK
#it only affects the elements of the class and not those of its instances. pretty useless
class depstatic(object):
   """ A decorator class for quantities that may depend
   from other quantities or be dependencies for other quantities.
   Usage is:
   @depend(value=<initial value>,deplist=<[list of dependencies]>,synclist=<[list of synced]>)
   def function(self):
   Whenever one of the dependencies get tainted, then function
   will be called on __get__, and a new value generated. Otherwise,
   the cached value is returned. 
   Synced objects are similarly cached, and marked as tainted when
   one explicitly sets their value, but it is assumed that call to 
   function will leave the synced objects in a consistent state. 
   Also, it is ensured that synchronization is a symmetric property.
   """
   
   def add_dependency(self,newdep):
      self.__deps.append(newdep)
      
   def add_synchro(self,newdep):
      self.__sync.append(newdep)
      # makes sure that syncing is symmetric
      if (not(self in newdep._depend__sync)): newdep.add_synchro(self) 
      
   #recursive tainting
   def taint(self,taintme=True):
      """Recursively sets tainted flag on dependent objects."""
      self.__tainted = taintme
      for item in self.__deps: item.taint()
      
   def unsync(self,unsyncme=True):
      """Signals unsync with the set of synced objects."""
      self.__unsync = unsyncme
      for item in self.__sync: item._depend__unsync=True
      
   def __init__(self,func=None,value=None,deplist=[],synclist=[]):
      print "initializing decorator object",
      print self
      self.__value=value
      self.__tainted=False;   self.__deps=[]
      self.__unsync=False;    self.__sync=[]
      self.__func=func
      for item in deplist:    item.add_dependency(self)
      for item in synclist:   item.add_synchro(self)

   def __call__(self, func):
      print "with function", func.func_name
      self.__func = func
      return self
      
   def __get__(self,obj,cls): 
      print "  inside decorator getter for", self.__func.func_name 
      if (self.__tainted or self.__unsync): 
         self.__value=self.__func(obj)
         self.taint(); self.__tainted=False; self.__unsync=False; 
      #return self.__value
      if (hasattr(self.__value,'__iter__')): 
         return dep_proxy(self.__value,dep=self)
      else: 
         return self.__value
      
   def __set__(self,obj,value): 
      print "  inside decorator setter for", self.__func.func_name
      self.taint(taintme=False);  #taints dependencies but not self
      self.unsync(unsyncme=False); 
      self.__value=value     
