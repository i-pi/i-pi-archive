import numpy
# An object class which allows instantiable descriptors and 
# which provides a method to bypass the descriptor machinery
class dobject(object):
    def getdescriptor(self, name):
        if name in self.__class__.__dict__ : return self.__class__.__dict__[name]
        else: return self.__dict__[name]

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

class dep_proxy(numpy.ndarray): 
   def __new__(cls, input_array, dep, info=None):
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
      
   def __init__(self,input_array, dep):
      #print "init called", input_array[0], self.__getitem__(0)
      self.__dep=dep
            
   def __getitem__(self, index):
      #print "           proxy getter ",index, super(dep_proxy,self).__getitem__(index)
      return super(dep_proxy,self).__getitem__(index)
      
   def __setitem__(self, index, val):
      #print "           proxy setter ",index,super(dep_proxy,self).__getitem__(index)
      super(dep_proxy,self).__setitem__(index,val)
      self.__dep.unsync(unsyncme=False)
      self.__dep.taint(taintme=False)


class depend(object):
   """ A descriptor class for quantities that may depend
   from other quantities or be dependencies for other quantities.
   Usage is:
   x=depend(func=<function>,
            deplist=<[list of dependencies]>,synclist=<[list of synced]>,
            value=<initial value>,name=<name>)

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
      
   def __init__(self,func=None,deplist=[],synclist=[],value=None,name=None,):
      #print "initializing decorator object",
      self.__deps=[]
      self.__sync=[]
      self.__func=func
      self.__name=name
      if (name is None and not func is None) : self.__name=self.__func.func_name
      self.__value=value
      if (value is None) : self.__tainted=True; self.__unsync=True;    
      else:  self.__tainted=False;  self.__unsync=False;     
      # adds dependencies
      for item in deplist:    item.add_dependant(self)
      for item in synclist:   item.add_synchro(self)
      
   def __get__(self,obj,cls): 
      #print "  inside decorator getter for", self.__name
      if ((self.__tainted or self.__unsync) and not self.__func is None): 
         self.__value=self.__func()
         self.taint(); self.__tainted=False; self.__unsync=False; 
      #return self.__value or a proxy object if it is a ndarray instance
      if (hasattr(self.__value,'__iter__')): 
         return dep_proxy(self.__value,dep=self)
      else: 
         return self.__value
      
   def __set__(self,obj,value): 
      #print "  inside decorator setter for", self.__name
      self.taint(taintme=False);  #taints dependencies but not self
      self.unsync(unsyncme=False); 
      self.__value=value     
