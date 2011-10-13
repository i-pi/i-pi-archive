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
   
   def add_dependency(self,newdep):
      self.__deps.append(newdep)
      
   def add_synchro(self,newdep):
      self.__sync.append(newdep)
      # makes sure that syncing is symmetric
      if (not(self in newdep._depend__sync)): newdep.add_synchro(self) 
      
   #recursive tainting
   def taint(self):
      """Recursively sets tainted flag on dependent objects."""
      self.__tainted = True
      for item in self.__deps: item.taint()
      
   def __init__(self,value=None,deplist=[],synclist=[]):
      print "initializing decorator object",
      self.__value=value
      self.__tainted=False;   self.__deps=[]
      self.__unsync=False;    self.__sync=[]
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
      return self.__value
      
   def __set__(self,obj,value): 
      print "  inside decorator setter for", self.__func.func_name
      self.taint(); self.__tainted=False  #taints dependencies but not self
      self.__value=value
      self.__unsync=False; 
      for item in self.__sync: item._depend__unsync=True
