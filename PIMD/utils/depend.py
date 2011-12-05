import numpy

class depend(object):
   """ A descriptor class for quantities that may depend
   from other quantities or be dependencies for other quantities.
   Contains: __deps = all objects dependent on self, __depgrps = 
   many to one type dependent groups of self, __tainted = tainted flag,
   __func = function to recalculate object if tainted, name = object name,
   __value = object value

   Initialised by: obj = depend(func, deplist, value, name)
   func = function to recalculate object if tainted, default = None
   deplist = list of objects dependent on obj, default = []
   value = initial value of obj, default = None
   name = object name, default = None

   Get and set must be called explicitly, and it is possible to 
   directly taint the object.
   Whenever obj get tainted, then __func will be called upon calling obj.get(), 
   and a new value generated. Otherwise, the cached value is returned."""
   
   def add_dependant(self,newdep):
      """Makes newdep dependent on self"""

      self.__deps.append(newdep)
      newdep.taint(taintme=True)      

   def add_depgrp(self,newgrp):
      """Makes a group of objects dependent on self. These are used for 
         one<-->many dependencies, for example the system position vector is
         dependent on all the atom position vectors. Tainting one of the 
         group will only taint self, whereas tainting self will taint all of
         the group, as any or all of these objects may have been changed"""

      self.__depgrp.append(newgrp)

   def taint(self,taintme=True, tainter=None):
      """Recursively sets tainted flag on dependent objects."""
      self.__tainted = True     #this is to prevent circular dependencies
      for item in self.__deps: 
         if (not item.tainted()):
            item.taint(tainter=self)
      for grp in self.__depgrp: 
         if (not tainter in grp):
            for item in grp:
               if (not item.tainted()):
                  item.taint(tainter=self)
         
      self.__tainted = taintme
      
   def tainted(self):
      """Returns tainted flag"""

      return self.__tainted

   def __init__(self,func=None,deplist=[],value=None,name=None):
      self.__deps=[]
      self.__depgrp=[]
      self.__func=func
      self.name=name
      if (name is None and not func is None):
         self.name=self.__func.func_name
      self.__value=value
      if (value is None): 
         self.__tainted=True
      else:
         self.__tainted=False
      for item in deplist:
         item.add_dependant(self)
      
   def get(self):
      """Returns a reference to the object's value after recalculating it
         if the object has been tainted"""
 
      if (self.__tainted and not self.__func is None):  
         self.__value=self.__func()
         self.taint(taintme=False, tainter=self)
      self.__tainted=False
      return self.__value

   def get_array(self, getter=None):
      if (self.__tainted and not self.__func is None):
         self.__value[:]=self.__func()
         self.taint(taintme=False, tainter=self)
      elif (self.__tainted and not self.__depgrp == []):
         for grp in self.__depgrp:
            if (not getter in grp):
               for item in grp:
                  item.get_array(getter=self)
            
      self.__tainted=False
      return self.__value
      
   def set(self,value): 
      """Changes value, and taints dependents"""

      self.taint(taintme=False, tainter=self)
      self.__value=value     

