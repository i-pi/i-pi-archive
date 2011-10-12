class depend(object):
   def add_dependency(self,newdep):
      self.__deps.append(newdep)

   #recursive tainting
   def taint(self):
      self.__tainted = True
      for item in self.__deps: item.taint()
      
   def __init__(self,value=None,deplist=[]):
      print "initializing decorator object",
      self.__value=value
      self.__tainted=False
      self.__deps=[]
      for item in deplist:
         item.add_dependency(self)
   
   def __call__(self, func):
      print "with function", func.func_name
      self.__func = func
      return self
      
   def __get__(self,obj,cls): 
      print "  inside decorator getter for", self.__func.func_name 
      if (self.__tainted): self.__set__(obj,value=self.__func(obj))
      return self.__value
      
   def __set__(self,obj,value): 
      print "  inside decorator setter for", self.__func.func_name
      self.taint(); self.__tainted=False  #taints dependencies but not self
      self.__value=value
      for item in self.__deps:
         item.__tainted = True   
         
class test2(object):
   # two properties which do not need to be computed
   @depend()
   def velocity(self): pass

   @depend()
   def position(self): pass
   
   #dependent properties
   @depend(0.0,[position])
   def potential(self):
      print "    getting potential"
      return self.position**4
   
   @depend(0.0,[velocity])
   def kinetic(self):
      print "    getting kinetic"
      return self.velocity**2
      
   @depend(3,[potential,kinetic])
   def energy(self):
      print "    getting energy";
      return self.kinetic+self.potential


dp=test2()
dp.position=0.1
dp.velocity=1.0
print "Getting total energy (first evaluation): "
print dp.energy
print "Getting total energy (second evaluation): "
print dp.energy

print "Hard-coded setting of kinetic energy to 4"
dp.kinetic=4.0
print "Getting total energy: "
print dp.energy

print "Now we change the position to 3"
dp.position=3.0
print "Getting total energy: "
print dp.energy

#dp.prop1.setvalue(5)

#print "dp is ", dp.prop1.getvalue()
#print "PROP2 dp is ", dp.prop2.getvalue()
#      

