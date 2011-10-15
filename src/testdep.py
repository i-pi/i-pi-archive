#from utils.depend import *
import numpy
import inspect
from utils.depend import *
            
class myobj(dobject):
   def getkin(self):
      print " [re-computing kin] ",
      k=0.0; 
      for vi in self.v: k+=vi**2
      return k/(2.0*self.mass)

   def __init__(self, varr):
      self.v=depend(name='v')
      self.v=varr
      self.mass=depend(name='mass')
      self.mass=1.0
      self.kin=depend(name='kin',func=self.getkin)
      self.getdesc('mass').add_dependant(self.getdesc('kin'))
      self.getdesc('v').add_dependant(self.getdesc('kin'))

print "qui"
vall=numpy.array(range(1,10))   
o1=myobj(vall[0:5])
o2=myobj(vall[5:10])

print "O1 v is ",   o1.v
print "O1 kin is ", o1.kin

o1.v[4]=12
print "O1 v is ",   o1.v
print "O1 kin is ", o1.kin

print "O1 kin is ", o1.kin



exit()
class dep_inst_proxy(numpy.ndarray): 
   def __new__(cls, input_array, dep):
      # print "__new__ called!", input_array[0]
      # Input array is an already formed ndarray instance
      # We first cast to be our class type
      obj = numpy.asarray(input_array).view(cls)
      # add the new attribute to the created instance
      obj.__dep = dep
      # Finally, we must return the newly created object:
      return obj
   
   def __init__(self, input_array, dep):
      #print "init called", input_array[0], self.__getitem__(0)
      self.__dep=dep
      
   def __array_finalize__(self, obj):
      # see InfoArray.__array_finalize__ for comments
      if obj is None: return
      self.__dep = getattr(obj, 'info', None)
                  
   def __getitem__(self, index):
      #print "           proxy getter ",index, super(dep_proxy,self).__getitem__(index)
      return super(dep_inst_proxy,self).__getitem__(index)
      
   def __setitem__(self, index, val):
      #print "           proxy setter ",index,super(dep_proxy,self).__getitem__(index)
      super(dep_inst_proxy,self).__setitem__(index,val)
      self.__dep.taint(taintme=False)

class dep_instanceable(object):
   def add_dep(self, dep):
      self.__deps.append(dep)
      dep.taint() # if a dependency is added, the target should be marked as tainted
      
   def taint(self, taintme=True):
      self.__tainted=taintme
      for item in self.__deps: item.taint()  
      
   def __init__(self, func=None, val=None, deps=[]):
      self.__func=func
      self.__val=val
      self.__tainted=(val is None)
      self.__deps=[]
      for item in deps: item.add_dep(self)
   
   def getter(self):
      #print "getter called"
      if (self.__tainted and not (self.__func is None)): 
         self.__val=self.__func()
         self.taint(taintme=False);
      if (self.__val.__class__ is numpy.ndarray): 
         return dep_inst_proxy(self.__val, self)
      return self.__val
   
   def setter(self, value):
      #print "setter called"
      self.__val=value
      self.taint(taintme=False)

   def deleter(self):  # we must announce that we have been removed to all the things which depend on us.
      #print "deleter called"
      self.taint(taintme=False)
      
class use_inst(object):
   def get_kin(self):
      print " [re-computing kin] ",
      k=0.0; 
      for vi in self.v.getter(): k+=vi**2
      return k/(2.0*self.mass.getter())

   def __init__(self, varr):
      self.v=dep_instanceable(None,varr,[])
      self.mass=dep_instanceable(None,1.0,[])
      self.kin=dep_instanceable(self.get_kin,None,[self.v, self.mass])

vall=numpy.array(range(1,10))
obj1=use_inst(vall[0:5])
obj2=use_inst(vall[5:10])

print "v=", obj1.v.getter()
print "K=", obj1.kin.getter()

print "v=",  obj2.v.getter()
print "K=",  obj2.kin.getter()

print "Check caching: "
print "K=", obj1.kin.getter()
print "K=", obj2.kin.getter()

print "Check change of mass: "
obj1.mass.setter(2.0)
print "K=", obj1.kin.getter()

print "Check change of v: (proxy) "
obj2.v.getter()[0]=-100
print "v=",  obj2.v.getter()
print "K=", obj2.kin.getter()

class compound(object):
   def get_tot(self):
      print " [ computing tot] ",
      return self.a.kin.getter()+self.b.kin.getter()
      
   def __init__(self, a, b):
      self.a=a
      self.b=b
      self.totk=dep_instanceable(self.get_tot,None,[a.kin,b.kin])

objc=compound(obj1,obj2)
print "Total kin.", objc.totk.getter()

obj2.v.getter()[0]=0.0
obj1.v.getter()[0]=0.0
print "Total kin.", objc.totk.getter()

# test that we can add some (quite nonsensical in this case) external depenencies on the fly
dummy=dep_instanceable(None,12,[])
dummy.add_dep(objc.totk)
print "Total kin.", objc.totk.getter()
print "Total kin.", objc.totk.getter()

dummy.setter(11)
print "Total kin.", objc.totk.getter()

exit()

#HERE BEGIN THE EXPERIMENTS...
#class test2(object):
#   @depend(1.0)
#   def half(self):
#      print "     getting half"
#      return self.twice/2

#   @depend(2.0,[],[half])
#   def twice(self):
#      print "     getting twice"
#      return self.half*2
#      
#   @depend(array([1,2]))
#   def halftwice(self):pass

#   # two properties which do not need to be computed
#   @depend()
#   def velocity(self): pass

#   @depend()
#   def position(self): pass
#   
#   #dependent properties
#   @depend(0.0,[position])
#   def potential(self):
#      print "    getting potential"
#      return self.position**4
#   
#   @depend(0.0,[velocity])
#   def kinetic(self):
#      print "    getting kinetic"
#      return self.velocity**2

#      
#   @depend(3,[potential,kinetic])
#   @depend()
#   def energy(self): pass
#      print "    getting energy";
#      return self.kinetic+self.potential

   #@depend()
   #def energy(self): pass
      
#   
#   
#   @property
#   def p(self):
#      return self._p
#   @p.setter
#   def p(self,val):
#      self._p=val
#   
#   @depend
#   def energy(self):
#      return 12; 
#   
#   def __init__(self,store):
#      print "init tests"
#      self.halftwice=store[0:3]
      
      
#class newdecorator(object):   
#   def __init__(self):
#      print "init"
#      
#   def __call__(self, func):
#      self.__func=func
#      return  self
#      
#   def __get__(self,obj,cls): 
#      print "__get__"
#      obj.__dict__['_'+obj.__class__.__name__+'__'+self.__func.func_name+'__val']=self.__func(obj)
#      return obj.__dict__['_'+obj.__class__.__name__+'__'+self.__func.func_name+'__val']

#   def __set__(self,obj,val): 
#      print "set"
#      print obj.__class__.__name__
#      obj.__dict__['_'+obj.__class__.__name__+'__'+self.__func.func_name+'__val']=val
#   
#def ndsetup(obj, func, val=None, deps=[]):
#   print '_'+obj.__class__.__name__+'__'+func+'__val'
#   obj.__dict__['_'+obj.__class__.__name__+'__'+func+'__val']=val
##   obj.__dict__.append(('_'+obj.__class__.__name__+'__'+func+'__deps',deps))
#   
#      
#class newuser(object):

#   @newdecorator()
#   def f(self):
#      print "fun", self.__f__val
#      return self.__f__val   
#   @newdecorator()
#   def f(self):
#      print "fun", self.__f__val
#      return self.__f__val
#   
#   def myget(obj,cls):
#      print "instance _get"
#      return 1
#   
#   def __init__(self):
#      ndsetup(self,'f',3.0,[])
#   


         

class newuser(object):
   def kinetic(self):
      print "recomputing twice"
      return self.half.getter()*2
         
   def __init__(self, allvec):
      self.half=depprop()
      self.twice=depprop(self.get_twice,deps=[self.half])


#velarray=numpy.array(range(1,10))
#aa=newuser(velarray[:4])
#bb=newuser(velarray[4:])

#print "aa:", aa.v, aa.kinetic
#print "bb:", bb.v, bb.kinetic

#bb.v[0]=0.0
#print "bb:", bb.v, bb.kinetic

#print "bb:", bb.v[2], bb.kinetic

#exit()
class depsweet(object):
   def addep(self, dep):
      self.__deps.append(dep)
      
   def taint(self, instance, taintme=True):
      setattr(instance,self.__name+"_tainted",taintme)
      for item in self.__deps: 
         item.taint(instance)  
      
   def __init__(self, deps=[]):
      self.__deps=[]
      for item in deps: 
         item.addep(self)

   def __call__(self, func):
      self.__name=func.func_name
      self.__func=func   
      return self

   def __set__(self, instance, val):
      print "setter called", self.__name
      setattr(instance,self.__name+"_val",val)
      self.taint(instance,taintme=False)

   def __get__(self, instance, cls):    
      print "getter called", self.__name
      if ( ( not hasattr(instance, self.__name+"_val") or getattr(instance,self.__name+"_tainted") )
         and not (self.__func is None)):
        setattr(instance,self.__name+"_val",self.__func(instance))
        self.taint(instance,taintme=False)
      return getattr(instance, self.__name+"_val")

class newuser(object):
   @depend()
   def half(self): pass
   
   @depend([half])
   def twice(self):
      print "recomputing twice"
      return self.half*2
      
#   @depstatic()
#   def v(self): pass
#   
#   @depstatic([v])
#   def kinetic(self):
#      print "computing k"
#      rk=0.0
#      for vi in self.v: rk+=vi**2
#      return rk
      

#   half=depsweet('half',func=None)
#   twice=depsweet('twice',func=get_twice,deps=[half])   

   def __init__(self, varray=[]):
      self.v=varray
     

velarray=numpy.array(range(1,10))
aa=newuser(velarray[:4])
bb=newuser(velarray[4:])

aa.half=4.
bb.half=-4
print aa.half, bb.half

exit()
@depend([aa.v])
def extkin(self):
   print "in extkin"
   return aa.v[0]+bb.v[0]

print "aa:", aa.v, aa.kinetic
print "bb:", bb.v, bb.kinetic

bb.v[0]=0.0
print "bb:", bb.v, bb.kinetic

print "bb:", bb.v[2], bb.kinetic

print "extkin:", extkin
exit()
      
   
   

#@depend()
#def ss(): pass

#store=dep_proxy([1.0,2.0,3.0,4.0],ss)
#print store[1]
#store[1]=5
#exit()


store=numpy.array([1.0,2.0,3.0,4.0,5.0,6.0])

aa=test2(store[0:3])
bb=test2(store[4:6])
aa.energy=5.0
bb.energy=8.0

print "aa", aa.energy, "bb", bb.energy

aa.p=5
bb.p=8
print "aa", aa.p, "bb", bb.p

print aa.energy is bb.energy
print aa.p is bb.p
exit()

ff=depend()
print "ecco"
dpp=[dep_proxy(store[0:4],ff),dep_proxy(store[3:6],ff)]
print dpp[0] 
print dpp[1]

print store[1], dpp[0][1], dpp[1][1]

dpp[0][1]=20.
print store[1], dpp[0][1], dpp[1][1]

exit()
dp=test2(store[0:3])
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

print "Now testing synced objects!"
print "   half: \n", dp.half, "\n   twice: \n",dp.twice

print "Setting half to 10"
dp.half=10.0
print "   half: \n", dp.half, "\n   twice: \n",dp.twice

print "Setting twice to 5"
dp.twice=5.0
print "   half: \n", dp.half, "\n   twice: \n",dp.twice

print "Setting both twice and half to 1"
dp.half=1.0
dp.twice=1.0
print "Result depends on whether we first set half"
print "   half: \n", dp.half, "\n   twice: \n",dp.twice

dp.twice=1.0
dp.half=1.0
print "...or twice"
print "   half: \n", dp.half, "\n   twice: \n",dp.twice

print "Slice references"
print  store[2], dp.halftwice[2]

store[2]=10.0
print  store[2], dp.halftwice[2]

dp.halftwice[2]=-1.0
print  store[2], dp.halftwice[2]

