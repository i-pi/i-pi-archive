from depend import *
import timeit

class testdep(object):

   def do_sumval(self):
      tmp=0.0
      for i in self.val.get():
         tmp+=i
      return tmp

   def do_sumrep(self):
      tmp=0.0
      for i in self.rep.get():
         tmp+=i
      return tmp

      
   def v2r(self):
      tmp=numpy.array(self.val.get(),copy=True) # creates a temporary
      for ind,val in enumerate(tmp): 
         tmp[ind]=self.factor.get()/val
      
      return tmp[:]
   
   def r2v(self):
      tmp=numpy.array(self.rep.get(),copy=True) # creates a temporary
      for ind,val in enumerate(tmp): 
         tmp[ind]=self.factor.get()/val
      return tmp[:]

   def __init__(self, factor, array):
      self.factor=depend_value(name="factor",value=factor)   

      self.sumval=depend_value(name="sumval",deps=depend_func(func=self.do_sumval))
      self.sumrep=depend_value(name="sumrep",deps=depend_func(func=self.do_sumrep))
      
      sync=synchronizer()
      self.rep=depend_array(name="rep",input_array=numpy.zeros(10,float), deps=depend_sync(func={  "val" : self.v2r },synchro=sync))
      self.val=depend_array(name="val",input_array=numpy.zeros(10,float), deps=depend_sync(func={  "rep" : self.r2v },synchro=sync))
            
      self.factor.deps.add_dependant(self.rep.deps)
      self.factor.deps.add_dependant(self.val.deps)
      
      self.val.deps.add_dependant(self.sumval.deps)
      self.rep.deps.add_dependant(self.sumrep.deps)      
   
class dtestdep(dobject):

   def do_sumval(self):
      tmp=0.0
      for i in self.val:
         tmp+=i
      return tmp

   def do_sumrep(self):
      tmp=0.0
      for i in self.rep:
         tmp+=i
      return tmp

      
   def v2r(self):
      tmp=numpy.array(self.val,copy=True) # creates a temporary
      for ind,val in enumerate(tmp): 
         tmp[ind]=self.factor/val
      
      return tmp[:]
   
   def r2v(self):
      tmp=numpy.array(self.rep,copy=True) # creates a temporary
      for ind,val in enumerate(tmp): 
         tmp[ind]=self.factor/val
      return tmp[:]

   def __init__(self, factor, array):
      self.factor=depend_value(name="factor",value=factor)   

      self.sumval=depend_value(name="sumval",deps=depend_func(func=self.do_sumval))
      self.sumrep=depend_value(name="sumrep",deps=depend_func(func=self.do_sumrep))
      
      sync=synchronizer()
      self.rep=depend_array(name="rep",input_array=numpy.zeros(10,float), deps=depend_sync(func={  "val" : self.v2r },synchro=sync))
      self.val=depend_array(name="val",input_array=numpy.zeros(10,float), deps=depend_sync(func={  "rep" : self.r2v },synchro=sync))
            
      depget(self,"factor").add_dependant(depget(self,"rep"))
      depget(self,"factor").add_dependant(depget(self,"val"))
      
      depget(self,"val").add_dependant(depget(self,"sumval"))
      depget(self,"rep").add_dependant(depget(self,"sumrep"))

def test_1():
   td1=testdep(2.0,numpy.zeros(10))      

   erep=td1.rep[0:5]
   td1.val[:]=numpy.array(range(10),float)
   td1.val[:]+=1.0
   for i in range(1000):
      grep=td1.val[3:6]
      grep[1:2]-=td1.sumval.get()/10
      erep[:]+=(td1.sumval.get()-td1.sumrep.get()*2)/100

def test_2():
   td1=dtestdep(2.0,numpy.zeros(10))      

   erep=td1.rep[0:5]
   td1.val=numpy.array(range(10),float)
   td1.val+=1.0
   for i in range(1000):
      grep=td1.val[3:6]
      grep[1:2]-=td1.sumval/10
      erep[:]+=(td1.sumval-td1.sumrep*2)/100

a = timeit.Timer("test_1()", "gc.enable(); from __main__ import test_1")
      
zeroth = a.timeit(number = 100)

print "Timing for test1: ", zeroth

a = timeit.Timer("test_2()", "gc.enable(); from __main__ import test_2")
      
zeroth = a.timeit(number = 100)

print "Timing for test2: ", zeroth
exit()
      
import numpy, gc, sys, types
def get_refcounts():
   d = {}
   sys.modules
   for m in sys.modules.values():   
      for sym in dir(m):
         o = getattr(m, sym)
         if type(o) is types.ClassType:
            d[o] = sys.getrefcount(o)

   pairs = map (lambda x: (x[1], x[0]), d.items())
   pairs.sort()
   pairs.reverse()
   return pairs

print gc.garbage
no = get_refcounts()
print "REFCOUNTS ", no


exit();


def type_1():
   ee=multi(4, numpy.zeros(10))
   
   ee.base[:]=1.0
   
   
   
   print "getting once", ee.mult[:]
   print "getting twice", ee.mult[:]
   
   ee.base[1]=0.0
   print "getting once", ee.mult[:]
   print "getting twice", ee.mult[:]
   
   ee.multi.set(5)
   ss=ee.base[4:8]
   tt=ss[0:2]
   tt[0]=12
   print "getting once", ee.mult[:]
   print "getting twice", ee.mult[:]
 
class dmulti(dobject):

   def multiply(self):
      print "@@ MULTIPLYING"
      return self.base*self.multi
   
   def __init__(self, multi, array):
      self.multi=depend_value(name="multi",value=multi)
      self.base=depend_array(array)
      self.mult=depend_array(numpy.zeros(array.size), deps=depend_func(func=self.multiply))      
      depget(self,"multi").add_dependant(depget(self,"mult"))
      depget(self,"base").add_dependant(depget(self,"mult"))

def type_2():
   dee=dmulti(4, numpy.zeros(10))
   
   dee.base=1.0
   
   
   
   print "getting once", dee.mult
   print "getting twice", dee.mult
   
   dee.base[1] =0.0
   print "getting once", dee.mult
   print "getting twice", dee.mult
   
   dee.multi = 5
   dss = dee.base[4:8]
   dtt = dss[0:2]
   dtt[0]=12.0
   print "getting once", dee.mult
   print "getting twice", dee.mult

class scalar(object):
   def multiply(self):
      return self.a.get()*self.b.get()

   def add(self):
      return (self.a.get() + self.b.get())*numpy.ones(12)

   def __init__(self):
      self.a = depend_value(name="a", value = 12)
      self.b = depend_value(name="b", value = 13)
      self.c = depend_value(name="c", value = 0, deps = depend_func(func = self.multiply))
      self.a.deps.add_dependant(self.c.deps)
      self.b.deps.add_dependant(self.c.deps)

def type_4():
   ee = scalar()
   print ee.a.get(), ee.b.get(), ee.c.get()

   ee.a.set(10)
   print ee.a.get(), ee.b.get(), ee.c.get()

class arr(dobject):
   def __init__(self):
      self.d = depend_array(numpy.ones(12))

def type_3():
   ee = arr()
   print
   print ee.d
   del ee.d.deps
   

import timeit

#print type_1()
#print
#print type_2()

#exit()

a = timeit.Timer("type_3()", "gc.enable(); from __main__ import type_3")
zeroth = a.timeit(number = 1)
#zeroth = a.timeit(number = 580000)
#exit()
type_4()

a = timeit.Timer("type_1()", "gc.enable(); from __main__ import type_1")
first = a.timeit(number = 1)
a = timeit.Timer("type_2()", "gc.enable(); from __main__ import type_2")
second = a.timeit(number = 1)
print first, second



for i in range(100):
   type_2()
import gc

print gc.garbage
no = get_refcounts()
print "REFCOUNTS ", no


exit(1)

class multisync(object):
   def multiply(self):
      print "@@ MULTIPLYING"
      return self.frac.get()*self.multi.get()
   def divide(self):
      print "@@ DIVIDING"
      return self.mult.get()/self.multi.get()
        
   def __init__(self, factor, size):
      self.multi=depend_value(name="multi",value=factor)
      sync=synchronizer()
      self.mult=depend_value(name="mult",value=0, deps=depend_sync(func={  "frac" : self.multiply },name="mult",synchro=sync))
      self.frac=depend_value(name="frac",value=0, deps=depend_sync(func={  "mult" : self.divide },name="frac",synchro=sync))
      
      self.multi.deps.add_dependant(self.mult.deps)
      self.multi.deps.add_dependant(self.frac.deps)
      
      
see=multisync(2.0, 10)

print "@@@ SYNCOBJECT"
see.frac.set(1.0)

print "multi ", see.mult.get()
print "multi ", see.mult.get()
print "fract ", see.frac.get()

see.mult.set(10.0)
print "multi ", see.mult.get()
print "fract ", see.frac.get()


see.multi.set(5.0)
print "multi ", see.mult.get()
print "fract ", see.frac.get()
      
exit(1)
ee=multi(4)
ee.multi.set(2)

exit()

for i in range(10000):
   type_1()
import gc
print gc.garbage
exit()

for i in range(10000):
   type_2()

exit(1)
print t1.d1.get()
t1.a1.set("manb2");
print t1.d1.get()

#now we break the thing

# creates a chain of dependencies 
#     A1 <--- B1  <--- C1
#      ^-----\    <--- C2  <----\
#     A2 <--- B2  <--- C3  <--- D1

t1.a1.taint()
print t1.b1.get()

t1.a2.taint()
print t1.d1.get()

