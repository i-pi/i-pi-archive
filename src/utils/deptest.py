from depend import *

# creates a chain of dependencies 
#     A1 <--- B1  <--- C1
#      ^-----\    <--- C2  <----\
#     A2 <--- B2  <--- C3  <--- D1

#class test1(object):  
#   def get_b1(self):  
#      print "Getting a1", self.a1.get(),             
#      print "Updating b1",
#      return "B1"
#   def get_b2(self):  
#      print "Getting a1", self.a1.get(),
#      print "Getting a2", self.a2.get(),                          
#      print "Updating b2",
#      return "B2"
#   def get_c1(self):  
#      print "Getting b1", self.b1.get(),          
#      print "Updating c1",
#      return "C1"
#   def get_c2(self):  
#      print "Getting b1", self.b1.get(),          
#      print "Updating c2",
#      return "C2"
#   def get_c3(self):  
#      print "Getting b2", self.b2.get(),    
#      print "Updating c3",
#      return "C3"
#   def get_d1(self):  
#      print "Getting c3", self.c3.get(),
#      print "Getting c2", self.c2.get(),      
#      print "Updating d1",
#      return "D1"
#      
#   def __init__(self): 
#      self.a1=depend_value(value="A1",name="A1")
#      self.a2=depend_value(value="A2",name="A2")      
#      self.b1=depend_calc(func=self.get_b1,name="B1")
#      self.b2=depend_calc(func=self.get_b2,name="B2")
#      self.c1=depend_calc(func=self.get_c1,name="C1")                  
#      self.c2=depend_calc(func=self.get_c2,name="C2")
#      self.c3=depend_calc(func=self.get_c3,name="C3")
#      self.d1=depend_calc(func=self.get_d1,name="D1")
#      self.a1.add_dependant(self.b1)   
#      self.a1.add_dependant(self.b2)
#      self.a2.add_dependant(self.b2)        
#      self.c3.add_dependant(self.d1)      
#      self.c2.add_dependant(self.d1)
#      self.b1.add_dependant(self.c1)      
#      self.b1.add_dependant(self.c2)      
#      self.b2.add_dependant(self.c3)   

   
class multi(object):

   def multiply(self):
      print "@@ MULTIPLYING"
      return self.base.get()*self.multi.get()
   
   def __init__(self, multi, array):
      self.multi=depend_value(name="multi",value=multi)      
      self.base=depend_array(array)
      self.mult=depend_array(numpy.zeros(array.size), deps=depend_func(func=self.multiply))      
      self.multi.deps.add_dependant(self.mult.deps)
      self.base.deps.add_dependant(self.mult.deps)
      
   
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

import numpy, gc, sys, types

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
exit()

#no = get_refcounts()
#print no


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

