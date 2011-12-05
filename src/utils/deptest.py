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
      
   
import numpy

ee=multi(4, numpy.zeros(10))

ee.base[:]=1.0

print "getting once", ee.mult[:]
print "getting twice", ee.mult[:]

ee.base[1]=0.0
print "getting once", ee.mult[:]
print "getting twice", ee.mult[:]

ee.multi.set(5)
ss=ee.base[4:8]
ss[0]=12
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

dee=dmulti(4.0,numpy.zeros(5))
print "getting once", dee.mult

dee.base=10.0
print "getting twice", dee.mult



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

print ee.synctest.mul
ee.synctest.mul=6
print "mul", ee.synctest.mul
print "div", ee.synctest.div



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

