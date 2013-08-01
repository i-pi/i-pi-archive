import utils.depend as dp
import numpy as np

print "## Creation test"
a=dp.depend_array(name="a",value=np.zeros((2,2),float))
b=dp.depend_array(name="b",value=np.zeros((2,2),float))

print "## Slicing test"
c=a[0]
print type(c)

print "## Addition test"
c=a+b
print type(c)

print "## Increment test"
c=np.zeros((2,2)); c+=a
print type(c)

print "## Dot test"
c=np.dot(a,b)
print type(c)

rdot=np.dot
def fdot(a,b):
   return rdot(a,b).view(np.ndarray)
#np.dot=fdot

print "## Dot-f test"
c=np.dot(a,b)

print "## check bases"
c=a[0:1]
c.printbases()
print type(c)
