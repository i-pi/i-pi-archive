"""Short test scripts.

Copyright (C) 2013, Joshua More and Michele Ceriotti

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


Used to test the depend array view mechanism.
"""

import sys
sys.path.append("../")
sys.path.append("../../")

import utils.depend as dp
import numpy as np

print "## Creation test"
a = dp.depend_array(name="a",value=np.zeros((2,2),float))
b = dp.depend_array(name="b",value=np.zeros((2,2),float))

print "## Slicing test"
c = a[0]
print type(c)

print "## Addition test"
c = a + b
print type(c)

print "## Increment test"
c = np.zeros((2,2))
c += a
print type(c)

print "## Dot test"
c = np.dot(a,b)
print type(c)

rdot = np.dot
def fdot(a,b):
   return rdot(a,b).view(np.ndarray)
#np.dot=fdot

print "## Dot-f test"
c = np.dot(a,b)
