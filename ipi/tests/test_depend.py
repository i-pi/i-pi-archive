"""Short test scripts.

Copyright (C) 2013, Joshua More and Michele Ceriotti

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http.//www.gnu.org/licenses/>.


Used to test the depend array view mechanism.
"""

import ipi.utils.depend as dp
import numpy as np

a = dp.depend_array(name="a",value=np.zeros((2,2),float))
b = dp.depend_array(name="b",value=np.zeros((2,2),float))

def test_slicing():
   """Depend: Slicing test"""
   c = a[0]
   print type(c)
   assert(type(c) == dp.depend_array)

def test_addition():
   """Depend: Addition test"""
   c = a + b
   print type(c)
   assert(type(c) == np.ndarray)

def test_increment():
   """Depend: Increment test"""
   c = np.zeros((2,2))
   c += a
   print type(c)
   assert(type(c) == np.ndarray)

def test_dot():
   """Depend: Dot test"""
   c = np.dot(a,b)
   print type(c)
   assert(type(c) == dp.depend_array)

def test_dotf():
   """Depend: Dot-f test"""
   rdot = np.dot
   def fdot(a,b):
      return rdot(a,b).view(np.ndarray)
   np.dot=fdot

   c = np.dot(a,b)
   assert(type(c) == np.ndarray)
   np.dot = rdot
