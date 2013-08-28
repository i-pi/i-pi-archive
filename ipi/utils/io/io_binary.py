"""Contains the functions used to print the trajectories and read input
configurations (or even full status dump) as unformatted binary.

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


Functions:
   print_bin: Prints an atomic configuration.
"""

__all__ = ['print_bin']

import os
import numpy as np
import math, sys
from ipi.utils.depend import depstrip

def print_bin(atoms, cell, filedesc = sys.stdout, title=""):
   """Prints the centroid configurations, into a binary file.

   Args:
      beads: An atoms object giving the centroid positions.
      cell: A cell object giving the system box.
      filedesc: An open writable file object. Defaults to standard output.
      title: This gives a string to be appended to the comment line.
   """

   buff = filedesc # .buffer
   cell.h.tofile(buff)
   nat = np.asarray([atoms.natoms])
   nat.tofile(buff)
   atoms.names.tofile(buff)
   atoms.q.tofile(buff)

