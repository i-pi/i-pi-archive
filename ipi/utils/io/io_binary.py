"""Contains the functions used to print the trajectories and read input
configurations (or even full status dump) as unformatted binary.

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

