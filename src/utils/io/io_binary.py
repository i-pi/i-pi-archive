"""Contains the functions used to print the trajectories and read input
configurations (or even full status dump) as unformatted binary

Functions:
   print_bin: Prints an atomic configuration.
   read_bin:  Reads the cell parameters and atom configurations from a xyz file.
"""

__all__ = ['print_bin', 'read_bin']

import os
import numpy as np
import math, sys
from utils.depend import depstrip

def print_bin(atoms, cell, filedesc = sys.stdout, title=""):
   buff=filedesc # .buffer 
   cell.h.tofile(buff)
   nat=np.asarray([atoms.natoms])
   nat.tofile(buff)
   atoms.names.tofile(buff)
   atoms.q.tofile(buff)

