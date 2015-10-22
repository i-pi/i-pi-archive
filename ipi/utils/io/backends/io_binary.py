"""Functions used to print the trajectories and read input configurations
(or even full status dump) as unformatted binary.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import sys
import numpy as np


__all__ = ['print_bin']


def print_bin(atoms, cell, filedesc=sys.stdout, title=""):
    """Prints an atomic configuration into a binary file.

    Args:
        beads: An atoms object giving the centroid positions.
        cell: A cell object giving the system box.
        filedesc: An open writable file object. Defaults to standard output.
        title: This gives a string to be appended to the comment line.
    """

    buff = filedesc
    cell.h.tofile(buff)
    nat = np.asarray([atoms.natoms])
    nat.tofile(buff)
    atoms.names.tofile(buff)
    atoms.q.tofile(buff)
