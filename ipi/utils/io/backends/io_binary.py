"""Functions used to print the trajectories and read input configurations
(or even full status dump) as unformatted binary.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import sys
import numpy as np


__all__ = ['print_binary', 'read_binary']


def print_binary(atoms, cell, filedesc=sys.stdout, title=""):
    """Prints an atomic configuration into a binary file.

    Args:
        beads: An atoms object giving the centroid positions.
        cell: A cell object giving the system box.
        filedesc: An open writable file object. Defaults to standard output.
        title: This gives a string to be appended to the comment line.
    """

    buff = filedesc    
    cell.h.flatten().tofile(buff)
    nat = np.asarray([atoms.natoms])
    nat.tofile(buff)
    atoms.q.tofile(buff)
    names = "|".join(atoms.names)  # concatenates names, assuming none contains '|'
    nat[0] = len(names)
    nat.tofile(buff)
    np.asarray([names]).tofile(buff)


def read_binary(filedesc, **kwarg):
    try: 
        cell = np.fromfile(filedesc, dtype=float, count=9)
        cell.shape = (3,3)
        nat = np.fromfile(filedesc, dtype=int, count=1)[0]
        qatoms = np.fromfile(filedesc, dtype=float, count=3*nat)
        nat = np.fromfile(filedesc, dtype=int, count=1)[0]
        names = np.fromfile(filedesc, dtype='|S1', count=nat)    
        names = "".join(names)
        names = names.split('|')
        masses = np.zeros(len(names))
    except (StopIteration,ValueError):
        raise EOFError    
    return ("", cell, qatoms, names, masses)
