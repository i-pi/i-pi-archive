"""Functions used to read input configurations and print trajectories
in the JSON format.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2015 i-PI developers
# See the "licenses" directory for full license information.


import sys
import json

import numpy as np

import ipi.utils.mathtools as mt
from ipi.utils.depend import depstrip
from ipi.engine.atoms import Atoms
from ipi.engine.cell import Cell
from ipi.utils.units import Elements


__all__ = ['print_json_path', 'print_json', 'read_json', 'iter_json']


def print_json_path(beads, cell, filedesc=sys.stdout):
    """Prints all the bead configurations into a JSON formatted file.

    Prints all the replicas for each time step separately, rather than all at
    once.

    Args:
        beads: A beads object giving the bead positions.
        cell: A cell object giving the system box.
        filedesc: An open writable file object. Defaults to standard output.
    """
    raise NotImplementedError("print_json_path is not implemented yet.")


def print_json(atoms, cell, filedesc=sys.stdout, title=""):
    """Prints an atomic configuration into an XYZ formatted file.

    Args:
        atoms: An atoms object giving the centroid positions.
        cell: A cell object giving the system box.
        filedesc: An open writable file object. Defaults to standard output.
        title: This gives a string to be appended to the comment line.
    """

    a, b, c, alpha, beta, gamma = mt.h2abc_deg(cell.h)

    natoms = atoms.natoms
    # direct access to avoid unnecessary slow-down
    qs = depstrip(atoms.q)
    lab = depstrip(atoms.names)
    filedesc.write(json.dumps([natoms, a, b, c, alpha, beta, gamma, title,
                   qs.tolist(), lab.tolist()]))
    filedesc.write("\n")


def read_json(filedesc, **kwargs):
    """Reads a JSON-style file with i-pi style comments and creates an Atoms and Cell object.

    Args:
        filedesc: An open readable file object from a json formatted file with i-PI header comments.

    Returns:
        An Atoms object with the appropriate atom labels, masses and positions.
        A Cell object.
    """

    try:
        line = json.loads(filedesc.readline())
    except ValueError:
        raise EOFError("The file descriptor hit EOF.")
    atoms = Atoms(line[0])
    atoms.q = np.asarray(line[8])
    atoms.names = np.asarray(line[9], dtype='|S4')
    atoms.m = np.asarray(map(Elements.mass, atoms.names))

    a = float(line[1])
    b = float(line[2])
    c = float(line[3])
    alpha = float(line[4]) * np.pi/180
    beta = float(line[5]) * np.pi/180
    gamma = float(line[6]) * np.pi/180
    h = mt.abc2h(a, b, c, alpha, beta, gamma)
    cell = Cell(h)

    return {
        "atoms": atoms,
        "cell": cell
    }


def iter_json(filedesc):
    """Takes a json-style file and yields one Atoms object after another.

    Args:
        filedesc: An open readable file object from a json formatted file.

    Returns:
        Generator over the json trajectory, that yields
        Atoms objects with the appropriate atom labels, masses and positions.
    """

    try:
        while True:
            yield read_json(filedesc)
    except EOFError:
        pass
