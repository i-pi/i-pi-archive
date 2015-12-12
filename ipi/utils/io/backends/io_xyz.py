"""Functions used to read input configurations and print trajectories
in the XYZ format.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import sys
import re

import numpy as np

import ipi.utils.mathtools as mt
from ipi.utils.depend import depstrip
from ipi.engine.atoms import Atoms
from ipi.engine.cell import Cell
from ipi.utils.units import Elements


__all__ = ['print_xyz_path', 'print_xyz', 'read_xyz', 'iter_xyz']


def print_xyz_path(beads, cell, filedesc=sys.stdout):
    """Prints all the bead configurations into a XYZ formatted file.

    Prints all the replicas for each time step separately, rather than all at
    once.

    Args:
        beads: A beads object giving the bead positions.
        cell: A cell object giving the system box.
        filedesc: An open writable file object. Defaults to standard output.
    """

    a, b, c, alpha, beta, gamma = mt.h2abc_deg(cell.h)

    fmt_header = "%d\n# bead: %d CELL(abcABC): %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f \n"
    natoms = beads.natoms
    nbeads = beads.nbeads
    for j in range(nbeads):
        filedesc.write(fmt_header % (natoms, j, a, b, c, alpha, beta, gamma))
        for i in range(natoms):
            qs = depstrip(beads.q)
            lab = depstrip(beads.names)
            filedesc.write("%8s %12.5e %12.5e %12.5e\n" % (lab[i], qs[j][3*i], qs[j][3*i+1], qs[j][3*i+2]))


def print_xyz(atoms, cell, filedesc=sys.stdout, title=""):
    """Prints an atomic configuration into an XYZ formatted file.

    Args:
        atoms: An atoms object giving the centroid positions.
        cell: A cell object giving the system box.
        filedesc: An open writable file object. Defaults to standard output.
        title: This gives a string to be appended to the comment line.
    """

    a, b, c, alpha, beta, gamma = mt.h2abc_deg(cell.h)

    natoms = atoms.natoms
    fmt_header = "%d\n# CELL(abcABC): %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %s\n"
    filedesc.write(fmt_header % (natoms, a, b, c, alpha, beta, gamma, title))
    # direct access to avoid unnecessary slow-down
    qs = depstrip(atoms.q)
    lab = depstrip(atoms.names)
    for i in range(natoms):
        filedesc.write("%8s %12.5e %12.5e %12.5e\n" % (lab[i], qs[3*i], qs[3*i+1], qs[3*i+2]))


def read_xyz(filedesc, **kwargs):
    """Readss an XYZ-style file with i-pi style comments and creates an Atoms and Cell object

    Args:
        filedesc: An open readable file object from a xyz formatted file with i-PI header comments.

    Returns:
        An Atoms object with the appropriate atom labels, masses and positions.
        A Cell object.
    """

    natoms = filedesc.readline()
    if natoms == "":
        raise EOFError("The file descriptor hit EOF.")
    natoms = int(natoms)
    comment = filedesc.readline()
    reabc = re.compile('# CELL.abcABC.: ([-0-9.Ee ]*) ').search(comment)
    regenh = re.compile('# CELL.GENH.: ([-0-9.Ee ]*)').search(comment)
    reh = re.compile('# CELL.H.: ([-0-9.Ee ]*)').search(comment)
    usegenh = False
    if reabc is not None:
        a, b, c, alpha, beta, gamma = reabc.group(1).split()
        a = float(a)
        b = float(b)
        c = float(c)
        alpha = float(alpha) * np.pi/180
        beta = float(beta) * np.pi/180
        gamma = float(gamma) * np.pi/180
        h = mt.abc2h(a, b, c, alpha, beta, gamma)
    elif reh is not None:
        h = np.array(reh.group(1).split(), float)
        h.resize((3,3))
    elif regenh is not None:
        genh = np.array(regenh.group(1).split(), float)
        genh.resize((3,3))
        invgenh = np.linalg.inv(genh)
        # convert back & forth from abcABC representation to get an upper triangular h
        h = mt.abc2h(*mt.genh2abc(genh))
        usegenh = True
    else:
        # defaults to unit box
        h = mt.abc2h(1.0, 1.0, 1.0, np.pi/2, np.pi/2, np.pi/2)
    cell = Cell(h)

    qatoms = []
    names = []
    masses = []
    iat = 0
    while (iat < natoms):
        body = filedesc.readline()
        if body.strip() == "":
            break
        body = body.split()
        name = body[0]
        names.append(name)
        masses.append(Elements.mass(name))
        x = float(body[1])
        y = float(body[2])
        z = float(body[3])

        if usegenh:
            # must convert from the input cell parameters to the internal convention
            u = np.array([x,y,z])
            us = np.dot(u, invgenh)
            u = np.dot(h, us)
            x, y, z = u

        qatoms.append(x)
        qatoms.append(y)
        qatoms.append(z)
        iat += 1

    if natoms != len(names):
        raise ValueError("The number of atom records does not match the header of the xyz file.")

    atoms = Atoms(natoms)
    atoms.q = np.asarray(qatoms)
    atoms.names = np.asarray(names, dtype='|S4')
    atoms.m = np.asarray(masses)

    return {
        "atoms": atoms,
        "cell": cell,
    }


def iter_xyz(filedesc):
    """Takes a xyz-style file and yields one Atoms object after another.

    Args:
        filedesc: An open readable file object from a xyz formatted file.

    Returns:
        Generator over the xyz trajectory, that yields
        Atoms objects with the appropriate atom labels, masses and positions.
    """

    try:
        while 1:
            yield read_xyz(filedesc)
    except EOFError:
        pass
