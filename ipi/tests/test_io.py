"""Deals with testing the io system.
"""

import filecmp
import os
import numpy as np
from numpy.testing import assert_equal
from common import local

from engine.cell import Cell

from utils.io import io_xyz
from utils.io import io_pdb

pos = np.array([i for i in range(3*3)])

def test_read_xyz():
    with open(local("test.pos_0.xyz"), "r") as f:
        atoms = io_xyz.read_xyz(f)
        assert(len(atoms) == 3)
        assert_equal(pos, atoms.q)


def test_iter_xyz():
    with open(local("test.pos_0.xyz"), "r") as f:
        for num, atoms in enumerate(io_xyz.iter_xyz(f)):
            assert(len(atoms) == 3)
            assert_equal(pos*(num+1), atoms.q)


def test_read_pdb():
    with open(local("test.pos_0.pdb"), "r") as f:
        atoms, cell = io_pdb.read_pdb(f)
        assert(len(atoms) == 3)
        assert_equal(pos, atoms.q)
        # TODO: test cell


def test_iter_pdb():
    with open(local("test.pos_0.pdb"), "r") as f:
        for num, (atoms, cell) in enumerate(io_pdb.iter_pdb(f)):
            assert(len(atoms) == 3)
            assert_equal(pos*(num+1), atoms.q)


def test_print_pdb():
    with open(local("test.pos_0.pdb"), "r") as f:
        with open(local("test.pos_1.xyz"), "w") as out:
            for num, (atoms, cell) in enumerate(io_pdb.iter_pdb(f)):
                assert(len(atoms) == 3)
                assert_equal(pos*(num+1), atoms.q)
                io_xyz.print_xyz(atoms, Cell(h=np.identity(3, float)), filedesc=out)
    assert(filecmp.cmp(local("test.pos_0.xyz"), local("test.pos_1.xyz")))
    os.unlink(local("test.pos_1.xyz"))


def test_print_xyz():
    with open(local("test.pos_0.pdb"), "r") as f:
        with open(local("test.pos_1.pdb"), "w") as out:
            for num, (atoms, cell) in enumerate(io_pdb.iter_pdb(f)):
                assert(len(atoms) == 3)
                assert_equal(pos*(num+1), atoms.q)
                io_pdb.print_pdb(atoms, Cell(h=np.identity(3, float)), filedesc=out)
    assert(filecmp.cmp(local("test.pos_0.pdb"), local("test.pos_1.pdb")))
    os.unlink(local("test.pos_1.pdb"))
