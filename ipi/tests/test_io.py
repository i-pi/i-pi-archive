"""Deals with testing the io system.

Note that this will only run if you have Python version 2.5 or later.
Otherwise, replace all the with statements with f = filestream.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import filecmp
import os
import sys

import numpy as np
from numpy.testing import assert_equal

from common import local
from ipi.engine.cell import Cell
from ipi.utils.io import io_xyz
from ipi.utils.io import io_pdb
import ipi.utils.io as io


pos = np.array([i for i in range(3*3)])


def test_read_xyz():
   """Tests that xyz files are read correctly."""

   with open(local("test.pos_0.xyz"), "r") as f:
      atoms = io_xyz.read_xyz(f)
      assert(len(atoms) == 3)
      assert_equal(pos, atoms.q)

def test_iter_xyz():
   """Tests that xyz files with multiple frames are read correctly."""

   with open(local("test.pos_0.xyz"), "r") as f:
      for num, atoms in enumerate(io_xyz.iter_xyz(f)):
         assert(len(atoms) == 3)
         assert_equal(pos*(num+1), atoms.q)

def test_read_pdb():
   """Tests that pdb files are read correctly."""

   with open(local("test.pos_0.pdb"), "r") as f:
      atoms, cell = io_pdb.read_pdb(f)
      assert(len(atoms) == 3)
      assert_equal(pos, atoms.q)
      # TODO: test cell

def test_iter_pdb():
   """Tests that pdb files with multiple frames are read correctly."""

   with open(local("test.pos_0.pdb"), "r") as f:
      for num, (atoms, cell) in enumerate(io_pdb.iter_pdb(f)):
         assert(len(atoms) == 3)
         assert_equal(pos*(num+1), atoms.q)

def test_print_pdb():
   """Tests that pdb files are printed correctly."""

   with open(local("test.pos_0.pdb"), "r") as f:
      with open(local("test.pos_1.xyz"), "w") as out:
         for num, (atoms, cell) in enumerate(io_pdb.iter_pdb(f)):
            assert(len(atoms) == 3)
            assert_equal(pos*(num+1), atoms.q)
            io_xyz.print_xyz(atoms, Cell(h=np.identity(3, float)), filedesc=out)

   assert(filecmp.cmp(local("test.pos_0.xyz"), local("test.pos_1.xyz")))
   os.unlink(local("test.pos_1.xyz"))

def test_print_xyz():
   """Tests that xyz files are printed correctly."""

   with open(local("test.pos_0.pdb"), "r") as f:
      with open(local("test.pos_1.pdb"), "w") as out:
         for num, (atoms, cell) in enumerate(io_pdb.iter_pdb(f)):
            assert(len(atoms) == 3)
            assert_equal(pos*(num+1), atoms.q)
            io_pdb.print_pdb(atoms, Cell(h=np.identity(3, float)), filedesc=out)

   assert(filecmp.cmp(local("test.pos_0.pdb"), local("test.pos_1.pdb")))
   os.unlink(local("test.pos_1.pdb"))

def test_read_xyz2():
   """Tests that mode/xyz files are read correctly."""

   with open(local("test.pos_0.xyz"), "r") as f:
      atoms = io.read_file("xyz", f)
      assert(len(atoms) == 3)
      assert_equal(pos, atoms.q)

def test_iter_xyz2():
   """Tests that mode/xyz files with multiple frames are read correctly."""

   with open(local("test.pos_0.xyz"), "r") as f:
      for num, atoms in enumerate(io.iter_file("xyz", f)):
         assert(len(atoms) == 3)
         assert_equal(pos*(num+1), atoms.q)

def test_read_pdb2():
   """Tests that mode/pdb files are read correctly."""

   with open(local("test.pos_0.pdb"), "r") as f:
      atoms, cell = io.read_file("pdb", f)
      assert(len(atoms) == 3)
      assert_equal(pos, atoms.q)
      # TODO: test cell

def test_iter_pdb2():
   """Tests that mode/pdb files with multiple frames are read correctly."""

   with open(local("test.pos_0.pdb"), "r") as f:
      for num, (atoms, cell) in enumerate(io.iter_file("pdb", f)):
         assert(len(atoms) == 3)
         assert_equal(pos*(num+1), atoms.q)

def test_print_pdb2():
   """Tests that mode/pdb files are printed correctly."""

   with open(local("test.pos_0.pdb"), "r") as f:
      with open(local("test.pos_1.xyz"), "w") as out:
         for num, (atoms, cell) in enumerate(io.iter_file("pdb", f)):
            assert(len(atoms) == 3)
            assert_equal(pos*(num+1), atoms.q)
            io.print_file("xyz", atoms, Cell(h=np.identity(3, float)), filedesc=out)

   assert(filecmp.cmp(local("test.pos_0.xyz"), local("test.pos_1.xyz")))
   os.unlink(local("test.pos_1.xyz"))

def test_print_xyz2():
   """Tests that mode/xyz files are printed correctly."""

   with open(local("test.pos_0.pdb"), "r") as f:
      with open(local("test.pos_1.pdb"), "w") as out:
         for num, (atoms, cell) in enumerate(io.iter_file("pdb", f)):
            assert(len(atoms) == 3)
            assert_equal(pos*(num+1), atoms.q)
            io.print_file("pdb", atoms, Cell(h=np.identity(3, float)), filedesc=out)

   assert(filecmp.cmp(local("test.pos_0.pdb"), local("test.pos_1.pdb")))
   os.unlink(local("test.pos_1.pdb"))
