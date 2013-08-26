"""Deals with testing the io system.

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


Note that this will only run if you have Python version 2.5 or later.
Otherwise, replace all the with statements with f = filestream.
"""

import sys
sys.path.append("../")
sys.path.append("../../")

import filecmp
import os, sys
import numpy as np
from numpy.testing import assert_equal
from common import local

from ipi.engine.cell import Cell

from ipi.utils.io import io_xyz
from ipi.utils.io import io_pdb

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
