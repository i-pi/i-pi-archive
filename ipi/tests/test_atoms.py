"""Deals with testing the Atoms object.

Copyright (C) 2014, i-PI Development Team

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http.//www.gnu.org/licenses/>.
"""

from common import local

from ipi.utils.io import io_xyz

def get_atoms(fin):
   """Reads atoms object from file @fin."""

   with open(local(fin), "r") as f:
      atoms = io_xyz.read_xyz(f)
   return atoms

def test_names():
   """Tests names of Atoms object."""
   atoms = get_atoms("test.pos_0.xyz")
   expected = ["O", "H", "H"]
   assert(len(atoms.names) == 3)
   for i, name in enumerate(atoms.names):
      assert(name == expected[i])
      assert(name == atoms[i].name)

   # Same test with iterator instead
   for i, atom in enumerate(atoms):
      assert(atom.name == expected[i])
      assert(atom.name == atoms.names[i])
