#!/usr/bin/env python2

import pytest
import numpy as np
import numpy.testing as npt

import ipi_tests.utils.io.backends.io_xyz_utils as xyz_gen
import ipi.utils.io.backends.io_xyz as io_xyz


test_data = [
    (10, 1, 'positions{angstrom}', 10)
]

@pytest.mark.parametrize('natoms,frames,comment,precision', test_data)
def test_read_xyz(natoms, frames, comment, precision):
    filedesc, xyz, atom_names = xyz_gen.xyz_traj_filedesc(natoms, frames, comment)

    filedesc.seek(0)
    tcomment, tcell, tqatoms, tnames, tmasses = io_xyz.read_xyz(filedesc)
    assert tcomment.strip() == comment
    npt.assert_array_almost_equal(tqatoms, xyz, decimal=precision)
    npt.assert_array_equal( np.array(atom_names, dtype='|S4'), tnames)
