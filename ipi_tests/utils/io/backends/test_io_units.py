#!/usr/bin/env python2

import pytest
import ipi_tests.utils.io.backends.io_xyz_utils as xyz_gen
import ipi.utils.mathtools as mt
import ipi.utils.io.backends.io_units as testing

from ipi.utils.units import Elements
import numpy as np
import numpy.testing as npt

# @pytest.fixture()
# def create_random_xyz_traj_to_read(request):

#     natoms, frames, comment, expected_cell, precision = request.param

#     filedesc, xyz, atom_names = xyz_gen.xyz_traj_filedesc(natoms,
#                                                           frames, comment)
#     filedesc.seek(0)

def test_process_units_noobj():

    natoms, frames, comment, output = 1, 1, 'nocomment', 'asd'

    filedesc, xyz, atom_names = xyz_gen.xyz_traj_filedesc(natoms,
                                                          frames, comment)

    cell = mt.abc2h(1.0, 1.0, 1.0, np.pi/2.0, np.pi/2.0, np.pi/2.0)

    masses = []
    _count = 0
    print 'atom_names', atom_names
    for _at in atom_names:
        print Elements.mass(_at), _count, _at
        masses.append(Elements.mass(_at))
        _count += 1

    masses = np.array(masses)
    res = testing.process_units(comment, cell, xyz, atom_names, masses, output=output)

    npt.assert_array_almost_equal(res['data'], xyz, 5)
    npt.assert_array_almost_equal(res['masses'], masses, 5)
    npt.assert_array_almost_equal(res['cell'], cell, 5)
    npt.assert_array_equal(res['names'], atom_names, 5)
    assert res['natoms'] == natoms
