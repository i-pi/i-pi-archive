#!/usr/bin/env python2
# pylint: disable=C0111,W0621,R0914,C0301
#+easier to find important problems

import re
import pytest
import numpy as np
import numpy.testing as npt

import ipi_tests.utils.io.backends.io_xyz_utils as xyz_gen
import ipi.utils.io.backends.io_xyz as io_xyz
import ipi.utils.mathtools as mt


#######################
# Testing reading xyz #
#######################

deg2rad = np.pi/180.0
cell_string = ' '.join([str(x) for x in mt.abc2h(5.1, 5.2, 5.0,
                                                 91*deg2rad,
                                                 89*deg2rad,
                                                 90*deg2rad).flatten()])
default_cell_mat = mt.abc2h(1.0, 1.0, 1.0, np.pi/2.0, np.pi/2.0, np.pi/2.0)

# natoms, frames, comment, expected_cell, precision
test_xyz = [
    (1, 1, 'just a string plus few numbers: 1.10 2 .1',
     default_cell_mat, 10),
    (2, 1, 'another random comment',
     default_cell_mat, 10),
    (1, 2, 'random comment', default_cell_mat, 10),
    (2, 2, 'random comment', default_cell_mat, 10),
    (10, 2, 'random comment', default_cell_mat, 10),
    (10, 10, 'random comment', default_cell_mat, 10),
    (2, 3, '# CELL(abcABC): 5.1 5.2 5.0 91.0  89  90',
     mt.abc2h(5.1, 5.2, 5.0, 91*deg2rad, 89*deg2rad, 90*deg2rad), 10),
    (2, 3, '# CELL[GENH]: '+ cell_string,
     mt.abc2h(5.1, 5.2, 5.0, 91*deg2rad, 89*deg2rad, 90*deg2rad), 10),
    (2, 3, '# CELL{H}: '+ cell_string,
     mt.abc2h(*mt.genh2abc(mt.abc2h(5.1, 5.2, 5.0, 91*deg2rad, 89*deg2rad, 90*deg2rad))), 10),
]


@pytest.fixture(params=test_xyz)
def create_random_xyz_traj(request):

    natoms, frames, comment, expected_cell, precision = request.param

    filedesc, xyz, atom_names = xyz_gen.xyz_traj_filedesc(natoms,
                                                          frames, comment)
    filedesc.seek(0)

    check_comment = re.match('# CELL[\(\[\{]H[\)\]\}]: ([-0-9\.?Ee ]*)\s*', comment) #pylint: disable=W1401

    if check_comment:
        genh = np.array(check_comment.group(1).split()[:9], float).reshape((3, 3))
        invgenh = np.linalg.inv(genh)
        for _ui in xrange(natoms*frames):
            _uu = np.array([xyz[3*_ui], xyz[3*_ui+1], xyz[3*_ui+2]])
            _us = np.dot(_uu, invgenh)
            _uu = np.dot(expected_cell, _us)
            xyz[3*_ui], xyz[3*_ui+1], xyz[3*_ui+2] = _uu


    return (filedesc, xyz, atom_names, natoms, frames,
            comment, expected_cell, precision)


def test_read_xyz(create_random_xyz_traj):

    filedesc, xyz, atom_names, \
        natoms, frames, comment, expected_cell, precision = create_random_xyz_traj

    for _fr in xrange(frames):

        tcomment, tcell, tqatoms, tnames, tmasses = io_xyz.read_xyz(filedesc)

        assert tcomment.strip() == comment
        npt.assert_array_almost_equal(tqatoms, xyz[_fr*natoms*3:_fr*3*natoms+3*natoms], decimal=precision)
        npt.assert_array_equal(np.array(atom_names[_fr*natoms:_fr*natoms+natoms], dtype='|S4'), tnames)
        npt.assert_array_almost_equal(tcell, expected_cell, decimal=precision)

def test_iter_xyz(create_random_xyz_traj):
    filedesc, xyz, atom_names, \
        natoms, junk, comment, expected_cell, precision = create_random_xyz_traj

    _fr = 0
    for _io in io_xyz.iter_xyz(filedesc):
        tcomment, tcell, tqatoms, tnames, tmasses = _io
        assert tcomment.strip() == comment
        npt.assert_array_almost_equal(tqatoms, xyz[_fr*natoms*3:_fr*3*natoms+3*natoms], decimal=precision)
        npt.assert_array_equal(np.array(atom_names[_fr*natoms:_fr*natoms+natoms], dtype='|S4'), tnames)
        npt.assert_array_almost_equal(tcell, expected_cell, decimal=precision)
        _fr += 1


# if __name__ == '__main__':

#     for st in test_cell:
#         string, float_n = st
#         floats = 100.0 * np.random.random_sample(float_n)

#         # Generate also a random number of spaces between numbers
#         spaces = ' ' * np.random.random_integers(1, 12)
#         string += spaces
#         string += spaces.join([str(x) for x in floats.tolist()])
#         print string
