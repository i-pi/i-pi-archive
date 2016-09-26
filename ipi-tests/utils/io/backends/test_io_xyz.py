#!/usr/bin/env python2

import numpy as np
from ipi.utils.units import Elements


def generate_rand_xyz(natoms):
    all_elem = Elements.mass_list.keys()
    xyz = np.random.random((natoms, 3))
    xyz = (xyz * 10 - 5) # To have both, positive and negative numbers
    names = [all_elem[i] for i in
             np.random.randint(0, len(all_elem), natoms)]
    return (xyz, names)


def generate_xyz_traj(natoms, nframe):
    pass
