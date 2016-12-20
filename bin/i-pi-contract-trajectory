#!/usr/bin/env python2

description = """
Read individual beads from a set of trajectory files (positions or other data)
and contract them to a different number of beads using ring polymer
contraction. Trajectory file formats are inferred from file extensions, the
number of input beads is given by the number of input files and the number of
output beads is given as an argument. Cell dimensions are taken from the first
input file.

As a special case, this can calculate the centroid.
"""

import os
import sys
import argparse

import numpy as np

from ipi.utils.nmtransform import nm_rescale
from ipi.utils.io import open_backup, iter_file_name, print_file


def contract_trajectory(fns_in, fn_out_template, n_new):

    n = len(fns_in)

    # Generate output file names.
    if n_new == 1:
        fns_out = [fn_out_template]
    else:
        fns_out = [fn_out_template.format(i) for i in range(n_new)]

    print "Contracting {:d} beads to {:d} beads.".format(n, n_new)
    print

    print "input file names:"
    for fn in fns_in:
        print fn
    print

    print "output file names:"
    for fn in fns_out:
        print fn
    print

    # Open input trajectory iterators.
    trjs_in = [iter_file_name(fn) for fn in fns_in]

    # Open output files.
    fs_out = [open_backup(fn, "w") for fn in fns_out]
    mode_out = os.path.splitext(fn_out_template)[1]

    # prepare ring polymer rescaler
    rescale = nm_rescale(n, n_new)

    # Loop over all frames.
    i_frame = 0
    while True:

        try:
            # Get the frames for all beads.
            frames = [trj.next() for trj in trjs_in]
        except StopIteration:
            # Stop when any of the trajectories runs out of frames.
            break

        # Consistency check.
        h = frames[0]["cell"].h
        natoms = frames[0]["atoms"].natoms
        for i in range(n):

            # Check that all the cells are the same.
            if (frames[i]["cell"].h != h).any():
                msg = "Cell for beads {:d} and {:d} differ in frame {:d}."
                raise ValueError(msg.format(0, i, i_frame))

            # Check that the numbers of atoms are the same.
            if frames[i]["atoms"].natoms != natoms:
                msg = "Different numbers fo atoms for beads {:d} and {:d} in frame {:d}."
                raise ValueError(msg.format(0, i, i_frame))

        # Reuse the first frame for output.
        cell = frames[0]["cell"]
        atoms = frames[0]["atoms"]

        # Compose the ring polymer.
        q = np.vstack([frame["atoms"].q for frame in frames])

        # Contract the coordinates to `n_new` beads.
        q_c = rescale.b1tob2(q)

        # Save the output data.
        for i, f_out in enumerate(fs_out):
            atoms.q = q_c[i,:]
            print_file(mode_out, atoms, cell, f_out)

        # Count frames and print information on progress.
        i_frame += 1
        if i_frame % 100 == 0:
            print "\rframe {:d}".format(i_frame),
        sys.stdout.flush()

    for f_out in fs_out:
        f_out.close()

    print
    print
    print "Processed {:d} frames.".format(i_frame)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("filenames", type=str, nargs="+",
                        help="Bead trajectories to process.")
    parser.add_argument("--filename-out", type=str,
                        help="Template for output trajectory file names. "
                             "Must contain one formatting field for an integer if n>1. "
                             "Example: 'positions-contracted-bead-{:02d}.xyz'")
    parser.add_argument("--n", type=int,
                        help="Number of beads to contract to.")

    args = parser.parse_args()

    # Process everything.
    contract_trajectory(args.filenames, args.filename_out, args.n)
