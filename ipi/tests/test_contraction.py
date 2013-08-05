from ipi.utils import nmtransform
import numpy as np
from numpy.testing import assert_almost_equal as assert_equals


def get_identicals(q):
    """Return bit mask, which particles coincide in all beads.
    """
    centroid = np.sum(q, axis=0)/float(q.shape[0])
    ret = []
    for i, c in enumerate(centroid):
        for bead in q:
            if abs(c-bead[i]) > 1e-5:
                ret.append(0)
                break
        else:
            ret.append(1)
    return ret


def check_up_and_down_scaling_pos(n, q):
    """Check if q can be upscaled to n beads and downscaled again.
    """
    rescale = nmtransform.mk_rs_matrix(q.shape[0], n)
    msk = get_identicals(q)
    print "Initial position of the bead system:"
    print q, q.shape, msk, rescale.shape

    # rescale up to the n beads
    beads_n = np.dot(rescale, q)#*n/q.shape[0]
    print "Upscaled to %d beads:"%n
    print beads_n, beads_n.shape

    rescale = nmtransform.mk_rs_matrix(n, q.shape[0])
    beads_final = np.dot(rescale, beads_n)
    print "Final position of the bead system:"
    print beads_final

    for bead in beads_n:
        for i, m in enumerate(msk):
            if m:
                assert_equals(bead[i], q[0][i])
    assert_equals(q, beads_final)


def check_up_and_down_scaling_frc(n, q):
    """Check if q can be upscaled to n beads and downscaled again.
    """
    rescale = nmtransform.nm_rescale(n, q.shape[0])
    msk = get_identicals(q)
    print "Initial position of the bead system:"
    print q, q.shape, msk

    # rescale up to the n beads
    beads_n = rescale.b2tob1(q)
    print "Upscaled to %d beads:"%n
    print beads_n, beads_n.shape

    beads_final = np.dot(rescale._b1tob2, beads_n)
    print "Final position of the bead system:"
    print beads_final

    for bead in beads_n:
        for i, m in enumerate(msk):
            if m:
                assert_equals(bead[i], q[0][i])
    assert_equals(q, beads_final)


numbers_to_check = range(10, 56, 9)


def test_1_to_n():
    for n in numbers_to_check:
        q = np.array([[0.0,0.0,0.0, 1.0,0.0,0.0]])
        yield check_up_and_down_scaling_pos, n, q
        yield check_up_and_down_scaling_frc, n, q


def test_2_to_n():
    for n in numbers_to_check:
        q = np.array([[0.0,0.0,0.0, 1.0,0.0,0.0],
                      [0.0,0.1,0.0, 1.0,0.1,0.0]])
        yield check_up_and_down_scaling_pos, n, q
        yield check_up_and_down_scaling_frc, n, q


def test_3_to_n():
    for n in numbers_to_check:
        q = np.array([[0.0, 0.0,0.0, 1.0, 0.0,0.0],
                      [0.0, 0.1,0.0, 1.0, 0.1,0.0],
                      [0.0,-0.1,0.0, 1.0,-0.1,0.0]])
        yield check_up_and_down_scaling_pos, n, q
        yield check_up_and_down_scaling_frc, n, q


def test_4_to_n():
    for n in numbers_to_check:
        q = np.array([[0.0, 0.0,0.0, 1.0, 0.0,0.0],
                      [0.0, 0.1,0.0, 1.0, 0.1,0.0],
                      [0.0, 0.2,0.0, 1.0, 0.2,0.0],
                      [0.0,-0.1,0.0, 1.0,-0.1,0.0]])
        yield check_up_and_down_scaling_pos, n, q
        yield check_up_and_down_scaling_frc, n, q
