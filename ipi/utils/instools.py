import numpy as np
from scipy import linalg


def banded_hessian(h,im,shift=0.001):
    """Given Hessian in the reduced format (h), construct
    the upper band hessian including the RP terms"""
    nbeads   = im.dbeads.nbeads
    natoms   = im.dbeads.natoms
    ii       = natoms * 3 * nbeads
    ndiag    = natoms * 3 + 1 #only upper diagonal form

    #np.set_printoptions(precision=6, suppress=True, threshold=np.nan, linewidth=1000)

    hnew = np.zeros((ndiag, ii))

    #add physical part
    for i in range(nbeads):
        h_aux = h[:,i*natoms*3:(i+1)*natoms*3] #Peaks one physical hessian
        for j in range(1,ndiag):
            hnew[j,(ndiag-1-j)+i*natoms*3:(i+1)*natoms*3] = np.diag(h_aux,ndiag-1-j)

    #add spring parts
    # Diagonal
    d_corner = im.dbeads.m3[0] * im.omega2
    d_0 = np.array([[d_corner * 2]]).repeat(im.dbeads.nbeads - 2, axis=0).flatten()
    diag_sp = np.concatenate((d_corner, d_0, d_corner))
    hnew[-1, :] += diag_sp

    # Non-Diagonal
    d_out = - d_corner
    ndiag_sp = np.array([[d_out]]).repeat(im.dbeads.nbeads-1, axis=0).flatten()
    hnew[0, :] = np.concatenate((np.zeros(natoms*3),ndiag_sp ))

    # Add safety shift value
    hnew[-1, :] += shift

    return hnew

def sym_band(A):
    """Return symmetric banded matrix from just upper banded."""
    u = len(A) - 1
    l = u
    M = A.shape[1]
    newA = np.empty((u+l+1,M))
    newA[:u+1] = A
    for i in xrange(1,l+1):
        newA[u+i,:M-i] = A[-1-i,i:]
    return newA

def invmul_banded(A, B, posdef=False):
    """A is in upper banded form
	Solve H.h = -G for Newton - Raphson step, h
    using invmul_banded(H, -G) take step x += h
    to  find minimum or transition state """

    if posdef:
        return linalg.solveh_banded(A, B)
    else:
        u = len(A) - 1
        l = u
        newA = sym_band(A)
        #np.set_printoptions(precision=6, suppress=True, threshold=np.nan, linewidth=1000)
        #print linalg.eigvals_banded(A)
        #sys.exit(0)
        return linalg.solve_banded((l,u), newA, B)