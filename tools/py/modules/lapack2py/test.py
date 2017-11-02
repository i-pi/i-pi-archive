#!/usr/bin/python

import lapack2py
from scipy import linalg
import numpy as np
from linalg import band

lower = True
eigvals_only = False

a = np.random.normal(0,1,(3,3))
a += a.T
a[0,-1] = a[-1,0] = 0
print a
a_band = band(a,1,lower)
print a_band

b = np.random.normal(0,1,(3,3)) + 3*np.identity(3)
b += b.T
b[0,-1] = b[-1,0] = 0
#b = np.identity(3)
print b
b_band = band(b,1,lower)
print b_band

print "scipy.linalg.eigh gives ", linalg.eigh(a,b,eigvals_only=eigvals_only)

if lower:
	UPLO = 'L'
else:
	UPLO = 'U'

A = np.array(a_band, order='F')
B = np.array(b_band, order='F')
print A
print B

info = lapack2py.wrap_dpbstf(B, UPLO)
print info
print B
info = lapack2py.wrap_dsbgst(A, B, UPLO)
print info
print A
b = linalg.eig_banded(A, lower=lower, eigvals_only=eigvals_only)
print "lapack2py.geneigvals_banded gives ", b
