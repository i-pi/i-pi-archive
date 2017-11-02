#-*- coding: utf-8 -*-

import numpy as np
from scipy import linalg
try:
	import lapack2py
except:
	print 'WARNING linalg.py: unable to import lapack2py'
import warnings

def band(a, u, lower=False):
	"""Return matrix a in banded form where u in the number of bands above the diagonal."""
	M = a.shape[0]
#	a_band = np.empty((u+1,M))
	a_band = np.zeros((u+1,M))
	if lower: # lower diagonal ordered form
		a_band[0] = np.diag(a)
		for i in xrange(1,u+1):
			a_band[i,:-i] = np.diag(a,-i)
	else: # upper diagonal ordered form
		a_band[-1] = np.diag(a)
		for i in xrange(1,u+1):
			a_band[-i-1,i:] = np.diag(a,i)
	return a_band

def unband(a_band, lower=False, symmetric=True):
	"""Return matrix in standard form."""
	M = a_band.shape[1]
	u = a_band.shape[0] - 1
	if lower: # lower diagonal ordered form
		a = np.diag(a_band[0])
		for i in xrange(1,u+1):
			a += np.diag(a_band[i,:-i],-i)
	else: # upper diagonal ordered form
		a = np.diag(a_band[-1])
		for i in xrange(1,u+1):
			a += np.diag(a_band[-i-1,i:],i)
	if symmetric:
		a += a.T - np.diag(np.diag(a))
	return a

def band_from_block(A, lower=False):
	"""Return matrix in upper diagonal ordered form
	when presented with a set of square blocks, A with shape (N,d,d)"""
	assert lower == False # this could be generalized to do lower ordered form too, but hasn't yet
	N = len(A)
	d = A.shape[1]
	assert A.shape == (N,d,d)
	bm = np.zeros((d, N, d))
	for i in xrange(N):
		for j in xrange(d):
			bm[-1-j,i,j:] += np.diag(A[i], j)
	bm.shape = (d, N*d)
	return bm

def lu_factor_banded(a_band):
	# assume A was a square matrix with l==u
	l = u = (a_band.shape[0] - 1)//2
	M = N = a_band.shape[1]
	AB = np.empty((2*l+u+1,N), order='F')
	AB[l:] = a_band
	ipiv, info = lapack2py.wrap_dgbtrf(M,l,u,AB)
	return AB, ipiv

def lu_solve_banded(lu_and_piv, b):
	AB, ipiv = lu_and_piv
	trans = 'N' # no transpose
	LDAB, N = AB.shape
	l = u = (LDAB - 1)//3
	B = np.array(b, order='F')
	if B.ndim == 1:
		B.shape = (-1,1)
	LDB, NRHS = B.shape
	info = lapack2py.wrap_dgbtrs(trans,l,u,AB,ipiv,B)
	return B

UPLO = {False: 'U', True: 'L'}

def geneigvals_banded(a_band, b_band, lower=False):
	"""Solve generalized eigenvalue problem, A*x = lambda*B*x, for banded matrices
	and return eigenvalues lambda"""
	A = np.array(a_band, order='F')
	B = np.array(b_band, order='F')
	info = lapack2py.wrap_dpbstf(B, UPLO[lower])
	if info > 0:
		raise linalg.LinAlgError('ERROR from DPBSTF: the factorization could not be completed because the update element a(%i,%i) was negative; the matrix A is not positive definate.'%(info,info))
	info = lapack2py.wrap_dsbgst(A, B, UPLO[lower])
	return linalg.eigvals_banded(A, lower=lower)

def orthogonalize(A):
	"""Return orthogonal matrix Q where Q[:,0] is proportional to A[:,0] and has the same sign
	using QR decomposition: A = Q.R where Q is orthogonal and R is upper triangular."""
	Q, R = linalg.qr(A)
	if R[0,0] < 0:
		Q *= -1
		# R *= -1
	return Q

#def expm(A):
#	return linalg.expm(A)

def dexpm(A, dA):
	"""Return derivative of expm(A) = int_0^1 expm(αA) dA expm((1-α)A) dα
	where A is a symmetric matrix and dA is the derivative of A with respect to some variable
	"""
	W, U = linalg.eigh(A)
	tmp = np.dot(U.T, np.dot(dA, U))
	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		Kubo = np.where(W.reshape(-1,1)==W, np.diag(np.exp(W)), (np.exp(W.reshape(-1,1)) - np.exp(W)) / (W.reshape(-1,1) - W))
	return np.dot(U, np.dot(tmp*Kubo, U.T))

def dotdot(A, B, C):
	return np.dot(A, np.dot(B, C))

def invmul(A, B):
	"""Return x = dot(inv(A),B), equivalent to solving for A.x==B"""
	return linalg.solve(A, B)

def sym_band(A):
	"""Return symmetric banded matrix from just upper banded."""
	u = len(A) - 1
	l = u
	M = A.shape[1]
	newA = np.zeros((u+l+1,M))
	newA[:u+1] = A
	for i in xrange(1,l+1):
		newA[u+i,:M-i] = A[-1-i,i:]
	return newA

def invmul_banded(A, B, posdef=False):
	"""A is in upper banded form"""
	if posdef:
		return linalg.solveh_banded(A, B)
	else:
		u = len(A) - 1
		l = u
		newA = sym_band(A)
		return linalg.solve_banded((l,u), newA, B)

########################

if __name__ == "__main__":
	a = np.random.normal(0,1,(3,3))
	a += a.T
	a[0,-1] = a[-1,0] = 0
	print a
	b = np.random.normal(0,1,(3,3)) + 3*np.identity(3)
	b += b.T
	b[0,-1] = b[-1,0] = 0
	print b
	lower = False
	eigvals_only = False
	print "scipy.linalg.eigh gives ", linalg.eigh(a,b,eigvals_only=eigvals_only)

	a_band = band(a,1,lower)
	b_band = band(b,1,lower)
	print "geneigvals_banded gives ", geneigvals_banded(a_band,b_band,lower)

	print a_band
	print sym_band(a_band)

	A = np.random.randint(1, 10, size=(4,3,3))
	print A
	print band_from_block(A)

	print '\nexpm'

	import numderiv
	beta = 1.4
	x = np.random.normal(0,1)
	def pot(x): return np.array([[(x-1)**2/2, 0.1], [0.1, (x+1)**2/2]])
	def grad(x): return np.array([[x-1, 0], [0, x+1]])
	V = pot(x)
	print V
	print linalg.expm(-beta*V)
	print numderiv.first2(x, lambda x: linalg.expm(-beta*pot(x)), 1e-3)
	print dexpm(-beta*V,-beta*grad(x))
