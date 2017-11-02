#!/usr/bin/python

"""For imaginary-time trajectories only"""

import numpy as np
import linalg as mylinalg
from linalg import invmul
from scipy import linalg, integrate
import pathlength

if 1: # LU
	banded = True
	if banded:
		factor = mylinalg.lu_factor_banded
		solve = mylinalg.lu_solve_banded
	else:
		factor = linalg.lu_factor
		solve = linalg.lu_solve
else: # CHOLESKY (also assumes that it is positive definite, which it isn't!)
	banded = True
	if banded:
		factor = lambda x: (linalg.cholesky_banded(x), False)
		solve = linalg.cho_solve_banded
	else:
		factor = linalg.cho_factor
		solve = linalg.cho_solve

def Jmatrix(d2Vdx2, dt, m=1.):
	N = len(d2Vdx2) - 1
	f = d2Vdx2.shape[1]
	t = np.sum(dt)
	eps = dt / t
	J = np.zeros((N-1,f,N-1,f))
	for n in xrange(f):
		J[:,n,:,n] = np.identity(N-1) * (1./eps[:-1] + 1./eps[1:])
		J[:,n,:,n] -= np.diag(np.ones(N-2),1) / eps[:-1]
		J[:,n,:,n] -= np.diag(np.ones(N-2),-1) / eps[1:]
	for i in xrange(N-1):
		J[i,:,i,:] += (eps[i]+eps[i+1])*t**2/(2.*m) * d2Vdx2[i+1]
#	J += (np.diag((eps[:-1]+eps[1:])*t**2/(2.*m)) * d2Vdx2[1:-1]).reshape(N-1,f,1,f) * np.identity(N-1).reshape(1,1,N-1,1)
	J.shape = ((N-1)*f,(N-1)*f)
	return J / N

def Jmatrix_banded(d2Vdx2, dt, m=1.):
	"""Return J in upper diagonal form.
	N.B. this can only be done if J is symmetric which is only true if all dts are the same"""
	assert np.all(abs(dt-np.mean(dt))/np.mean(dt)<1e-8) # ensure all dts are the same
	N = len(d2Vdx2) - 1
	f = d2Vdx2.shape[1]
	t = np.sum(dt)
	eps = dt / t
	J = np.zeros((f+1,(N-1)*f))
	J[1:] = mylinalg.band_from_block((eps[:-1]+eps[1:]).reshape(-1,1,1)*t**2/(2.*m) * d2Vdx2[1:-1])
	J.shape = (f+1,N-1,f)
	J[0] -= (1./eps[1:]).reshape(-1,1)
	J[-1] += (1./eps[:-1] + 1./eps[1:]).reshape(-1,1)
	J.shape = (f+1,(N-1)*f)
	return J / N

class Trajectory(object):
	def __init__(self, x, V, dVdx, d2Vdx2, dt, m=1., linear=False, massweight=False):
		"""
		x.shape = (N+1,f)
		V.shape = (N+1,)
		dVdx.shape = (N+1,f)
		d2Vdx2.shape = (N+1,f,f)
		dt.shape = (N,) or float
		Compute derivatives w.r.t. xi=x[0] and xf=x[-1]

		If massweight is False, m must be a scalar
		However, if massweight is True, m can be an array of shape (f,)
		and x, dVdx and d2Vdx2 are mass weighted automatically;
		in this case all outputs (e.g. dS/dx) are also mass-weighted.
		"""

		self.N = len(x) - 1
		self.f = x.shape[1]
		self.x = x
		self.V = np.asarray(V)
		self.dVdx = np.asarray(dVdx)
		if np.asarray(d2Vdx2).ndim == 2:
			self.d2Vdx2 = np.array([d2Vdx2]*(self.N+1))
		elif np.asarray(d2Vdx2).ndim == 3:
			self.d2Vdx2 = np.asarray(d2Vdx2)
		if np.asarray(dt).ndim == 0:
			self.dt = np.ones(self.N) * dt
		else:
			self.dt = np.asarray(dt)
		self.t = np.sum(self.dt)
		self.eps = self.dt / float(self.t)
		assert self.x.shape == (self.N+1,self.f)
		assert self.V.shape == (self.N+1,)
		assert self.dVdx.shape == (self.N+1,self.f)
		assert self.d2Vdx2.shape == (self.N+1,self.f,self.f)
		assert self.dt.shape == (self.N,)
		if massweight:
			m = np.atleast_1d(m)
			assert m.shape == (self.f,)
			self.x *= np.sqrt(m)
			self.dVdx /= np.sqrt(m)
			self.d2Vdx2 /= np.sqrt(np.outer(m,m))
			self.m = 1
		else:
			self.m = m
			assert not hasattr(self.m, '__getitem__') # ensure that m is scalar
		if banded: # lower=False
			#J_band = mylinalg.band(self.J,self.f)
			J_band = Jmatrix_banded(self.d2Vdx2, self.dt, self.m)
			self.JF = factor(mylinalg.sym_band(J_band))
		else:
			self.J = Jmatrix(self.d2Vdx2, self.dt, self.m)
			self.JF = factor(self.J)
		self.linear = linear
		if linear:
			self.dx = self.x[1:] - self.x[:-1]
			self.adx = np.sqrt(np.sum(self.dx**2,1))
			self.kappa = (self.V[1:] - self.V[:-1]) / self.adx

	def delta(self, i, j):
		d = np.zeros((self.N-1,self.f))
		d[i,j] = 1
		return d.ravel()
	
	def S(self):
		S = self.m*np.sum((self.x[1:]-self.x[:-1])**2/self.dt.reshape(-1,1))/2 + np.sum(self.dt*(self.V[1:] + self.V[:-1]))/2
		if self.linear:
			S -= np.sum(self.kappa**2 * self.dt**3) / (24. * self.m)
		return S

	def dSdx(self):
		dSdxi = self.dt[0]/2*self.dVdx[0] + self.m/self.dt[0]*(self.x[0] - self.x[1])
		dSdxf = self.dt[-1]/2*self.dVdx[-1] + self.m/self.dt[-1]*(self.x[-1] - self.x[-2])
		if self.linear:
			dSdxi -= self.kappa[0]*self.dt[0]**3/(12.*self.m) * (self.adx[0]*-self.dVdx[0] - (self.V[1]-self.V[0])*(self.x[0]-self.x[1])/self.adx[0])/self.adx[0]**2
			dSdxf -= self.kappa[-1]*self.dt[-1]**3/(12.*self.m) * (self.adx[-1]*self.dVdx[-1] - (self.V[-1]-self.V[-2])*(self.x[-1]-self.x[-2])/self.adx[-1])/self.adx[-1]**2
		return np.concatenate([dSdxi, dSdxf])
	
	def dxdx(self, i, end=0):
		"""Return derivative of bead i with respect to either end-point 0 or -1"""
		assert end in [0,-1]
		N = self.N
		f = self.f
		a = np.zeros((f,f))
		for j in xrange(f): a[j] = solve(self.JF, self.delta(end,j)/self.eps[end]).reshape(N-1,f)[i-1] / N
		return a.T

	def d2Sdx2(self):
		m = self.m
		N = self.N
		f = self.f
		d2Sdxif2 = np.zeros((2*f,2*f))
		a = np.zeros((f,f))
		# d2S/dx0^2
#		for j in xrange(f): a[j] = linalg.solve(N*self.J, self.delta(0,j)/self.eps[0]).reshape(N-1,f)[0]
		for j in xrange(f): a[j] = solve(self.JF, self.delta(0,j)/self.eps[0]).reshape(N-1,f)[0] / N
		d2Sdxif2[:f,:f] = self.dt[0]/2*self.d2Vdx2[0] + m/self.dt[0] * (np.identity(f) - a)
		# d2S/dx0dxN
#		for j in xrange(f): a[j] = linalg.solve(N*self.J, self.delta(-1,j)/self.eps[-1]).reshape(N-1,f)[0]
		for j in xrange(f): a[j] = solve(self.JF, self.delta(-1,j)/self.eps[-1]).reshape(N-1,f)[0] / N
#		print a
#		p = np.zeros((f,f))
#		for j in xrange(f):
#			for k in xrange(f):
#				minor = np.delete(np.delete(self.J, 2*(N-2)+j, 0), k, 1)
#				p[j,k] = linalg.det(self.eps[-1]*minor)
#				a[j,k] = linalg.det(self.eps[-1]*minor) / linalg.det(self.eps[-1]*self.J) # Cramer's rule
#		print a
#		print linalg.det(a), linalg.det(p) / linalg.det(self.eps[-1]*self.J)**f
#		T = self.eps[-1]*self.J[:-f,f:]
#		print linalg.det(self.eps[-1]*self.J)**(f-1), linalg.det(p)/linalg.det(T)
#		print linalg.det(T), np.prod(-self.eps[-1]/self.eps[:])**f
#		print linalg.det(a), np.prod(-self.eps[-1]/self.eps[:])**f / linalg.det(self.eps[-1]*self.J)
#		print linalg.det(a), np.prod(-1/self.eps[:])**f * self.eps[-1]**(N*f) / linalg.det(self.J) / self.eps[-1]**((N-1)*f)
#		print linalg.det(a), np.prod(-1/self.eps[:])**f * self.eps[-1]**f / linalg.det(self.J)
#		exit(0)
		d2Sdxif2[:f,f:] = - m/self.dt[0]*a.T
		d2Sdxif2[f:,:f] = - m/self.dt[0]*a
		# d2S/dxN^2
#		for j in xrange(f): a[j] = linalg.solve(N*self.J, self.delta(-1,j)/self.eps[-1]).reshape(N-1,f)[-1]
		for j in xrange(f): a[j] = solve(self.JF, self.delta(-1,j)/self.eps[-1]).reshape(N-1,f)[-1] / N
		d2Sdxif2[f:,f:] = self.dt[-1]/2*self.d2Vdx2[-1] + m/self.dt[-1] * (np.identity(f) - a)
#		if self.linear:
#			d2Sdxif2[:f,:f] -= self.dt[0]**3/(12.*self.m) * (-(self.V[1]-self.V[0])*self.dVdx[0]/self.adx[0]**2 - (self.V[1]-self.V[0])**2*(self.x[0]-self.x[1])/self.adx[0]**4)
#			d2Sdxif2[f:,f:] -= self.dt[-1]**3/(12.*self.m) * ((self.V[-1]-self.V[-2])*self.dVdx[-1]/self.adx[-1]**2 - (self.V[-1]-self.V[-2])**2*(self.x[-1]-self.x[-2])/self.adx[-1]**4)
		return d2Sdxif2
	
	def dSdt(self):
		dx = self.x[1:] - self.x[:-1]
		dSdt = - self.m*np.sum(dx**2/self.dt.reshape(-1,1))/(2*self.t) + np.sum(self.eps*(self.V[1:] + self.V[:-1]))/2
		if self.linear:
			dSdt -= np.sum(self.kappa**2 * self.eps**3) / (8. * self.m) * self.t**2
		return dSdt

	def E(self): return self.dSdt()
	
	def dxdt(self, i):
		m = self.m
		c = - (self.dt[:-1]+self.dt[1:]).reshape(-1,1)/m * self.dVdx[1:-1]
		b = solve(self.JF, c.ravel()).reshape(self.N-1,self.f) / self.N
		return b[i-1]

	def d2Sdxdt(self):
		m = self.m
		dx = self.x[1:] - self.x[:-1]
#		c = m/self.t * (dx[:-1]/self.dt[:-1].reshape(-1,1) - dx[1:]/self.dt[1:].reshape(-1,1)) - (self.dt[1:] + self.dt[:-1]).reshape(-1,1)/(2*self.t)*self.dVdx[1:-1]
		c = - (self.dt[:-1]+self.dt[1:]).reshape(-1,1)/m * self.dVdx[1:-1]
#		b = linalg.solve(self.N*self.J, c.ravel()).reshape(self.N-1,self.f)
		b = solve(self.JF, c.ravel()).reshape(self.N-1,self.f) / self.N
		d2Sdxidt = self.dt[0]/(2*self.t)*self.dVdx[0] - self.m*(self.x[0] - self.x[1]) / (self.dt[0]*self.t) - self.m/self.dt[0]*b[0]
		d2Sdxfdt = self.dt[-1]/(2*self.t)*self.dVdx[-1] - self.m*(self.x[-1] - self.x[-2]) / (self.dt[-1]*self.t) - self.m/self.dt[-1]*b[-1]
		d2Sdxdt = np.concatenate([d2Sdxidt, d2Sdxfdt])
		if self.linear:
			d2Sdxdt[:self.f] -= self.kappa[0]*self.eps[0]**3*self.t**2/(4.*self.m) * (self.adx[0]*-self.dVdx[0] - (self.V[1]-self.V[0])*(self.x[0]-self.x[1])/self.adx[0])/self.adx[0]**2
			d2Sdxdt[self.f:] -= self.kappa[-1]*self.eps[-1]**3*self.t**2/(4.*self.m) * (self.adx[-1]*self.dVdx[-1] - (self.V[-1]-self.V[-2])*(self.x[-1]-self.x[-2])/self.adx[-1])/self.adx[-1]**2
		return d2Sdxdt
	
	def d2Sdt2(self):
		m = self.m
		dx = self.x[1:] - self.x[:-1]
		eps = (self.dt / self.t).reshape(-1,1)
#		print 'dS/dt =', - m*np.sum(dx**2/eps)/(2*self.t**2) + np.sum(eps.ravel()*(self.V[1:] + self.V[:-1]))/2
#		print 'dS/dt =', - m*np.sum(dx**2/eps)/(2*self.t**2) + eps[-1,0]*self.V[-1]/2 + np.sum(eps[:-1].ravel()*self.V[1:-1])/2 + eps[0,0]*self.V[0]/2 + np.sum(eps[1:].ravel()*self.V[1:-1])/2
#		c = m/self.t * (dx[:-1]/self.dt[:-1].reshape(-1,1) - dx[1:]/self.dt[1:].reshape(-1,1)) - (self.dt[1:] + self.dt[:-1]).reshape(-1,1)/(2*self.t)*self.dVdx[1:-1]
		c = - (self.dt[:-1]+self.dt[1:]).reshape(-1,1)/m * self.dVdx[1:-1]
#		b = linalg.solve(self.J, c.ravel()).reshape(self.N-1,self.f) / self.N
		b = solve(self.JF, c.ravel()).reshape(self.N-1,self.f) / self.N
		ddx = np.zeros_like(dx)
		ddx[:-1] += b
		ddx[1:] -= b
		d2Sdt2 = self.m*np.sum(dx**2/eps)/self.t**3 - self.m*np.sum(dx*ddx/eps)/self.t**2 + np.sum(eps[:-1]*self.dVdx[1:-1]*b + eps[1:]*self.dVdx[1:-1]*b)/2
		if self.linear:
			d2Sdt2 -= np.sum(self.kappa**2 * self.eps**3) / (4. * self.m) * self.t
		return d2Sdt2
	
	def Gelfand_Yaglom(self, BC=0, retall=False):
		"""Return fxf matrix computed using Gelfand-Yaglom method.
		The determinant of this is equal to |J| and hence related to C = |-d2Sdx'dx''|
		BC is 0 or 1, for initial value of D.
		"""
		dt = np.mean(self.dt)
		assert np.all(abs(self.dt-dt)/dt<1e-8) # ensure all dts are the same
		D = np.zeros((self.N-1,self.f,self.f))
		I = np.identity(self.f)
		if BC == 0:
			D[0] = 2*I + dt**2*self.d2Vdx2[1]
			if self.N==2: return D[0]
			D[1] = np.dot(2*I + dt**2*self.d2Vdx2[2], D[0]) - I # N.B. this order is correct, not as in Kleinert
		elif BC == 1:
			# I derived these so that they return cosh(N*dt*wN) and dt*w/2=sinh(dt*wN/2) for the harmonic system d2Vdx2=w^2
			# the same pattern is applied to each D, even to D[0] and D[1]
			D0 = I + dt**2*self.d2Vdx2[0]/2
			D[0] = np.dot(2*I + dt**2*self.d2Vdx2[1], D0) - I
			if self.N==2: return D[0]
			D[1] = np.dot(2*I + dt**2*self.d2Vdx2[2], D[0]) - D0
		else:
			print 'ERROR: unknown boundary condition: %s' % BC
			exit(1)
		for i in xrange(2,self.N-1):
			D[i] = np.dot(2*I + dt**2*self.d2Vdx2[i+1], D[i-1]) - D[i-2]
		if retall:
			return D
		else:
			return D[-1]

	def monodromy(self, method='RK4'):
		"""Return monodromy matrix
		assuming all dts are equal, and that N is even"""
		dt = np.mean(self.dt)
		# np.set_printoptions(precision=6)
		assert np.all(abs(self.dt-dt)/dt<1e-8) # ensure all dts are the same
		assert self.N % 2 == 0
		R = np.identity(2*self.f, np.complex) # at time 0
		F = np.zeros((2*self.f,2*self.f))
		F[:self.f,self.f:] = - np.identity(self.f) / self.m
		if method == 'Euler': # Euler method
			for i in xrange(self.N):
				F[self.f:,:self.f] = self.d2Vdx2[i]
				R -= 1j*self.dt[i] * np.dot(F, R)
			# print "Euler: \n", R, "\n"
		elif method == 'Midpoint': # midpoint method
			for i in xrange(0,self.N,2):
				F[self.f:,:self.f] = self.d2Vdx2[i]
				k1 = -2j*dt * np.dot(F, R)
				F[self.f:,:self.f] = self.d2Vdx2[i+1]
				k2 = -2j*dt * np.dot(F, R+0.5*k1)
				R = R + k2
			# print "Midpoint: \n", R, "\n"
		elif method == 'RK4': # Runge-Kutta
			for i in xrange(0,self.N,2):
				F[self.f:,:self.f] = self.d2Vdx2[i]
				k1 = -2j*dt * np.dot(F, R)
				F[self.f:,:self.f] = self.d2Vdx2[i+1]
				k2 = -2j*dt * np.dot(F, R+0.5*k1)
				k3 = -2j*dt * np.dot(F, R+0.5*k2)
				F[self.f:,:self.f] = self.d2Vdx2[i+2]
				k4 = -2j*dt * np.dot(F, R+k3)
				R = R + k1/6. + k2/3. + k3/3. + k4/6.
			# print "\nRungeKutta: \n", R, "\n"
		elif method == 'Ordexp':  # Analytical exponential
			R = np.zeros((self.N+1, 2*self.f, 2*self.f), np.complex)
			F[self.f:, :self.f] = self.d2Vdx2[0]
			# R[0] = linalg.expm(-F * 1j*self.dt[0])
			R[0] = np.identity(4)
			for i in xrange(1, self.N+1):
				F[self.f:, :self.f] = self.d2Vdx2[i]
				R[i] = np.dot(linalg.expm(-F * 1j*self.dt[0]), R[i-1])
			R = R[-1]
			# print "Exponential: \n", R, "\n"
		elif method == 'Action':  # Action derivative method
			d2Sdx2 = self.d2Sdx2()
			a = d2Sdx2[:self.f,:self.f]
			b = d2Sdx2[:self.f,self.f:]
			c = d2Sdx2[self.f:,self.f:]
			invb = np.linalg.inv(b)
			R = np.zeros((2 * self.f, 2 * self.f), np.complex)
			R[:self.f,:self.f] = -invmul(b, a) # equals -np.dot(invb, a)
			R[:self.f,self.f:] = -1j*invb
			R[self.f:,:self.f] = -1j*b.T + 1j*np.dot(c,invmul(b,a))
			R[self.f:,self.f:] = -invmul(b.T,c.T).T # is same as -np.dot(c, invb)
			# print "\nAction: \n", R, "\n"
		else:
			print 'WARNING: unknown method %s' % method
			exit(1)

#		print 'F', F
#		print 'R', R
		return R
	
	def stability(self, method='RK4', monodromy=False):
		"""Return stability parameters
		assuming all dts are equal, and that N is even.

		method for calculating monodromy matrix R:
		Euler = 'Euler'
		Midpoint = 'Midpoint'
		Runge-Kutta 4th Order = 'RK4'
		Ordered Exponential Wavefunction = 'Ordexp'
		d2Sdx2 Approach = 'Action'
		"""
		R = self.monodromy(method)
		eigvals = linalg.eigvals(R)
#		print 'eigvals', eigvals
#		print 'log eigvals', np.log(eigvals)
		params = np.log(eigvals).real
		params.sort()
		if monodromy:
			return params[self.f:], R, eigvals
		else:
			return params[self.f:]
	
	def W(self):
		"""Return W which should equal S-Et in infinite N limit"""
		r = pathlength.coord(self.x)
		E = self.E()
		tmp = self.V - E
		tmp = np.where(tmp > 0, tmp, 0)
		integrand = np.sqrt(2*self.m*tmp)
		return integrate.simps(integrand, r)

########################################

if __name__ == "__main__":
	import math
	N = 24
	def V(x): return 0.5*x**2
	def dVdx(x): return x
	def d2Vdx2(x): return np.ones_like(x)
	if 1:
		t = np.linspace(0,3.5,N+1)
	else:
		dt = 0.001 + abs(np.tan(np.linspace(0,2*math.pi,16)))
		dt *= 3.5 / sum(dt)
		print dt
		t = np.cumsum(dt)
	dt = np.diff(t)
	x0, p0 = 2., -1.9
	x = x0*np.cosh(t) + p0*np.sinh(t)
#	print min(x), V(min(x)), min(V(x))
	print 'end points', x[0], x[-1]
	if 0:
		turn = min(x)
		ends = x[0], x[-1]
		x = np.concatenate((np.linspace(ends[0],turn,N/2),np.linspace(turn,ends[1],N/2)[1:])) 
		E = V(turn)
		Vx = V(x)
		p = np.sqrt(2*abs(Vx - E))
#		dpdE = - self.mass / p * np.sign(V-self.E)
		dp3dE = - 3 * p * np.sign(Vx - E) # d(p^3)/dE = 3 p^2 dp/dE
		dx = np.diff(x, axis=0)
		adx = abs(dx)
		slope = (Vx[1:] - Vx[:-1]) / adx
		dt = - (dp3dE[1:] - dp3dE[:-1]) / (3*slope)
#		print dt
		t = np.concatenate([[0],np.cumsum(dt)])
		print 'time = %g' % np.sum(dt)
	############## exact #########
	E = - p0**2/2 + V(x0)
	print 
	sinh, cosh = math.sinh(t[-1]), math.cosh(t[-1])
	S = ((x[0]**2+x[-1]**2)*cosh - 2*x[0]*x[-1]) / (2*sinh)
	dSdt = (-(x[0]**2+x[-1]**2)/2 + x[0]*x[-1]*cosh) / sinh**2
	dEdx = np.array([-x[0] + x[-1]*cosh, -x[-1] + x[0]*cosh]) / sinh**2
	d2Sdt2 = (sinh**2*x[0]*x[-1]*sinh - (-(x[0]**2+x[-1]**2)/2 + x[0]*x[-1]*cosh)*2*sinh*cosh) / sinh**4
	dSdx = np.array([x[0]*cosh - x[-1], x[-1]*cosh - x[0]]) / sinh
	d2Sdx2 = np.array([[cosh,-1],[-1,cosh]]) / sinh
	############## numeric #########
	print 'trajectory without linear'
	Traj = Trajectory(x.reshape(-1,1), V(x), dVdx(x).reshape(-1,1), d2Vdx2(x).reshape(-1,1,1), dt, linear=False)
	print 'S = %g ~ %g' % (Traj.S(), S)
	print 'E = dS/dt = %g ~ %g = %g' % (Traj.dSdt(), dSdt, E)
	print 'dE/dt = d2S/dt2 = %g ~ %g' % (Traj.d2Sdt2(), d2Sdt2)
	print 'dS/dx = %s ~ %s' % (Traj.dSdx(), dSdx)
	print 'dE/dx = %s ~ %s' % (Traj.d2Sdxdt(), dEdx)
	print 'd2S/dx2 = %s ~ %s' % (Traj.d2Sdx2().ravel(), d2Sdx2.ravel())
	F = np.zeros((2,2))
	F[:1,1:] = -np.identity(1)
	F[1:,:1] = d2Vdx2(0)
	print 'F', F
	R = linalg.expm(-1j*F*Traj.t)
	print 'R', R
	eigvals = linalg.eigvals(R)
	print 'eigvals', eigvals
	print 'log eigvals', np.log(eigvals)
	print 'stability parameters', Traj.stability(), d2Vdx2(0)*t[-1]

	print '\ntrajectory with linear'
	Traj = Trajectory(x.reshape(-1,1), V(x), dVdx(x).reshape(-1,1), d2Vdx2(x).reshape(-1,1,1), dt, linear=True)
	print 'S = %g ~ %g' % (Traj.S(), S)
	print 'E = dS/dt = %g ~ %g = %g' % (Traj.dSdt(), dSdt, E)
	print 'dE/dt = d2S/dt2 = %g ~ %g' % (Traj.d2Sdt2(), d2Sdt2)
	print 'dS/dx = %s ~ %s' % (Traj.dSdx(), dSdx)
	print 'dE/dx = %s ~ %s' % (Traj.d2Sdxdt(), dEdx)
	print 'd2S/dx2 = %s ~ %s' % (Traj.d2Sdx2().ravel(), d2Sdx2.ravel())

	import pylab
	fig, (axx, axt) = pylab.subplots(2)
	grid = np.linspace(-1,2.5,100)
	axx.plot(grid,V(grid))
	axx.plot(x, E*np.ones_like(x), 'x-')
	axx.set_ylim(0,2*E)
	axt.plot(t,x,'x-')
	pylab.show()
