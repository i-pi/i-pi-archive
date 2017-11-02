import numpy as np

def coord(x, mass=1):
	"""Return (mass-weighted) integrated path-length.
	where r[0] = 0
	"""
	xm = x * np.sqrt(mass)
	dxm = np.diff(xm, axis=0)
	dr = np.sqrt(np.sum(dxm.reshape(len(dxm),-1)**2, 1))
	r = [0] + list(np.cumsum(dr))
	r = np.array(r)
	return r

def pathlength(x, mass=1):
	return coord(x, mass)[-1]

######################

if __name__ == "__main__":
	x = np.arange(10)
	print coord(x)
