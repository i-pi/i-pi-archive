"""
Calculate the accumulated average of the exponential estimator <exp(-beta*spr)>
for the isotope fractionation. 

accumulated.out:
first column : simulation time
second column : the instantanous exponential estimator
third column : the accumulated average of the exponential estimator

"""

import numpy as np

# the number of steps that allows the system to equilibrate
equilibriumstep = 25000  
 
table = np.loadtxt('no_rpc.out')
 
(asize, bsize) = np.shape(table)
print asize, bsize

out = np.zeros((asize-equilibriumstep, 3))

count = 0

for a in range(asize):
	
	if a >= equilibriumstep:
		out[count][0] = table[a][1]
		out[count][1] = table[a][9]
		if count == 0:
			out[count][2] = table[a][9]
		else:
			out[count][2] = out[count-1][2]*count/(count+1) +table[a][9]/(count+1)					
		count += 1

np.savetxt('accumulated.out', out)

