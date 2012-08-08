.PHONY: all

all: consolidate mom_est

consolidate: consolidate.f
	gfortran -O3 -o consolidate consolidate.f 
	
mom_est: mom_est.f
	gfortran -O3 -o mom_est mom_est.f
