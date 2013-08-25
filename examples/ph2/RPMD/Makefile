scripts: consolidate vel_est

consolidate: consolidate.f
	gfortran -O3 -o consolidate consolidate.f 
	
vel_est: vel_est.f
	gfortran -O3 -o vel_est vel_est.f
	
clean:
	rm -f consolidate vel_est
