# Makefile for the RPMD post-processing files
# 
# Copyright (C) 2013, Joshua More and Michele Ceriotti
# 
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
# 
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

.PHONY: clean scripts postproc

scripts: consolidate vel_est

consolidate: consolidate.f
	gfortran -O3 -o consolidate consolidate.f 
	
vel_est: vel_est.f
	gfortran -O3 -o vel_est vel_est.f

postproc: consolidate vel_est
	for r in `seq 1 2 3 4`; do cd run_$$r; for j in `seq 1 16`; do ../vel_est $$j &> /dev/null; done; cd ..; done; \
	echo "4 16\n" | ./consolidate

clean:
	rm -f consolidate vel_est vel_corr_tot.out */vel_corr*
	
