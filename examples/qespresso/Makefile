# Makefile for the Quantum Espresso examples
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

.PHONY: all water lithium clean
all: water lithium

PW:=pw.x
IPI:=../../../i-pi

define run_qe
  for i in `seq 1 $1`; do \
    mkdir -p run_$$i; cd run_$$i; \
    $(PW) < ../$2 > qe.out & cd ..; \
  done;
endef
-include make.in

water:
	cd h2o; $(IPI) h2o-ipi.xml & sleep 3; \
	$(call run_qe,4,h2o-pw.in) \
	wait

lithium:
	cd li4; $(IPI) li4-ipi.xml & sleep 3; \
	$(call run_qe,4,li4-pw.in) \
	wait
	
clean:
	rm -rf */*_output.* */*.log */run_* */RESTART */EXIT
