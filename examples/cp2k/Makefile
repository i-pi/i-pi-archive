# Makefile for the CP2K example
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

.PHONY: all clean npt-classical npt-piglet nvt-classical nvt-piglet
all: npt-classical npt-piglet nvt-classical nvt-piglet

CP2K:=cp2k.x
IPI:=../../../i-pi

define run_cp2k
  for i in `seq 1 $1`; do \
    mkdir -p run_$$i; cd run_$$i; \
    $(CP2K) -i ../$2 -o cp2k.out & cd ..; \
  done;
endef
-include make.in

npt-classical:
	cd npt-classical; $(IPI) nptcls_ipi.xml & sleep 5; \
	$(call run_cp2k,1,nptcls_cp2k.in) \
	wait
	
nvt-classical:
	cd nvt-classical; $(IPI) nvtcls_ipi.xml & sleep 5; \
	$(call run_cp2k,1,nvtcls_cp2k.in) \
	wait

npt-piglet:
	cd npt-piglet; $(IPI) nptpgl_ipi.xml & sleep 10; \
	$(call run_cp2k,4,nptpgl_cp2k.in) \
	wait

nvt-piglet:
	cd nvt-piglet; $(IPI) nvtpgl_ipi.xml & sleep 10; \
	$(call run_cp2k,4,nvtpgl_cp2k.in) \
	wait

clean:
	rm -rf */*_output.* */*.log */run_* */RESTART */EXIT
