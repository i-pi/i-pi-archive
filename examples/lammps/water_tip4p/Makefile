# Makefile for the q-TIP4P-F examples
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

.PHONY: all clean gle-piglet_2 gle-piglet_4 gle-piglet_8 rpc no_rpc
all: gle-piglet_2 gle-piglet_4 gle-piglet_8 rpc no_rpc

LAMMPS:=lmp_serial
IPI:=../../../../i-pi

define run_lammps
  for i in `seq 1 $1`; do \
    $(LAMMPS) < $2 & \
  done;
endef
-include make.in

gle-piglet_2:
	cd GLE/piglet_2; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,2,in.water) \
	wait
	
gle-piglet_4:
	cd GLE/piglet_4; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.water) \
	wait

gle-piglet_8:
	cd GLE/piglet_8; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.water) \
	wait

no_rpc:
	cd no_RPC; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.water) \
	wait

rpc:
	cd RPC; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.water_shortrange) \
        $(call run_lammps,1,in.water_longrange) \
	wait

clean:
	cd GLE; \
	rm -f */gle_lammps.* */RESTART */EXIT */log.lammps; cd ..; \
	cd no_RPC; \
	rm -f RESTART EXIT log.lammps no_rpc.*; cd ..; \
	cd RPC; \
	rm -f RESTART EXIT log.lammps rpc.*; cd ..;
