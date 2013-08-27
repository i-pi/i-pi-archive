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

.PHONY: all clean graphene h2o-piglet_2 h2o-piglet_4 h2o-piglet_8 rpc h2o-rpc h2o-pimd
all: graphene h2o-piglet_2 h2o-piglet_4 h2o-piglet_8 h2o-rpc h2o-pimd

LAMMPS:=lmp_serial
IPI:=../../../i-pi

define run_lammps
  for i in `seq 1 $1`; do \
    $(LAMMPS) < $2 & \
  done;
endef
-include make.in

graphene:
	cd graphene; $(IPI) input.xml & \
	sleep 40; \
        $(call run_lammps,4,in.kappa.Graphene) \
	wait

h2o-piglet_2:
	cd h2o-piglet.2; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,2,in.water) \
	wait
	
h2o-piglet_4:
	cd h2o-piglet.4; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.water) \
	wait

h2o-piglet_8:
	cd h2o-piglet.8; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.water) \
	wait

h2o-pimd:
	cd h2o-pimd; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.water) \
	wait

h2o-rpc:
	cd h2o-pimd-rpc; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.water_shortrange) \
        $(call run_lammps,1,in.water_longrange) \
	wait

clean:
	rm -f */gle_lammps.* */RESTART */EXIT */log.lammps; \
	cd h2o-pimd; \
	rm -f RESTART EXIT log.lammps no_rpc.*; cd ..; \
	cd h2o-pimd-rpc; \
	rm -f RESTART EXIT log.lammps rpc.*; cd ..; \
	cd graphene; \
	rm -f graph.* RESTART EXIT log.lammps graphene.lammpstraj nohup.out; cd ..
