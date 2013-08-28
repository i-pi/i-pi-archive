# Makefile for the para-hydrogen examples
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

.PHONY: all clean mht-25K mht-20K mht-17K mht-14K nvt rpmd-run_1 rpmd-run_2 rpmd-run_3 rpmd-run_4 rpmd-run_5
all: mht-25K mht-20K mht-17K mht-14K nvt rpmd-run_1 rpmd-run_2 rpmd-run_3 rpmd-run_4 rpmd-run_5

driver_dir:=../../../../drivers/
driver:=../../../drivers/driver.x
IPI:=../../../i-pi

define run_driver
  for i in `seq 1 $1`; do \
    $(driver) -u -h $2 -m sg -o 15.0 & \
  done
endef
-include make.in

nvt:
	cd nvt; $(IPI) input.xml & sleep 5; \
	$(call run_driver,4,ph2_nvt); \
	wait

mht-25K:
	cd MHT-25K; $(IPI) input_equilibrated.xml & sleep 5; \
	$(call run_driver,4,T_25K); \
	wait
	
mht-20K:
	cd MHT-20K; $(IPI) input_equilibrated.xml & sleep 5; \
	$(call run_driver,4,T_20K); \
	wait

mht-17K:
	cd MHT-17K; $(IPI) input_equilibrated.xml & sleep 5; \
	$(call run_driver,4,T_17K); \
	wait

mht-14K:
	cd MHT-14K; $(IPI) input_equilibrated.xml & sleep 5; \
	$(call run_driver,4,T_14K); \
	wait

# this rule creates a run_% subdirectory of RPMD, and executes 16 RPMD runs in it,
# always restarting from the previous run with randomized velocities.
RPMD/run_%/test16.out: 
	cd RPMD; mkdir -p run_$*; cp start-$*.chk run_$*/RESTART; \
	for x in `seq 1 16`; do \
	  sed "s/YYY/"$*"/; s/XXX/$$x/;" input.xml > run_$*/input.xml; \
	  $(IPI) run_$*/input.xml & sleep 25; \
	  $(call run_driver,8,run-$*-$$x); wait; \
	  mv RESTART run_$*/; \
	done  
	

rpmd: RPMD/run_1/test16.out RPMD/run_2/test16.out RPMD/run_3/test16.out RPMD/run_4/test16.out
	make -C RPMD scripts; make -C RPMD postproc
	
	
rpmd-run_1:
	cd RPMD/run_1; bash ../launch.sh 31415 & sleep 5; \
	cd $(driver_dir); bash launch_multi.sh 31415; \
	wait

rpmd-run_2:
	cd RPMD/run_2; bash ../launch.sh 31414 & sleep 5; \
	cd $(driver_dir); bash launch_multi.sh 31414; \
	wait

rpmd-run_3:
	cd RPMD/run_3; bash ../launch.sh 31413 & sleep 5; \
	cd $(driver_dir); bash launch_multi.sh 31413; \
	wait

rpmd-run_4:
	cd RPMD/run_4; bash ../launch.sh 31412 & sleep 5; \
	cd $(driver_dir); bash launch_multi.sh 31412; \
	wait

rpmd-run_5:
	cd RPMD/run_5; bash ../launch.sh 31411 & sleep 5; \
	cd $(driver_dir); bash launch_multi.sh 31411; \
	wait

clean:
	rm -f */ph2-* */RESTART */EXIT; \
	make -C RPMD clean; rm -rf RPMD/run_*
