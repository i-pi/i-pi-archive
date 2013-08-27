# Makefile for the para-hydrogen examples
# 
# Copyright (C) 2013, Joshua More and Michele Ceriotti
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http.//www.gnu.org/licenses/>.

.PHONY: all clean tuckerman_25K tuckerman_20K tuckerman_17K tuckerman_14K nvt rpmd-run_1 rpmd-run_2 rpmd-run_3 rpmd-run_4 rpmd-run_5
all: tuckerman_25K tuckerman_20K tuckerman_17K tuckerman_14K nvt rpmd-run_1 rpmd-run_2 rpmd-run_3 rpmd-run_4 rpmd-run_5

driver_dir:=../../../drivers/
driver:=../../../../drivers/driver.x
IPI:=../../../../i-pi

define run_driver
  for i in `seq 1 $1`; do \
    $(driver) -u -h $2 -m sg -o 15.0 & \
  done;
endef
-include make.in

nvt:
	cd nvt; $(IPI) input.xml & sleep 5; \
	$(call run_driver,4,ph2_nvt) \
	wait

tuckerman_25K:
	cd Tuckerman/25K; $(IPI) input_equilibrated.xml & sleep 5; \
	$(call run_driver,4,T_25K) \
	wait
	
tuckerman_20K:
	cd Tuckerman/20K; $(IPI) input_equilibrated.xml & sleep 5; \
	$(call run_driver,4,T_20K) \
	wait

tuckerman_17K:
	cd Tuckerman/17K; $(IPI) input_equilibrated.xml & sleep 5; \
	$(call run_driver,4,T_17K) \
	wait

tuckerman_14K:
	cd Tuckerman/14K; $(IPI) input_equilibrated.xml & sleep 5; \
	$(call run_driver,4,T_14K) \
	wait

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
	cd Tuckerman; rm -f */ph2-* */RESTART */EXIT; cd ..; \
	cd nvt ; rm -f ph2.* RESTART EXIT; cd ..; \
	cd RPMD; bash mkclean.sh
