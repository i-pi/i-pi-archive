# Makefile for the RPMD post-processing files
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

.PHONY: all scripts clean run_1 run_2 run_3 run_4 run_5

all: run_1 run_2 run_3 run_4 run_5
scripts: consolidate vel_est

driver_dir:=../../../../drivers

consolidate: consolidate.f
	gfortran -O3 -o consolidate consolidate.f 
	
vel_est: vel_est.f
	gfortran -O3 -o vel_est vel_est.f
	
clean:
	rm consolidate vel_est; bash mkclean.sh

run_1:
	cd run_1; bash ../launch.sh 31415 & sleep 5; \
	cd $(driver_dir); bash launch_multi.sh 31415; \
	wait

run_2:
	cd run_2; bash ../launch.sh 31414 & sleep 5; \
	cd $(driver_dir); bash launch_multi.sh 31414; \
	wait

run_3:
	cd run_3; bash ../launch.sh 31413 & sleep 5; \
	cd $(driver_dir); bash launch_multi.sh 31413; \
	wait

run_4:
	cd run_4; bash ../launch.sh 31412 & sleep 5; \
	cd $(driver_dir); bash launch_multi.sh 31412; \
	wait

run_5:
	cd run_5; bash ../launch.sh 31411 & sleep 5; \
	cd $(driver_dir); bash launch_multi.sh 31411; \
	wait
