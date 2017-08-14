# Makefile for the CP2K example
#
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all clean dia-npt-classical dia_big-npt-classical dia-npt-pimd.4 dia_big-npt-pimd.4 water lithium
all: dia-npt-classical dia_big-npt-classical dia-npt-pimd.4 dia_big-npt-pimd.4

DRIVER:=pw.x
IPI:=i-pi

-include make.in

define run_driver
  for i in `seq 1 $1`; do \
    mkdir -p run_$$i; cd run_$$i; \
    $(DRIVER) < ../$2 > pw.out --ipi $3 & cd ..; \
	sleep 3; \
  done;
endef

water:
	cd h2o; $(IPI) ipi.xml & sleep 3; \
	$(call run_driver,4,pw.in,localhost:3141) \
	wait

lithium:
	cd li4; $(IPI) ipi.xml & sleep 3; \
	$(call run_driver,4,pw.in,localhost:31415) \
	wait

dia-npt-classical:
	cd dia-npt-classical; $(IPI) ipi.xml & sleep 5; \
	$(call run_driver,1,pw.in, localhost:3142) \
	wait

dia-npt-pimd.4:
	cd dia-npt-pimd.4; $(IPI) ipi.xml & sleep 5; \
	$(call run_driver,4,pw.in, localhost:3143) \
	wait

dia_big-npt-classical:
	cd dia_big-npt-classical; $(IPI) ipi.xml & sleep 5; \
	$(call run_driver,1,pw.in, localhost:3144) \
	wait

dia_big-npt-pimd.4: 
	cd dia_big-npt-pimd.4; $(IPI) ipi.xml & sleep 5; \
	$(call run_driver,4,pw.in, localhost:3145) \
	wait

npt-piglet:
	cd npt-piglet; $(IPI) nptpgl_ipi.xml & sleep 10; \
	$(call run_driver,4,pw.in) \
	wait

nvt-piglet:
	cd nvt-piglet; $(IPI) nvtpgl_ipi.xml & sleep 10; \
	$(call run_driver,4,pw.in) \
	wait

clean:
	rm -rf */\#*\#* */example.* */*.out */*.log */run_* */RESTART */EXIT
