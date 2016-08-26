# Makefile for the CP2K example
#
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all clean npt-classical npt-piglet nvt-classical nvt-piglet
all: npt-classical npt-piglet nvt-classical nvt-piglet

CP2K:=cp2k.x
IPI:=i-pi

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
