# Makefile for the Quantum Espresso examples
#
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all water lithium clean
all: water lithium zundel

PW:=pw.x
CP:=cp.x
IPI:=i-pi

define run_pw
  for i in `seq 1 $1`; do \
    mkdir -p run_$$i; cd run_$$i; \
    $(PW) < ../$2 > qe.out & cd ..; \
  done;
endef

define run_cp
  for i in `seq 1 $1`; do \
    mkdir -p run_$$i; cd run_$$i; \
    $(CP) < ../$2 > qe.out & cd ..; \
  done;
endef

water:
	cd h2o; $(IPI) h2o-ipi.xml & sleep 3; \
	$(call run_pw,4,h2o-pw.in) \
	wait

lithium:
	cd li4; $(IPI) li4-ipi.xml & sleep 3; \
	$(call run_pw,4,li4-pw.in) \
	wait

zundel:
	cd zundel-cp; $(IPI) zundel-ipi.xml & sleep 3; \
	$(call run_cp,4,zundel-cp.in) \
	wait

clean:
	rm -rf */*_output.* */*.log */run_* */RESTART */EXIT
