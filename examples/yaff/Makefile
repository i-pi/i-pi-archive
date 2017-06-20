# Makefile for the yaff example
#
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all clean mil53
all: mil53

IPI:=i-pi

define run_yaff
  for i in `seq 1 $1`; do \
    python run.py & \
  done;
endef
-include make.in

mil53:
	cd MIL53; $(IPI) nst_ipi.xml & sleep 10; \
	$(call run_yaff,4) \
	wait

clean:
	rm */*.out* */*.log */*.restart* */*RESTART* */EXIT */*pos* */yaffdriver.pyc
