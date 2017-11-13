# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all clean harmonic qtip4pf
all: harmonic qtip4pf

driver:=i-pi-driver
IPI:=i-pi

define run_driver
  for i in `seq 1 $1`; do \
    $(driver) -u -h $2 -m harm -o 1 & \
  done;
endef
-include make.in

harmonic:
	cd harmonic; $(IPI) input.xml & sleep 5; \
	$(call run_driver,4, harmonic) \
	wait
	
qtip4pf:
	cd qtip4pf; $(IPI) input.xml & sleep 5; \
	$(call run_driver,4, qtip4pf) \
	wait

clean:
	rm -f */*simulation* */*RESTART* */EXIT */*.log*
 
