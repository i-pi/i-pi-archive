# Makefile for the para-hydrogen examples
#
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all clean mht-25K mht-20K mht-17K mht-14K nvt multi rpmd 
all: mht-25K mht-20K mht-17K mht-14K nvt multi rpmd 

driver:=i-pi-driver
IPI:=i-pi

define run_driver
  for i in `seq 1 $1`; do \
    $(driver) -u -h $2 -m sg -o 30.0 & \
  done
endef

nvt:
	cd nvt; $(IPI) input.xml & sleep 5; \
	$(call run_driver,4,ph2_nvt); \
	wait

multi:
	cd multi; $(IPI) input.xml & sleep 5; \
	$(call run_driver,4,ph2_multi); \
	wait

paratemp:
	cd paratemp; $(IPI) input.xml & sleep 2; \
	$(call run_driver,3,ph2_pt); \
	wait

parawte:
	cd parawte; $(IPI) input.xml & sleep 3; \
	$(call run_driver,3,ph2_wte); \
	wait

paraptwte:
	cd paraptwte; $(IPI) input.xml & sleep 3; \
	$(call run_driver,3,ph2_ptwte); \
	wait

mht-25K:
	cd MHT-25K; $(IPI) input.xml & sleep 5; \
	$(call run_driver,4,T_25K); \
	wait

mht-20K:
	cd MHT-20K; $(IPI) input.xml & sleep 5; \
	$(call run_driver,4,T_20K); \
	wait

mht-17K:
	cd MHT-17K; $(IPI) input.xml & sleep 5; \
	$(call run_driver,4,T_17K); \
	wait

mht-14K:
	cd MHT-14K; $(IPI) input.xml & sleep 5; \
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

clean:
	rm -f */ph2* */RESTART */EXIT; \
	rm -f para*/*_*ph2* para*/PARA*; \
	make -C RPMD clean; rm -rf RPMD/run_*
