# Makefile for the q-TIP4P-F examples
#
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all clean graphene h2o-piglet_2 h2o-piglet_4 h2o-piglet_8 h2o-pimd-rpc h2o-pimd ice-npt h2o-mts ice-nst isof-vapor isof-water isofsc-vapor isofsc-water 
all: graphene h2o-piglet_2 h2o-piglet_4 h2o-piglet_8 h2o-pimd-rpc h2o-pimd ice-npt h2o-mts ice-nst isof-vapor isof-water isofsc-vapor isofsc-water 

LAMMPS:=lmp_serial
IPI:=i-pi

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

ice-npt:
	cd ice-npt; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.ice) \
	wait

ice-nst:
	cd ice-npt; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.ice) \
	wait

isof-vapor:
	cd isof-vapor; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.vapor) \
	wait

isof-water:
	cd isof-water; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.water) \
	wait

isofsc-vapor:
	cd isofsc-vapor; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.vapor) \
	wait

isofsc-water:
	cd isofsc-water; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.water) \
	wait

h2o-mts:
	cd h2o-mts; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.water_longrange) \
        $(call run_lammps,4,in.water_shortrange) \
	wait

h2o-pimd-rpc:
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
	rm -f graph.* RESTART EXIT log.lammps graphene.lammpstraj nohup.out; cd ..; \
	cd h2o-mts; \
	rm -f simulation.pos_?.xyz simulation.restart; cd ..
