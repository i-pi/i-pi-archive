.PHONY: all clean graphene
all: graphene

LAMMPS:=lmp_serial
IPI:=../../../i-pi

define run_lammps
  for i in `seq 1 $1`; do \
    $(LAMMPS) < $2 & \
  done;
endef
-include make.in

gle-piglet_2:
	cd GLE/piglet_2; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,2,in.water) \
	wait
	
gle-piglet_4:
	cd GLE/piglet_4; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.water) \
	wait

gle-piglet_8:
	cd GLE/piglet_8; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.water) \
	wait

no_rpc:
	cd no_RPC; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.water) \
	wait

rpc:
	cd RPC; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.water_shortrange) \
        $(call run_lammps,1,in.water_longrange) \
	wait

clean:
	cd GLE; \
	rm -f */gle_lammps.* */RESTART */EXIT */log.lammps; cd ..; \
	cd no_RPC; \
	rm -f RESTART EXIT log.lammps no_rpc.*; cd ..; \
	cd RPC; \
	rm -f RESTART EXIT log.lammps rpc.*; cd ..;
