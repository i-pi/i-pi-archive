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

graphene:
	nohup $(IPI) input.xml & \
	sleep 25; \
        $(call run_lammps,4,in.kappa.Graphene) \
	wait

clean:
	rm -f graph.* RESTART EXIT log.lammps graphene.lammpstraj nohup.out
