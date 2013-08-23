.PHONY: all clean graphene
all: graphene

LAMMPS:=~/local/dpi_lammps/lammps-gle/src/lmp_serial
IPI:=../../../i-pi


graphene:
	nohup $(IPI) input.xml & \
	sleep 25; \
        $(LAMMPS) < in.kappa.Graphene; \
	wait

clean:
	rm -f graph.* RESTART EXIT log.lammps graphene.lammpstraj nohup.out
