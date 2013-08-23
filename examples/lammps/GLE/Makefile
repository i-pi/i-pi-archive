.PHONY: all clean piglet_2 piglet_4 piglet_8
all: piglet_2 piglet_4 piglet_8

LAMMPS:=~/local/dpi_lammps/lammps-gle/src/lmp_serial
IPI:=../../../../i-pi


piglet_2:
	cd piglet_2; $(IPI) input.xml & sleep 5; \
        $(LAMMPS) < in.water; \
	wait
	
piglet_4:
	cd piglet_4; $(IPI) input.xml & sleep 5; \
        $(LAMMPS) < in.water; \
	wait

piglet_8:
	cd piglet_8; $(IPI) input.xml & sleep 5; \
        $(LAMMPS) < in.water; \
	wait

clean:
	rm -f */gle_lammps.* */RESTART */EXIT */log.lammps
