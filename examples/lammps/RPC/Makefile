.PHONY: all clean rpc
all: rpc

LAMMPS:=~/local/dpi_lammps/lammps-gle/src/lmp_serial
IPI:=../../../i-pi


rpc:
	$(IPI) input.xml & sleep 5; \
        $(LAMMPS) < in.water_longrange & \
        $(LAMMPS) < in.water_shortrange; \
	wait
	
clean:
	rm -f RESTART EXIT log.lammps rpc.*
