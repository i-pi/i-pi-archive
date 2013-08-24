.PHONY: all clean rpc
all: rpc

LAMMPS:=lmp_serial
IPI:=../../../i-pi


rpc:
	$(IPI) input.xml & sleep 5; \
        $(LAMMPS) < in.water_longrange & \
        $(LAMMPS) < in.water_shortrange; \
	wait
	
clean:
	rm -f RESTART EXIT log.lammps rpc.*
