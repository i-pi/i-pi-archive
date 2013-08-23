.PHONY: all clean no_rpc
all: no_rpc

LAMMPS:=lmp_serial
IPI:=../../../i-pi


no_rpc:
	$(IPI) input.xml & sleep 5; \
        $(LAMMPS) < in.water; \
	wait
	
clean:
	rm -f RESTART EXIT log.lammps no_rpc.*
