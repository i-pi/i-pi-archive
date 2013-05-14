.PHONY: all

SOURCES=$(shell ls *.f)
OBJECTS=$(SOURCES:.f=.o)
all: objects ljengine.x oned.x sgengine.x gas_mol.x

objects: $(SOURCES)
	gfortran -O3 -c $(SOURCES)

sgengine.x: $(OBJECTS)
	gfortran -O3 -o sgengine.x ../src/driver/sockets.o $(OBJECTS) interface_SG.F 
	
%.o: %.f
	gfortran -O3 -c $<
	
system_read.o: xml_read.o sys_vars.o
interface.o: system_read.o LJ.o SG.o

oned.x: oned.F
	gfortran -O3 -o oned.x oned.F ../src/driver/sockets.o sys_vars.o

gas_mol.x: gas_mol.F 
	gfortran -O3 -o gas_mol.x gas_mol.F ../src/driver/sockets.o sys_vars.o

ljengine.x: interface_LJ.F
	gfortran -O3 -o ljengine.x ../src/driver/sockets.o $(OBJECTS) interface_LJ.F 

clean:
	rm *.o *.x
