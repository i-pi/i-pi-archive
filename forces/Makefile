.PHONY: all

SOURCES=$(shell ls *.f)
OBJECTS=$(SOURCES:.f=.o)
all: ljengine.x

ljengine.x: $(OBJECTS)
	gfortran -O3 -o ljengine.x ../src/driver/sockets.o $(OBJECTS)
	
%.o: %.f
	gfortran -O3 -c $<
	
system_read.o: xml_read.o sys_vars.o
interface.o: system_read.o LJ.o SG.o
