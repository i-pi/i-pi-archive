.PHONY: all

SOURCES=$(shell ls *.f90)
MODULES=distance.f90 LJ.f90 SG.f90
OBJECTS=$(SOURCES:.f90=.o) sockets.o
all: modules sockets.o driver.x

modules: $(MODULES)
	gfortran -O3 -c $(MODULES)

sockets.o: sockets.c
	gcc -c -o sockets.o sockets.c

driver.x: $(OBJECTS)
	gfortran -O3 -o driver.x $(OBJECTS)

%.o: %.f90
	gfortran -O3 -c $<

clean:
	rm *.o *.mod *.x
