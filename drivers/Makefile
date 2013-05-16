.PHONY: all

SOURCES=$(shell ls *.f90)
OBJECTS=$(SOURCES:.f90=.o) sockets.o
all: driver.x

sockets.o: sockets.c
	gcc -c -o sockets.o sockets.c

driver.x: $(OBJECTS)
	gfortran -O3 -o driver.x $(OBJECTS)

%.o: %.f90
	gfortran -O3 -c $<

clean:
	rm *.o *.mod *.x
