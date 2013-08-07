# Makefile for the driver tests
# 
# Copyright (C) 2013, Joshua More and Michele Ceriotti
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http.//www.gnu.org/licenses/>.

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
