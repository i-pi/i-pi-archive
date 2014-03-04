# Makefile for the driver tests
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
# 
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 
FC=gfortran
FFLAGS=-O3 -J build
CC=gcc

SRC_F=build/distance.f90 build/LJ.f90 build/SG.f90 build/driver.f90
OBJECTS=$(SRC_F:.f90=.o) build/sockets.o

.PHONY: all

all: build_dir $(OBJECTS) ../bin/i-pi-driver

build_dir:
	@mkdir -p build

build/sockets.o: sockets.c
	$(CC) -c $< -o $@

../bin/i-pi-driver: $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $(OBJECTS)

build/%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -rf build ../bin/i-pi-driver
