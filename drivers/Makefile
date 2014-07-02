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
 
.PHONY: all modules clean

FLAGS=-g -O3 -Wall -Wno-maybe-uninitialized
MODULES=distance.f90 LJ.f90 SG.f90
all: modules driver.x driver_pure.x

modules: $(MODULES)
	gfortran $(FLAGS) -c $(MODULES)

sockets.o: sockets.c
	gcc $(FLAGS) -c -o sockets.o sockets.c

driver.x: $(MODULES) sockets.o fsockets.o driver.o
	gfortran $(FLAGS) -o driver.x $^

driver_pure.x: $(MODULES) fsockets_pure.o driver.o 
	gfortran $(FLAGS) -o driver_pure.x $^
	
%.o: %.f90
	gfortran $(FLAGS) -c $<

clean:
	rm *.o *.mod *.x
