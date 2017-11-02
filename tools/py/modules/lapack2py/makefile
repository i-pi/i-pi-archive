MODULE=lapack2py
FC=gfortran
LIBS=/usr/lib/liblapack.so
PROJECTS=$(HOME)/projects

all: $(MODULE).so

$(MODULE).so: lapack.f90
	f2py -c -m $(MODULE) $^ $(LIBS)

%.o: %.f90
	$(FC) -c $^ -O2 -fPIC #$(INCLUDE)
%.o: %.F90
	$(FC) -c $^ -O2 -fPIC #$(INCLUDE)

install: $(MODULE).so
	install $^ $(PROJECTS)/lib/f2py

clean:
	$(RM) -f *.o *.so *.mod *.c

again: clean all
