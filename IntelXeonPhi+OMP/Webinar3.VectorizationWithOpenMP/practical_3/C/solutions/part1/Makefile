# The file containing the main program
main_src:=cheby.c
# Compiler optimization flags to be used
OPT_FLAGS:=-O3 -qopenmp -qopt-report-phase=vec -qopt-report=5

# No need to edit below here

FC:=ifort
CC:=icc
LD:=$(FC)

FFLAGS:=$(OPT_FLAGS)
CFLAGS:=-I . -std=c99 $(OPT_FLAGS)
LDFLAGS:=$(OPT_FLAGS) -nofor-main 

source:=$(main_src) eps.c epsilon.f90
main:=$(patsubst %.c,%.exe,$(main_src)) 
main:=$(patsubst %.f90,%.exe,$(main)) 
objects:=$(patsubst %.f90,%.o,$(source)) 
objects:=$(patsubst %.f,%.o,$(objects)) 
objects:=$(patsubst %.c,%.o,$(objects)) 
#objects:=$(sort $(objects))

$(main): $(objects)
	$(LD) -o $@ $^ $(LDFLAGS)

%.o: %.f90 Makefile
	$(FC) -c -o $@ $(FFLAGS) $<

%.o: %.f Makefile
	$(FC) -c -o $@ $(FFLAGS) $<

%.o: %.c Makefile
	$(CC) -c -o $@ $(CFLAGS) $<

clean:
	rm -f *.exe *.o *.mod *.modmic
