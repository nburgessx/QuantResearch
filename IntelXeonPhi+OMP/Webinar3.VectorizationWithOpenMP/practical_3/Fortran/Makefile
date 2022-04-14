# The file containing the main program
cheby_src:=cheby.f90
# Compiler optimization flags to be used
OPT_FLAGS:=-O3 -qopenmp

# No need to edit below here

FC:=ifort
CC:=icc
LD:=$(FC)

FFLAGS:=-stand=f03 $(OPT_FLAGS)
CFLAGS:=-I . -std=c99 $(OPT_FLAGS)
LDFLAGS:=$(OPT_FLAGS) 

source:=util.f90 $(cheby_src) main.f90
main:=$(patsubst %.c,%.exe,$(cheby_src)) 
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
