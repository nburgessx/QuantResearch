

CC := icpc
CFLAGS := --std=c++11 -O3 -qopenmp -xMIC-AVX512 -fp-model fast=2 -qopt-report-phase=vec -qopt-report=5

lj : lj.cpp
	$(CC) $(CFLAGS) lj.cpp -o lj



.PHONY : clean



clean :
	rm -rf  ./lj
