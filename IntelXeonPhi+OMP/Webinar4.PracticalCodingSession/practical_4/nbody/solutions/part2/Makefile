

CC := icpc
CFLAGS := --std=c++11 -O3 -qopenmp -xMIC-AVX512

lj : lj.cpp
	$(CC) $(CFLAGS) lj.cpp -o lj



.PHONY : clean



clean :
	rm -rf  ./lj
