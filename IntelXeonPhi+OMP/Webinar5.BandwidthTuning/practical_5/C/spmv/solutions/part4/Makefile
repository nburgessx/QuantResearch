

CC := icpc
CFLAGS := -fno-alias -xMIC-AVX512 -O2 -qopenmp --std=c++11 -I /fserver/reida/clones/memkind/include/ 

main : main.cpp sparse.h
	$(CC) $(CFLAGS) main.cpp -o main /usr/lib64/libmemkind.so.0.0.1




.PHONY : clean




clean :
	rm -rf ./main
