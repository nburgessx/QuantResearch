



CC := icpc
CFLAGS := -O3 -qopenmp --std=c++11 -xMIC-AVX512



main : main.cpp
	$(CC) $(CFLAGS) main.cpp -o main -lmemkind


.PHONY : clean



clean :
	rm -rf ./main
