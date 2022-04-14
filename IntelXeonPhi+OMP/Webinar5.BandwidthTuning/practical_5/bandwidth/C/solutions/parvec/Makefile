



CC := icpc
CFLAGS := -O3 -qopenmp --std=c++11 -xMIC-AVX512 -fno-alias



main : main.cpp
	$(CC) $(CFLAGS) main.cpp -o main


.PHONY : clean



clean :
	rm -rf ./main
