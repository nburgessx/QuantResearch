


FORT := ifort
FLAGS := -O3 -xMIC-AVX512 -qopenmp

main : main.f90 util.o
	$(FORT) $(FLAGS) main.f90 -o main util.o


util.o : util.f90
	$(FORT) -c util.f90



.PHONY : clean



clean :
	rm -rf ./main
	rm -rf ./util.o
	rm -rf ./util.mod
