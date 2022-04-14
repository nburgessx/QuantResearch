

all: vectoradd vectoradd_for 


vectoradd : vectoradd_parvec_solution.cpp
	icpc -qopenmp --std=c++11 -xHOST -qopt-report=4 vectoradd_parvec_solution.cpp -o vectoradd

vectoradd_for : vectoradd_parvec_solution.f90
	ifort -qopenmp -xHOST vectoradd_parvec_solution.f90 -o vectoradd_for



.PHONY : clean




clean : 
	rm -rf ./vectoradd
	rm -rf ./vectoradd_parvec_solution.optrpt
	rm -rf ./vectoradd_for
	rm -rf ./vectoradd_parvec_solution_for.optrpt
	rm -rf ./instructions.html
