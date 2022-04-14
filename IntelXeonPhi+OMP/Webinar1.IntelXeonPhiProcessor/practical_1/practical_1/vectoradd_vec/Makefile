PANDOC := pandoc
PANDOC_FLAGS := -s -f markdown -t html



all: vectoradd vectoradd_for instructions.html


instructions.html : INSTRUCTIONS.md
	$(PANDOC) $(PANDOC_FLAGS) -o instructions.html INSTRUCTIONS.md



vectoradd : vectoradd_vec_solution.cpp
	icpc -qopenmp --std=c++11 -xHOST -qopt-report=4 vectoradd_vec_solution.cpp -o vectoradd

vectoradd_for : vectoradd_vec_solution.f90
	ifort -qopenmp -xHOST vectoradd_vec_solution.f90 -o vectoradd_for



.PHONY : clean




clean : 
	rm -rf ./vectoradd
	rm -rf ./vectoradd_vec_solution.optrpt
	rm -rf ./vectoradd_for
	rm -rf ./vectoradd_vec_solution_for.optrpt
	rm -rf ./instructions.html
