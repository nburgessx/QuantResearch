PANDOC := pandoc
PANDOC_FLAGS := -s -f markdown -t html



all : vectoradd_unopt vectoradd_unopt_for 

vectoradd_unopt : vectoradd.cpp
	icpc -no-vec -qopenmp --std=c++11 -xHOST -qopt-report=4 vectoradd.cpp -o vectoradd_unopt


vectoradd_unopt_for : vectoradd.f90
	ifort -no-vec -qopenmp -xHOST vectoradd.f90 -o vectoradd_unopt_for



.PHONY : clean




clean : 
	rm -rf ./vectoradd.optrpt
	rm -rf ./vectoradd_unopt
	rm -rf ./vectoradd_for
	rm -rf ./vectoradd_for.optrpt
	rm -rf ./vectoradd_unopt_for
	rm -rf ./instructions.html
