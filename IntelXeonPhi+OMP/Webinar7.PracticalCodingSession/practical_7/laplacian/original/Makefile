


fd : fd.cpp
	icpc -qopenmp -O2 --std=c++11 -xMIC-AVX512 -qopt-report=5 -qopt-report-phase=vec fd.cpp -o fd






.PHONY : clean



clean :
	rm -rf ./fd
	rm -rf ./*.optrpt
