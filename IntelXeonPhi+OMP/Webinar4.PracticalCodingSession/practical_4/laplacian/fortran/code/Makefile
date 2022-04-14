



FORTRAN := gfortran -g
CC := gcc

fd : fd.f90 util.o fdutils.o coeffs.o
	$(FORTRAN) fd.f90 util.o coeffs.o fdutils.o -ofd

fdutils.o : fdutils.f90
	$(FORTRAN) -c fdutils.f90

util.o : util.f90
	$(FORTRAN) -c util.f90 

coeffs.o : coeffs.c coeffs.h
	$(CC) -c coeffs.c 




.PHONY : clean




clean :
	rm -rf fdutils.o
	rm -rf util.o
	rm -rf fd
	rm -rf fdutils.mod
	rm -rf util.mod
	rm -rf coeffs.o
