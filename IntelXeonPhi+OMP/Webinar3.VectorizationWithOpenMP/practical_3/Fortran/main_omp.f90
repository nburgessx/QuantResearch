program main

use cheby, only : gen_data, find_Tcoeffs, evaluate_Tpoly
use util, only : wp, debug, time
implicit none

integer :: n1, j
real (kind=wp), allocatable :: xbar(:), fr(:), coeffs(:), fe(:)
real (kind=wp) :: dif, var, nrm2, tol, macheps
double precision :: t1, t2

write(*,*) "Enter a value of n1 (>1)"
read(5,*) n1
if (n1 < 2) then
   write(0,*) "Error, n1 must be > 1."
   stop
end if

! Allocate space for the generated function and
! x-values
allocate(xbar(n1), fr(n1))

! Generate some data from a function
call gen_data(fr, xbar, n1)

write(*,*) "Input size = ",n1

if (debug) then
   write(*,*) "Input data (x, f(x)):"
   do j = 1, n1
      write(*,*) j, xbar(j), fr(j)
   end do
end if

! Allocate space for the evaluated function and
! the coefficients
allocate(fe(n1), coeffs(n1))

! Find the Chebyshev coefficients
t1 = time()
call find_Tcoeffs(coeffs, fr, n1)
t2 = time()
write(*,'(A,F10.6,A)') "Time to find coefficients = ",t2-t1," secs."

if (debug) then
   write(*,*) "Chebyshev coefficients:"
   do j = 1, n1
     write(*,*) coeffs(j)
   end do
end if

! Evaluate using full poly - should be exact
t1 = time()
call evaluate_Tpoly(fe, xbar, n1, coeffs, n1)
t2 = time()
write(*,'(A,F10.6,A)') "Time to evaluate poly = ",t2-t1," secs."

if (debug) write(*,*) "Original and evaluated values:"
var = 0.0_wp
!$omp parallel default(none) &
!$omp private(j,dif) &
!$omp shared(n1,fr,fe,var)
!$omp do schedule(static) reduction(+:var)
do j = 1, n1
   if (debug) write(*,*) fr(j), fe(j)
   dif = fr(j) - fe(j)
   var = var + dif*dif
end do
!$omp end do
!$omp end parallel

nrm2 = sqrt(var)
macheps = 0.5_wp*epsilon(var)
tol = macheps*real(n1,kind=wp)**2
if (nrm2 <= tol) then
   write(*,'(A,E13.6,A)') "Result is OK (Euclidean norm :",nrm2,")"
else
   write(*,*) "Error:"
   write(*,'(A,E13.6,A,E13.6)') "      Euclidean norm ",nrm2," greater than tolerance ",tol
end if

deallocate(xbar,fe,fr,coeffs)

end program main
