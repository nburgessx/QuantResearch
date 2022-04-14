!
!CWA - Feb 2014
!--------------
!Generate some (xbar, fr) values of a function, evaluate the
!coefficients of a Chebyshev polynomial which interpolates
!the generated values exactly at the given points (fe), evaluate
!the polynomial at the given points and then find the Euclidean
!norm of the difference between the original fr and the evaluated
!fr.
!
!OpenMP implementation with orphaned directives, for 1pr variant.

module cheby

use util, only : wp, pi

contains

subroutine evaluate_Tpoly(fe, xbar, nx, coeffs, nplus1)
   implicit none
   integer, intent(in) :: nx, nplus1
   real (kind=wp), intent(in) :: coeffs(nplus1), xbar(nx)
   real (kind=wp), intent(out) :: fe(nx)
   integer :: i, n, k, nplus2
   real (kind=wp) :: factor, dk, bk

   n = nplus1 - 1
   nplus2 = nplus1 + 1

   !$omp do schedule(static)
   do i = 1, nx
      factor = 2.0_wp*(1.0_wp + xbar(i))
      dk = 0.0_wp
      bk = 0.0_wp
      do k = nplus1, 2, -1
         dk = coeffs(k) - dk + factor*bk
         bk = dk - bk
      end do
      fe(i) = 0.5_wp*coeffs(1) - dk + 0.5_wp*factor*bk
   end do
   !$omp end do

end subroutine 



subroutine find_Tcoeffs(coeffs, fr, nplus1)
   implicit none
   integer, intent(in) :: nplus1
   real (kind=wp), intent(in) :: fr(nplus1)
   real (kind=wp), intent(out) :: coeffs(nplus1)
   integer :: i, k, n
   real (kind=wp) :: fln, piby2n, f0, halffn, fli, factor, bk, dk

   if (nplus1<2) return
   if (nplus1==2) then
      coeffs(1) = fr(1) +fr(2)
      coeffs(2) = 0.5_wp*(fr(1) - fr(2))
      return
   end if

   n = nplus1 - 1
   fln = real(n,kind=wp)
   piby2n = 0.5_wp*pi/fln
   f0 = fr(1)
   halffn = 0.5_wp*fr(nplus1)

   !$omp do schedule(static)
   do i = 1, nplus1
      fli = real(i-1,kind=wp)
      factor = 4.0_wp*sin(piby2n*fli)**2
      bk = halffn
      dk = halffn
      do k = n, 2, -1
         dk = fr(k) + dk - factor*bk
         bk = bk + dk
      end do
      coeffs(i) = (f0+2.0E0_wp*dk-factor*bk)/fln
   end do
   !$omp end do
   coeffs(nplus1) = 0.5_wp*coeffs(nplus1)

end subroutine find_Tcoeffs



subroutine gen_data(fr, xbar, n1)
   implicit none
   integer, intent(in) :: n1
   real (kind=wp), intent(out) :: fr(n1)
   real (kind=wp), intent(out) :: xbar(n1)
   integer :: n, i
   n = n1 - 1
   !$omp do schedule(static)
   do i = 1, n1
      xbar(i) = cos(pi*real(i-1,kind=wp)/real(n,kind=wp))
      fr(i) = exp(xbar(i))
   end do
   !$omp end do
end subroutine gen_data

end module cheby
