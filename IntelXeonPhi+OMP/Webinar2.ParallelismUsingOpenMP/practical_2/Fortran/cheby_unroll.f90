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
!Unrolled loops implementation. The unrolling allows vectorization
!on the innermost loops because they now have no iteration dependencies -
!iterations are independent and therefore parallelisable with vector
!instructions. Unrolling also allows better reuse of the arrays indexed
!by k.
 
module cheby

use util, only : wp, pi
integer, parameter :: urf = 64

contains

subroutine evaluate_Tpoly(fe, xbar, nx, coeffs, nplus1)
   implicit none
   integer, intent(in) :: nx, nplus1
   real (kind=wp), intent(in) :: coeffs(nplus1), xbar(nx)
   real (kind=wp), intent(out) :: fe(nx)
   integer :: i, n, k, nplus2, j
   real (kind=wp) :: factor(urf), dk(urf), bk(urf)

   n = nplus1 - 1
   nplus2 = nplus1 + 1

   do i = 1, nx-(urf-1), urf
      do j = 1, urf
         factor(j) = 2.0_wp*(1.0_wp + xbar(i+j-1))
         dk(j) = 0.0_wp
         bk(j) = 0.0_wp
      end do
      do k = nplus1, 2, -1
         do j = 1, urf
            dk(j) = coeffs(k) - dk(j) + factor(j)*bk(j)
            bk(j) = dk(j) - bk(j)
         end do
      end do
      do j = 1, urf
         fe(i+j-1) = 0.5_wp*coeffs(1) - dk(j) + 0.5_wp*factor(j)*bk(j)
      end do
   end do
   j = urf*(nx/urf) + 1;
   do i = j, nx
      factor(1) = 2.0_wp*(1.0_wp + xbar(i))
      dk(1) = 0.0_wp
      bk(1) = 0.0_wp
      do k = nplus1, 2, -1
         dk(1) = coeffs(k) - dk(1) + factor(1)*bk(1)
         bk(1) = dk(1) - bk(1)
      end do
      fe(i) = 0.5_wp*coeffs(1) - dk(1) + 0.5_wp*factor(1)*bk(1)
   end do

end subroutine 



subroutine find_Tcoeffs(coeffs, fr, nplus1)
   implicit none
   integer, intent(in) :: nplus1
   real (kind=wp), intent(in) :: fr(nplus1)
   real (kind=wp), intent(out) :: coeffs(nplus1)
   integer :: i, k, n, j
   real (kind=wp) :: fln, piby2n, f0, halffn, fli(urf), factor(urf), bk(urf), dk(urf)

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

   do i = 1, nplus1-(urf-1), urf
      do j = 1, urf
         fli(j) = real(i-1+j-1,kind=wp)
         factor(j) = 4.0_wp*sin(piby2n*fli(j))**2
         bk(j) = halffn
         dk(j) = halffn
      end do
      do k = n, 2, -1
         do j = 1, urf
            dk(j) = fr(k) + dk(j) - factor(j)*bk(j)
            bk(j) = bk(j) + dk(j)
         end do
      end do
      do j = 1, urf
         coeffs(i+j-1) = (f0+2.0E0_wp*dk(j)-factor(j)*bk(j))/fln
      end do
   end do
   j = urf*(nplus1/urf) + 1
   do i = j, nplus1
      fli(1) = real(i-1,kind=wp)
      factor(1) = 4.0_wp*sin(piby2n*fli(1))**2
      bk(1) = halffn
      dk(1) = halffn
      do k = n, 2, -1
         dk(1) = fr(k) + dk(1) - factor(1)*bk(1)
         bk(1) = bk(1) + dk(1)
      end do
      coeffs(i) = (f0+2.0E0_wp*dk(1)-factor(1)*bk(1))/fln
   end do

   coeffs(nplus1) = 0.5_wp*coeffs(nplus1)

end subroutine find_Tcoeffs



subroutine gen_data(fr, xbar, n1)
   implicit none
   integer, intent(in) :: n1
   real (kind=wp), intent(out) :: fr(n1)
   real (kind=wp), intent(out) :: xbar(n1)
   integer :: n, i
   n = n1 - 1
   do i = 1, n1
      xbar(i) = cos(pi*real(i-1,kind=wp)/real(n,kind=wp))
      fr(i) = exp(xbar(i))
   end do
end subroutine gen_data

end module cheby
