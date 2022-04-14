module fdutils
  implicit none
  private :: order
  private :: half
  public :: get_stencil
  public :: coeffs
  public :: get_coeffs
  public :: fd12
  public :: fill_input
  public :: exact_laplacian
  public :: error
  real, dimension(13,13) :: coeffs
  integer, parameter :: order = 12
  integer, parameter :: half = 6
  contains
    integer function get_stencil(n,i)
      implicit none
      integer, intent(in) :: n 
      integer, intent(in) :: i
      integer             :: j
      integer             :: m
      integer             :: tmp

      j=i-1
      if (j<half) then
        get_stencil=j
        return
      end if
      if (j>=n-half) then
        tmp=n-half
        get_stencil=(half+(j-tmp)+1)
        return
      end if
      get_stencil=half
      return
    end function get_stencil

    subroutine get_coeffs
      call fill_coeffs(coeffs)
    end subroutine get_coeffs

    subroutine fd12(in_u,out_u,dx,dy,dz,nx,ny,nz)
      real, dimension(:,:,:), allocatable, intent(in) :: in_u
      real, dimension(:,:,:), allocatable             :: out_u
      real, intent(in) :: dx
      real, intent(in) :: dy
      real, intent(in) :: dz
      integer, intent(in) :: nx
      integer, intent(in) :: ny
      integer, intent(in) :: nz



      integer ix,iy,iz,i
      integer jx,jy,jz;
      real d2fdx2,d2fdy2,d2fdz2

      do ix=1,nx
        do iy=1,ny
          do iz=1,nz
            jx=get_stencil(nx,ix)+1
            jy=get_stencil(ny,iy)+1
            jz=get_stencil(nz,iz)+1
            d2fdx2=0.0
            d2fdy2=0.0
            d2fdz2=0.0
            do i=-jx+1,order-jx+1
              d2fdx2 = d2fdx2 + coeffs(i+jx,jx)*in_u(ix+i,iy,iz)/(dx*dx)
            enddo

            do i=-jy+1,order-jy+1
              d2fdy2 = d2fdy2 + coeffs(i+jy,jy)*in_u(ix,iy+i,iz)/(dy*dy)
            enddo

            do i=-jz+1,order-jz+1
              d2fdz2 = d2fdz2 + coeffs(i+jz,jz)*in_u(ix,iy,iz+i)/(dz*dz)
            enddo

            out_u(ix,iy,iz)=d2fdx2+d2fdy2+d2fdz2

          enddo
        enddo
      enddo
    end subroutine fd12



    subroutine exact_laplacian(out_u,x0,y0,z0,dx,dy,dz,nx,ny,nz)
    real, dimension(:,:,:), allocatable :: out_u
    real, intent(in) :: x0
    real, intent(in) :: y0
    real, intent(in) :: z0
    real, intent(in) :: dx
    real, intent(in) :: dy
    real, intent(in) :: dz
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: nz

    integer ix,iy,iz;
    real x,y,z;
    real d1,d2,d3;

    do ix=1,nx
      do iy=1,ny
        do iz=1,nz
          x=x0+ix*dx
          y=y0+iy*dy
          z=z0+iz*dz
          d1=-sin(x+y+z)
          d2=-sin(x+y+z)
          d3=-sin(x+y+z)
          out_u(ix,iy,iz)=d1+d2+d3
        enddo
      enddo
    enddo
    end subroutine exact_laplacian


    subroutine fill_input(out_u,x0,y0,z0,dx,dy,dz,nx,ny,nz)
    real, dimension(:,:,:), allocatable :: out_u
    real, intent(in) :: x0
    real, intent(in) :: y0
    real, intent(in) :: z0
    real, intent(in) :: dx
    real, intent(in) :: dy
    real, intent(in) :: dz
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: nz


    integer ix,iy,iz;
    real x,y,z;

    do ix=1,nx
      do iy=1,ny
        do iz=1,nz
          x=x0+ix*dx
          y=y0+iy*dy
          z=z0+iz*dz
          out_u(ix,iy,iz)=sin(x+y+z)
        enddo
      enddo
    enddo


    end subroutine fill_input


    real function error(exact,approx)
      real,dimension(:,:,:),intent(in) :: exact
      real,dimension(:,:,:),intent(in) :: approx

      integer i,j,k;
      real merr
      real tmp
      real nrm

      merr = 0.0
      nrm  = 0.0
      do i=1,size(exact,1)
        do j=1,size(exact,2)
          do k=1,size(exact,3)
            nrm=nrm+exact(i,j,k)*exact(i,j,k)
            tmp=exact(i,j,k)-approx(i,j,k)
            merr=merr+tmp*tmp
          enddo
        enddo
      enddo

      

      error=sqrt(merr)/sqrt(nrm)
      return
    end function error

end module fdutils

