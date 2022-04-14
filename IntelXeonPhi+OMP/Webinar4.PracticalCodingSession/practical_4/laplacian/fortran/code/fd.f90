program FD

  use util, only : time
  use fdutils
  !, only : get_stencil, get_coeffs, coeffs,j
  !fd12,fill_input,exact_laplacian,error

  implicit none

  integer, parameter :: order = 12
  integer, parameter :: half = 6



  !Grid dimensions
  integer :: nx = 512
  integer :: ny = 512
  integer :: nz = 512
  real :: dx = 0.5
  real :: dy = 0.5
  real :: dz = 0.5
  real :: x0 = 0.0
  real :: y0 = 0.0
  real :: z0 = 0.0
  real :: merr
  real :: t1
  real :: t2


  real, dimension(:,:,:), allocatable :: u_in
  real, dimension(:,:,:), allocatable :: u_out
  real, dimension(:,:,:), allocatable :: u_exact
  allocate(u_in(nx,ny,nz))
  allocate(u_out(nx,ny,nz))
  allocate(u_exact(nx,ny,nz))




  !Initialize finite difference coefficients.
  call get_coeffs()
  call fill_input(u_in,x0,y0,z0,dx,dy,dz,nx,ny,nz)
  call exact_laplacian(u_exact,x0,y0,z0,dx,dy,dz,nx,ny,nz)
  t1 = time()
  call fd12(u_in,u_out,dx,dy,dz,nx,ny,nz)
  t2 = time()
  merr=error(u_exact,u_out)
  print *,' Time(s): ' ,t2-t1 
  print *,' Error:   ' ,merr




  deallocate(u_in)
  deallocate(u_out)
  deallocate(u_exact)

end program FD




