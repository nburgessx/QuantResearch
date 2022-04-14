program MYSTREAM
!DEC$ ATTRIBUTES NOINLINE :: add

      use util, only : time
      implicit none

interface
  subroutine add(z,x,y,n)
  double precision, allocatable, dimension(:)             :: z
  double precision, allocatable, dimension(:), intent(in) :: x
  double precision, allocatable, dimension(:), intent(in) :: y
  integer(8), intent(in) :: n
  end subroutine add
  
  subroutine fill(z,x,y,n)
  double precision, allocatable, dimension(:) :: z
  double precision, allocatable, dimension(:) :: x
  double precision, allocatable, dimension(:) :: y
  integer(8), intent(in) :: n
  end subroutine fill
    
end interface



      integer(8), parameter :: n      = 500000000
      integer(8), parameter :: ntimes = 20
      integer(8), parameter :: bytes_per_unit = 8
      integer(8), parameter :: nreads = 2
      integer(8), parameter :: nwrites = 1
      double precision :: t1
      double precision :: t2
      double precision :: gbs
      
      integer(8)            :: i
      double precision, dimension(:), allocatable :: x
      double precision, dimension(:), allocatable :: y
      double precision, dimension(:), allocatable :: z


      allocate(x(n))
      allocate(y(n))
      allocate(z(n))


      call fill(z,x,y,n)

      !Call add first time to "warm up"
      call add(z,x,y,n)

      t1=time()
      do i=1,ntimes
        call add(z,x,y,n)      
      enddo
      t2=time()

      print *,'Effective GB: ',ntimes*(bytes_per_unit*(nreads+nwrites)*n*1e-09)
      print *,'Time(s)     : ',(t2-t1)
      print *,'Bandwidth   :', ntimes*(bytes_per_unit*(nreads+nwrites)*n*1e-09)/(t2-t1)



      deallocate(x)
      deallocate(y)
      deallocate(z)

end program MYSTREAM


subroutine add(z,x,y,n)
  implicit none
  double precision, allocatable, dimension(:)             :: z
  double precision, allocatable, dimension(:), intent(in) :: x
  double precision, allocatable, dimension(:), intent(in) :: y
  integer(8), intent(in) :: n
  integer(8) :: i

!$OMP PARALLEL DO num_threads(68)
  do i=1,n
    z(i)=x(i)+y(i)
  enddo
!$OMP END PARALLEL DO 



end subroutine add


subroutine fill(z,x,y,n)
  implicit none
  double precision, allocatable, dimension(:) :: z
  double precision, allocatable, dimension(:) :: x
  double precision, allocatable, dimension(:) :: y
  integer(8), intent(in) :: n
  integer(8) :: i

!$OMP PARALLEL DO num_threads(68)
  do i=1,n
    z(i)=0.0
    x(i)=1.0
    y(i)=1.0
  enddo
!$OMP END PARALLEL DO 



end subroutine fill
