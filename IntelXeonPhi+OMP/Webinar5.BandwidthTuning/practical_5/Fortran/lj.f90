! A program to calculate the Lennard-Jones Forces between a set of water molecules 

program LJ

  use util, only : time

  implicit none

  ! constants for water molecule 
  double precision, parameter :: sigma = 0.3165555 ! nanometers 
  double precision, parameter :: epsilon = 0.6501696 ! kJ/mol 

  ! numbers of molecules 
  integer, parameter ::  n = 100000

  double precision dx, dy, dz, r, mag, t1, t2;
  double precision, dimension(n) :: x, y, z, fx, fy, fz;
  integer  i, j;
  
  do i = 1, n
     ! x, y, z are positions of water molecules in space (nanometers)
     call random_number(x(i))
     call random_number(y(i))
     call random_number(z(i))
     x(i) = 10.*x(i)
     y(i) = 10.*y(i)
     z(i) = 10.*z(i)
     fx(i) = 0.
     fy(i) = 0.
     fz(i) = 0.
  end do

  t1 = time()
  
  do i = 1, n
     do j = 1, n
        
        if (i .ne. j) then
           ! r_j - r_i  vector components of vector difference
           dx = x(j) - x(i)
           dy = y(j) - y(i)
           dz = z(j) - z(i)
           
           ! calculate distrance between particles i and j 
           r = (dx*dx + dy*dy + dz*dz)
           r = sqrt(r)
           
           ! force magnitude 
           mag = -24.*epsilon*((2.*sigma**12/r**13) - (sigma**6/r**7))
           ! force components 
           fx(i) = fx(i) + mag*dx
           fy(i) = fy(i) + mag*dy
           fz(i) = fz(i) + mag*dz
        end if
        
     end do
  end do

  t2 = time()

  print *,' TIME TO COMPUTE LJ FORCES ',t2-t1 

end program LJ
