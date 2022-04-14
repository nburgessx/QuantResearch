program vectoradd
  use omp_lib 
  integer, parameter :: n = 500000000
  double precision, allocatable, dimension(:) :: x, y, z
  double precision :: t1, t2
  allocate(x(n))
  allocate(y(n))
  allocate(z(n))
  t1 = omp_get_wtime()
!$OMP PARALLEL DO num_threads(4)
  do i = 1, n
     z(i) = x(i)+y(i)
  end do
!$OMP END PARALLEL DO 
  t2 = omp_get_wtime()
  print *,' time to compute vector add = ',t2 - t1,' s '
  deallocate(x)
  deallocate(y)
  deallocate(z)
end program vectoradd


