program hello_world_solution
  use omp_lib
  integer :: id
!$OMP PARALLEL num_threads(4) private(id)
  id = omp_get_thread_num()
  print *,'Hello from thread ',id
!$OMP END PARALLEL
end program hello_world_solution
