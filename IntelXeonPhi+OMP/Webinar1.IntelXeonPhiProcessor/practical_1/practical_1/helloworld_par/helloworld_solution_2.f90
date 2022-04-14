program hellworld_solution_2
  use omp_lib
  integer :: id
!$OMP PARALLEL num_threads(4) private(id)
!$OMP CRITICAL
    id=omp_get_thread_num()
    print *,'Hello from thread ',id
!$OMP END CRITICAL
!$OMP END PARALLEL 
end program hellworld_solution_2
