module util

integer, parameter :: wp = kind(0.d0)
real (kind=wp), parameter :: pi = 3.14159265358979323846_wp
logical, parameter :: debug = .false.

contains

function time()
!$ use omp_lib
   double precision :: time
   integer, dimension(8) :: current_time
!$ time = omp_get_wtime()
!$ return
    call date_and_time(values=current_time)
    time = 60.d0*60.d0*current_time(5) + 60.d0*current_time(6) + &
               current_time(7) + current_time(8)/1000.d0
    return
end function time

end module util
