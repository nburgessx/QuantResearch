subroutine epsf(x,e)
 real :: x, e
 e = 0.5*epsilon(x)
end subroutine epsf

subroutine eps(x,e)
 double precision :: x, e
 e = 0.5d0*epsilon(x)
end subroutine eps
