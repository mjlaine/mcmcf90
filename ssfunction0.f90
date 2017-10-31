!!! 
!!! dummy ssfunction for libmcmcrun
!!!
function ssfunction(theta,npar,ny) result(ss)
  implicit none
  integer*4 npar, ny
  real*8 theta(npar)
  real*8 ss(ny)

  ss = 0.0d0
  write(*,*) 'ERROR: you are using the built in ssfunction'
  write(*,*) 'please do your own and relink your program'
  write(*,*) 'stopping now'
  stop

end function ssfunction
