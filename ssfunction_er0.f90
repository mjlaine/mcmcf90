!!! default early rejection version of ssfunction
!!! calls ssfunction without er
function ssfunction_er(theta,npar,ny,sscrit) result(ss)

  use mcmcmod, only : verbosity

  implicit none
  integer(kind=4), intent(in) :: npar, ny
  real(kind=8), intent(in) :: theta(npar), sscrit
  real(kind=8) :: ss(ny)

  logical, save :: once = .true.

  interface 
     function ssfunction(theta,npar,ny)
       integer(kind=4) :: npar, ny
       real(kind=8) theta(npar)
       real(kind=8) ssfunction(ny)
     end function ssfunction
  end interface

  if (once) then
     if (verbosity>0) write(*,*) 'note: the default er ssfunction, so no er for ss'
     once = .false.
  end if

  ss = ssfunction(theta,npar,ny)

end function ssfunction_er
