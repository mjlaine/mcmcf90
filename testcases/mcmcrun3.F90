!!! ------------------------------------------------------------------------
!!! Example that uses its own initialize n:o 3
!!!
!!! Marko Laine <marko.laine@fmi.fi>
!!! ------------------------------------------------------------------------

!!! Defines main, ssfunction, modelfunction, checkbounds

!!!
!!! Main program
!!! 
program mcmcmain
  use mcmcprec
  use mcmcmod, only : MCMC_setpar0, MCMC_setsigma2nobs, MCMC_setcmat0
  implicit none

  !! Set initial values for par0, cmat0, sigma2, nobs
  call MCMC_setpar0((/10.0_dbl,0.2_dbl/))
  call MCMC_setpar0('mcmcpar.dat') ! read file
  call MCMC_setcmat0(0.1_dbl) ! pct
  call MCMC_setcmat0((/0.2_dbl,0.01_dbl/)) ! or as std
  call MCMC_setcmat0('mcmccov.dat')
  call MCMC_setsigma2nobs(0.5_dbl,11) ! sigma2, nobs

  call mcmc_main()   ! call mcmc_main in the library libmcmcrun.a

end program mcmcmain

!!!
!!! -2*log(p(obs|params)) for
!!! general Gaussian likelihood, data are read from file
!!!
function ssfunction(theta,npar,ny) result(ss)
  use mcmcprec
  use matutils, only : loaddata

  implicit none
  integer(kind=ik4) :: npar, ny
  real(kind=dbl) :: theta(npar)
  real(kind=dbl) :: ss(ny)

  real(kind=dbl), save, pointer :: data(:,:)
  logical, save :: first = .true.

  interface
     function modelfunction(theta,x) result(y)
       use mcmcprec
       implicit none
       real(kind=dbl) :: theta(:), x(:), y(size(x))
     end function modelfunction
  end interface

  !! during the first call, data is loaded
  if (first) then
     call loaddata('data.dat',data)
     first = .false.
  endif
 
  !! calculate the actual sum of squares
  ss(1) = sum((data(:,2) - modelfunction(theta,data(:,1)))**2)
end function ssfunction

!!!
!!! model function used by ssfunction
!!!
function modelfunction(theta,x) result(y)
  use mcmcprec
  implicit none
  real(kind=dbl) :: theta(:), x(:), y(size(x))
  if (size(theta) .ne. 2) stop 'something wrong with sizes in modelfunction'
  !! simple exponential model
  y = theta(1)*exp(-theta(2)*x)
end function modelfunction

!!!
!!! this function returns false if any theta(i) is out of bounds
!!!
function checkbounds(theta)
  implicit none
  real*8 theta(:)
  logical checkbounds  
!! example: all thetas must be positive
  checkbounds = .true.
  if (any(theta<=0.0)) checkbounds = .false.
  return
end function checkbounds

