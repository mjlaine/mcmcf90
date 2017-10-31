!!! testcode for mcmc library
!!! Gaussian standalone test
!!!

program mcmcmain
  use mcmcprec
  use mcmcmod, only : MCMC_setpar0, MCMC_setsigma2nobs, MCMC_setcmat0
  implicit none

  !! Set initial values for par0, cmat0, sigma2, nobs
  call MCMC_setpar0(5, 0.0_dbl)
  call MCMC_setcmat0(0.1_dbl) ! pct

  call mcmc_main()   ! call mcmc_main in the library libmcmcrun.a

end program mcmcmain


!!!
!!! general gaussian target, parameters are read from files
!!!
function ssfunction(theta,npar,ny) result(ss)

  use mcmcprec
  use matutils, only : loaddata

  implicit none
  integer*4 npar, ny
  real*8 theta(npar)
  real*8 ss(ny)

  real(kind=dbl), save, pointer :: mu(:), lam(:,:)
  logical, save :: first = .true.

  !! first call
  if (first) then
     call loaddata('mcmctest_mu.dat',mu) ! mean
     call loaddata('mcmctest_lam.dat',lam) ! inverse covariance matrux
     if (size(mu) /= npar) then
        write(*,*) 'error in sizes'
        stop
     end if
     first = .false.
  endif
 
  !! (theta-mu)'*lam*(theta-mu) 
  ss(1) = dot_product(matmul(lam,theta-mu),theta-mu)

end function ssfunction
