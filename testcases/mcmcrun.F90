!!! ------------------------------------------------------------------------
!!! mcmc library
!!! File: mcmcrun.F90
!!! Purpose: example main program for the standalone MCMC code
!!!
!!! Marko Laine <marko.laine@fmi.fi>
!!! ------------------------------------------------------------------------
!!
!! This file provides
!!
!! 1. main program that calls mcmc_main in the mcmc library
!!
!! 2. subroutine ssfunction that calculates the minus two times the log
!! likelihood function -2*log(y|theta) aka sum-of-squares function
!!
!!  function ssfunction(theta,npar,ny) result(ss)
!!    integer*4 npar, ny
!!    real*8 theta(npar)
!!    real*8 ss(ny)
!!  end function ssfunction
!!
!! The function defined below load data from file "data.dat" and then
!! calculates sum((ydata-model(theta,xdata))**2). The model is defined
!! in function modelfunction.
!!
!! 3. checkbounds function to check for out of bounds mcmc proposals
!!
!!
!! In addtion to "data.dat" we need the following files
!!
!! 1. The run time mcmc parameters in fortran namelist "mcmcinit.nml"
!!
!! 2. The following ascii .dat files
!!
!!    mcmcpar.dat  - initial value of the parameter
!!    mcmccov.dat  - initial proposal covariance matrix
!!
!!   optionally also
!!
!!    mcmcsigma2.dat - Gaussian error variance and number of observations,
!!      optional, needed id updatesigma = 1 in mcmcinit.nml
!!    mcmcnycol.dat - optional, number of values ssfunction returns (usually=1)
!! 

!!!
!!! Main program
!!! 
program mcmcmain
  implicit none
  !! call mcmc_main in the library libmcmcrun.a
  call mcmc_main()
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
