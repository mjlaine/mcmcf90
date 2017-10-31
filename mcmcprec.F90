!!! $Id: mcmcprec.F90,v 1.15 2012/06/27 10:10:37 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! this file is part of mcmc library
!!! File: mcmcprec.F90
!!! Purpose: machine constants for the mcmc module
!!!
!!! Marko Laine 2001 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!!
module mcmcprec

!  use, intrinsic :: ieee_arithmetic ! only in f2003
  implicit none

  private
  public :: Byte, Short, Long 
  public :: Single, Double 
  public :: ik, rk 
  ! Integer kinds 
  integer, parameter :: Byte  = selected_int_kind(1)
  integer, parameter :: Short = selected_int_kind(4)
  integer, parameter :: Long  = selected_int_kind(8)
  ! Floating point kinds 
  integer, parameter :: Single = selected_real_kind(6)
  integer, parameter :: Double = selected_real_kind(15)
  ! Generic kinds 
  integer, parameter :: ik = Short  ! Generic integer 
  integer, parameter :: rk = Double ! Generic float

!! need specific real*8 and integer*4 for old fortran code compatibility
  integer, parameter, public :: dbl=Double
  integer, parameter, public :: ik4=Long
  !! log(realmin) is  -708.396418532264  = log(tiny(0.0d0))
#if defined(__GFORTRAN__)
  real(kind=dbl), parameter, public :: log_realmin = log(tiny(0.0e0_dbl))
#elif defined(__PGI)
  real(kind=dbl), parameter, public :: log_realmin = Z'c086232bdd7abcd2'
#else
  real(kind=dbl), parameter, public :: log_realmin = log(tiny(0.0e0_dbl))
#endif

end module mcmcprec
