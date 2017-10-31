!!! $Id: MCMC_signal_handler.F90,v 1.6 2012/07/03 10:57:34 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mcmc library
!!! File: MCMC_signal_handler.F90 
!!! Purpose: signal handler modulde
!!!
!!! Marko Laine 2008 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!! 
!!! works with ifort, gfortran and VisualFortran, at least
!!!

subroutine signal_handler_init()

#ifdef _WIN32  
  use dflib     ! for signal handler
#endif

#ifdef __IFC
  use iflport   ! Intel fortran portability library
!!! you may have to copy and compile iflport.f90 in the current directory
  !  include 'iflport.f90'
#endif

  implicit none
#ifdef _WIN32  
  integer(4) iret
#endif
#ifdef __IFC
  integer(4) iret
#endif

!!! install signal handler for ctrl-c
#ifdef _WIN32
  iret = signalqq(SIG$INT, cc_handler)
#endif
#ifdef __IFC
  iret = signalqq(SIG$INT, cc_handler)
#endif

#if defined(__GFORTRAN__) || defined(__PGI)
  interface
     subroutine  signalqq(h)
       external h
     end subroutine signalqq
  end interface

  call signalqq(cc_handler)
#endif

end subroutine signal_handler_init

!!!
!!! control-c handler
!!!
#ifdef _WIN32
function cc_handler(signum)
  implicit none
  !dec$attributes c :: cc_handler
  integer(4) cc_handler
  integer(2) signum
  integer status

!!! save chain
  write(*,*) 'Saving chain upto ', simuind
  call MCMC_writechains(status)

  cc_handler = 1
  stop 'Coltrol-c interrupt'
end function cc_handler
#endif

#ifdef __IFC
function cc_handler(signum)
  implicit none
  !dec$attributes c :: cc_handler
  integer(4) cc_handler
  integer(2) signum
  integer status

  if (MCMC_running == 1) then
!!! save chain
     write(*,*) 'Saving chain upto ', simuind
     call MCMC_writechains(status)
  end if

  cc_handler = 1
  stop 'Coltrol-c interrupt'
end function cc_handler
#endif


#if defined(__GFORTRAN__) || defined(__PGI)
subroutine cc_handler()
  implicit none

  integer :: status

  if (MCMC_running == 1) then
     write(*,*) 'Saving chain upto ', simuind
     call MCMC_writechains(status)
  end if

  stop 'Coltrol-c interrupt'

end subroutine cc_handler
#endif

