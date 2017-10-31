!!! $Id: mcmc_main.F90,v 1.19 2012/07/03 10:57:34 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mcmc library
!!! File: mcmc_main.F90
!!! Purpose: main program (as a subroutine) for standalone MCMC code
!!!
!!! Marko Laine 2008 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!!
!!!
subroutine mcmc_main()

  use mcmcmod
  implicit none

  integer(4) fstatus

  write(*,*) 'MCMC code version: ',Mcmc_Code_Version

  call MCMC_init()

  if ( nsimu <= 0 ) then
     write(*,*) 'nsimu <= 0 stopping'
     stop
  end if

  MCMC_running = 1
  if (method == 'scam') then
     call MCMC_run_scam()
  elseif (method == 'er') then
     call MCMC_run_er()
  elseif (method == 'ram') then
     call MCMC_run_ram()
  else
     call MCMC_run()
  end if
  MCMC_running = 0

  call MCMC_writechains(fstatus)

  call MCMC_cleanup(fstatus)

end subroutine mcmc_main

!!!
!!! For consecutive one time runs with MCMC_run1
!!!
subroutine mcmc_main_one()

  use mcmcmod

  implicit none

  integer(4) fstatus

  write(*,*) 'MCMC code version: ',Mcmc_Code_Version

  call MCMC_init()

  if ( nsimu <= 0 ) then
     stop 'nsimu <= 0'
  end if

  if (method == 'er') then
     call MCMC_run1_er()
  else
     call MCMC_run1()
  end if

  call MCMC_writechains(fstatus)

  call MCMC_cleanup(fstatus)

end subroutine mcmc_main_one
