!!! $Id: mcmcrun1.F90,v 1.6 2012/06/27 10:10:37 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mcmc library
!!! File: mcmcrun1.F90
!!! Purpose: module to save info between multiple calls for MCMC_run1
!!!
!!! Marko Laine 2010 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.

module mcmcrun1

  use mcmcprec ! floating point precision
  use matutils  ! utilities

  implicit none

  public
  real(kind=dbl) :: alpha12, sscrit
  integer :: drstage, isimu, nrej, ieval

  namelist /mcmcrun/ drstage, isimu, ieval, nrej, alpha12, sscrit

contains

!!! inital values for variables in  mcmcrun.nml namelist
  subroutine init_mcmcrun_namelist()
    implicit none

    drstage = 1
    isimu = 1
    nrej = 0
    ieval = 0
    alpha12 = 0.0_dbl
    sscrit = -1.0_dbl

  end subroutine init_mcmcrun_namelist

!!! read mcmcrun.nml namelist
  subroutine read_mcmcrun_namelist()
    implicit none
    integer :: fstat


    !! read mcmcnmlfile 
    open(unit=10, file='mcmcrun.nml', status='old', iostat=fstat)
    if (fstat /= 0) then
       write(*,*) 'ERROR: File mcmcrun.nml not found'
       stop
    else
       read(10,nml=mcmcrun,iostat=fstat)
       if (fstat /= 0) then
          write(*,*) 'ERROR: Error reading parameters namelist from file mcmcinit.nml'
          write(*,*) '   status:',fstat
          close(unit=10)
          stop
       end if
       close(unit=10)
    end if
  end subroutine read_mcmcrun_namelist


!!! save parameter values to mcmcrun.nml namelist 
  subroutine write_mcmcrun_namelist()
    implicit none
    integer :: fstat
    !! save the control variables to a namelist
    open(11, file='mcmcrun.nml', status='replace',iostat=fstat)
    write(11,nml=mcmcrun)
    close(11)

  end subroutine write_mcmcrun_namelist

end module mcmcrun1
