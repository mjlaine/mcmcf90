!!! $Id: mcmc.F90,v 1.40 2012/09/04 14:22:54 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mcmc library
!!! File: mcmc.F90 
!!! Purpose: the main mcmc module
!!!
!!!
!!! Marko Laine 2001 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!!
module mcmcmod

  use mcmcprec  ! floating point precision
  use matutils  ! utilities
  use matfiles  ! writing mat files
  use mcmcrand  ! random numbers
  use mcmcinit  ! mcmcinit.nml namelist variables

  implicit none

!!!
  public

#include "version.h"

  !! other global, but (some) private to the module, variables
  integer, public, save :: npar       ! n:o of parameters
  integer, private, dimension(:), allocatable, save :: nobs  ! n:o of observations
  integer, public, save :: nycol      ! n:o of y variables
  real(kind=dbl), public, dimension(:,:), allocatable, save :: chain
  real(kind=dbl), public, dimension(:,:), allocatable, save :: s2chain
  real(kind=dbl), public, dimension(:,:), allocatable, save :: sschain
  real(kind=dbl), private, dimension(:), allocatable, save :: par0 ! par removed
  real(kind=dbl), private, dimension(:,:), allocatable, save :: R, R2, iC
  real(kind=dbl), private, dimension(:), allocatable, save :: qcovstd
  real(kind=dbl), private, dimension(:,:), allocatable, save :: cmat0
  real(kind=dbl), private, dimension(:,:), allocatable, save :: chaincmat
  real(kind=dbl), private, dimension(:), allocatable, save :: chainmean
  real(kind=dbl), private, save :: chainwsum ! sum of weights, (change name !!!)
  real(kind=dbl), public, dimension(:), allocatable, save :: sigma2  ! error variance
  real(kind=dbl), private, dimension(:,:), allocatable, save :: usrfunmat

  integer, private, save :: inited = 0 ! have we initialized?
  integer, public, save :: MCMC_running = 0 ! are we doing the real thing yet
  integer, private, save :: stayed     ! 
  integer, private, save :: chainisstuck
  integer, private, save :: bndstayed  ! stayed because of bounds
  integer, private, save :: erstayed  ! stayed because of er
  integer, private, save :: ncolchain  ! columns in chain
  integer, public, save :: simuind    ! loop index
  integer, public, save :: chainind   ! indet to current chain index
  integer, private, save :: vecss      ! vectorized ss and s2

  integer, private, save :: draccepted, drtries ! DR

  logical, private, save :: nparok  = .false.
  logical, private, save :: par0ok  = .false.
  logical, private, save :: cmat0ok = .false.
  logical, private, save :: sigma2ok = .false.

  interface MCMC_setpar0
     module procedure MCMC_setpar0_vec, MCMC_setpar0_n, MCMC_setpar0_file
  end interface
  interface MCMC_setcmat0
     module procedure MCMC_setcmat0_mat, MCMC_setcmat0_pct, MCMC_setcmat0_std, &
          MCMC_setcmat0_file
  end interface
  interface MCMC_setsigma2nobs
     module procedure MCMC_setsigma2nobs_sca, MCMC_setsigma2nobs_vec
  end interface

contains

!!! the module subroutines are included from external files

#include "MCMC_init.F90"

#include "MCMC_run.F90"

#include "MCMC_run_er.F90"

#include "MCMC_run1.F90"

#include "MCMC_run1_er.F90"

#include "MCMC_run_scam.F90"

#include "MCMC_run_ram.F90"

#include "MCMC_adapt.F90"

#include "MCMC_DRAM.F90"

#include "MCMC_aux.F90"

#include "MCMC_dump.F90"

#include "MCMC_userfun.F90"

#include "MCMC_signal_handler.F90"

end module mcmcmod
