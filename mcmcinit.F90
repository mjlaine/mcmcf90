!!! $Id: mcmcinit.F90,v 1.9 2012/07/03 10:57:34 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mcmc library
!!! File: mcmcinit.F90
!!! Purpose: handle mcmcinit.nml namelist parameters
!!!
!!! Marko Laine 2008 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!!
!!! 
module mcmcinit

  use mcmcprec

  implicit none

  public
  !! mcmc namelist variables
  integer, save :: nsimu      ! length of the chain
  integer, save :: doadapt    ! do we adapt?
  integer, save :: doburnin   ! do we have burn in time?
  integer, save :: adaptint   ! adaptation interval
  integer, save :: adapthist  ! history size
  integer, save :: initcmatn  ! n for initial cmat
  integer, save :: burnintime ! time to burn in
  integer, save :: badaptint  ! burn-in adapt interval
  integer, save :: adaptend   ! end adaptation
  integer, save :: greedy     ! greedy burn in adaptation
  real(kind=dbl), save :: scalelimit  ! when to scale during burn in
  real(kind=dbl), save :: scalefactor ! scale this much
  real(kind=dbl), save :: drscale  ! scale in dr
  integer, save :: svddim ! not used, use condmax instead
  real(kind=dbl), save :: condmax    ! use svd with condmax
  real(kind=dbl), save :: condmaxini ! modest jtj reqularization
  real(kind=dbl), save :: N0, S02 ! error variance prior
  integer, save :: sstype ! type of the likelihood
  real(kind=dbl), save :: sstrans ! transform in mdstlsqs
  integer, save :: filepars   ! read initial values from files
  integer, save :: printint   ! how often to give statistics
  integer, save :: updatesigma ! do we update sigma2
  integer, save :: usrfunlen ! how many extra values to calculate?
  integer, save :: dumpint ! dump interval
  character(len=128), save :: chainfile ! name of the chain file
  character(len=128), save :: s2file ! name of the s2chain file
  character(len=128), save :: ssfile ! name of the sschain file
  character(len=128), save :: priorsfile ! name of the priors file
  character(len=256), save :: cov0file ! name of the initial covmat file
  character(len=256), save :: covffile ! name of the final covmat file
  character(len=256), save :: covnfile ! file to hold the weight of covmat
  character(len=256), save :: meanfile ! file to the mean
  character(len=256), save :: nmlffile ! file to save namelist at the end

  character(len=256), save :: parfile ! file to read initial parameter from
  character(len=256), save :: parffile ! file to save final parameter to
  character(len=256), save :: sigma2file ! file to read sigma2 and nobs
  character(len=256), save :: sigma2ffile ! etc

  integer, save :: verbosity ! how much to print
  character(len=10), save :: method ! the method used (dram, scam)
  real(kind=dbl), save :: alphatarget ! alpha target in ram
  real(kind=dbl), save :: nuparam ! ram parameter

  !! some extra control flags
  character(len=128), save :: s2dmpfile ! ascii dump
  character(len=128), save :: chdmpfile ! ascii dump
!!! common control variables
  logical, save :: dodr       ! do use DR step?
  logical, save :: doscam     ! do scam
  integer, save :: usesvd     ! use svd instead of chol


!!! namelist to hold initialization and other parameters
  namelist /mcmc/ nsimu, doadapt, doburnin, adaptint, adapthist, &
       scalelimit, scalefactor, drscale, &
       badaptint, adaptend, initcmatn, N0, S02, &
       filepars, burnintime, greedy, printint, updatesigma, &
       usrfunlen, chainfile, s2file, ssfile, svddim, condmax, &
       cov0file, covffile, covnfile, meanfile, &
       nmlffile, parfile, parffile, sigma2file, sigma2ffile, &
       condmaxini, sstype, sstrans, dumpint, priorsfile, verbosity, &
       method, alphatarget, nuparam

contains

  subroutine read_mcmcinit_namelist(status)
    implicit none
    integer, intent(out) :: status
    integer :: fstat
    character(len=128) :: nmlfile ! name of the 'mcmcinit.nml' file
    logical :: fexist

    status = 0

    call MCMC_init_namelist() ! default values for namelist parameters

    !! open nml-file for options
    !! first check if mcmcinit.nml exists and use it
    !! if not then read the name of the nml file from mcmcnml.txt
    !! which is made by the program nmlio (2007 version)
    inquire(file='mcmcinit.nml',exist=fexist)
    if (fexist) then
       nmlfile = 'mcmcinit.nml'
    else
       inquire(file='mcmcnml.txt',exist=fexist)
       if (fexist) then
          open(unit=10, file='mcmcnml.txt', status='old', iostat=fstat)
          read(10,*) nmlfile
          close(10)
          nmlfile=trim(nmlfile)
          if (len_trim(nmlfile) == 0) then
             nmlfile = 'mcmcinit.nml'       
          end if
       else
          nmlfile = 'mcmcinit.nml'
       end if
    end if
!!! read some parameters from namelist file
    open(unit=10, file=nmlfile, status='old', iostat=fstat)
    if (fstat /= 0) then
       write(*,*) ''
       write(*,*) 'File ',trim(nmlfile),' not found, no MCMC run'
       status = -1
       return
       !      write(*,*) 'Error opening file mcmcinit.nml'
       !      write(*,*) 'Using default values'
       !      stop
    else
       read(10,nml=mcmc,iostat=fstat)
       if (fstat /= 0) then
          write(*,*) 'Error reading mcmc namelist from file ',trim(nmlfile)
          write(*,*) ' status:',fstat
          write(*,*) 'No mcmc run'
          close(unit=10)
          status = -2
          return
       end if
       close(unit=10)
    end if
    if (verbosity>0) write(*,*) 'note: using nmlfile ',trim(nmlfile)

  end subroutine read_mcmcinit_namelist

!!!
!!! save current values to namelist
!!!
  subroutine write_mcmcinit_namelist(status,file)
    implicit none
    integer, intent(out) :: status
    character(len=*), intent(in), optional :: file
    integer :: fstat

    character(len=256) :: nmlfile !! fixed length here, please fix me !!!

    if (present(file)) then
       nmlfile = trim(file)
    else
       nmlfile = 'mcmcinit.nml'
    end if
    status = 0
    ! write mcmcinit.nml, need to have delim= with GFORTRAN if
    ! namelist contains character variables
    open(unit=10, file=nmlfile, status='replace', delim='APOSTROPHE', iostat=fstat)
    if (fstat /= 0) then
       write(*,*) 'ERROR: File ',trim(nmlfile),' not opened for writing'
       status = -1
       return
    else
       write(10,nml=mcmc,iostat=fstat)
       if (fstat /= 0) then
          write(*,*) 'ERROR: Error writing parameters namelist to file ',trim(nmlfile)
          write(*,*) '   status:',fstat
          close(unit=10)
          status = -2
          return
       end if
       close(unit=10)
    end if
  end subroutine write_mcmcinit_namelist

!!!
!!! default values for the namelist contron variables
!!! 
  subroutine MCMC_init_namelist()
    implicit none

    nsimu       = 0 ! 10000
    doadapt     = 1     ! do we adapt?
    doburnin    = 0     ! do we have burn in?
    burnintime  = 0     ! simulations before (real) adaptation
    badaptint   = -1    !
    greedy      = 0     !
    scalelimit  = 0.05_dbl  !
    scalefactor = 2.5_dbl   !
    drscale     = 0.0_dbl   ! DR scale (0 = no DR)
    adaptint    = 100   ! interval between adaptations
    adapthist   = 0     ! adaptation history size, AP
    adaptend    = 0     ! end adaptation at this time (if > 0)
    initcmatn   = 0     ! n for init cmat
    N0          = 1.0_dbl   ! sigma2 prior
    S02         = 0.0_dbl   ! sigma2 prior
    filepars    = 1     ! read init data from files?
    printint    = 500
    updatesigma = 1     ! update sigma2
    usrfunlen   = 0
    dumpint     = 0
    chainfile   = 'chain.dat'
    s2file      = 's2chain.dat'
    ssfile      = 'sschain.dat'
    priorsfile  = '' ! empty or 'priors.dat'
    cov0file    = 'mcmccov.dat'
    covffile    = 'mcmccovf.dat'
    covnfile    = ''
    meanfile    = 'mcmcmean.dat'
    nmlffile    = ''
    parfile     = 'mcmcpar.dat'
    parffile    = 'mcmcparf.dat'
    sigma2file  = 'mcmcsigma2.dat'
    sigma2ffile = 'mcmcsigma2f.dat'
    svddim      = 0
    condmax     = 0.0_dbl
    condmaxini  = 1.0e15_dbl
    sstype      = 0     ! type of log(likelihood) in mdstlsqs
    sstrans     = -1.0_dbl  ! transformation in mdstlsqs
    verbosity   = 1
    method      = 'dram'
    alphatarget = 0.234_dbl
    nuparam     = 0.7_dbl

  end subroutine MCMC_init_namelist

!!!
!!! sanity check on mcmcinit parameters, and some global variables
!!!
  subroutine check_mcmcinit_parameters()
    implicit none

    !! nsimu is checked in the calling routine
!    if (nsimu < 1) then
!       return
!     end if
    !! adapthist, etc
    if (adapthist < 0) then
       adapthist = 0
    end if
    if (adaptint < 0) then
       adaptint = 0
       doadapt = 0
    end if
    if (burnintime < 0) then
       burnintime = 0
    end if
    if (badaptint <= 0) then
       badaptint = adaptint
    end if
    if (badaptint == 0) then
       doburnin = 0
    end if
    if (initcmatn < 0) then
       initcmatn = 0
    end if
    if (scalelimit < 0.0_dbl .or. scalelimit > 0.5_dbl) then
       write(*,*) 'ERROR: Scalelimit control variable should be between [0,0.5]'
       stop
    end if
    if (scalefactor < 0.0_dbl) then
       write(*,*) 'Scalefactor < 0?? Changed to 1'
       scalefactor = 1.0_dbl
    end if

    !! sstype, sstrans are not used in the default standalone library
    !! sstrans, if not set defaults to -1
    !! then try if (old option) sstype is set
    !! sstype = 0 
    if (sstype <= 0) then
       sstype  = 1
       sstrans = 1.0_dbl
    elseif (sstype == 1) then
       if (sstrans < 0.0_dbl) then
          !          write(*,*) 'negative sstrans, will set it to 1'
          sstrans = 1.0_dbl 
       endif
    elseif (sstype == 2) then
       !       write(*,*) 'Warning sstype = 2, not defined yet'
       !       write(*,*) 'will use log normal now'
       sstype = 1
       sstrans = 0.0_dbl  ! 
    elseif (sstype == 3) then ! Poisson likelihood
       sstrans = -999.0e0_dbl
       updatesigma = 0
       ! sigma2 = 0.0 ! not allocated yet !!!
       S02 = 1.0_dbl
       !    elseif (sstype == 4) then ! Binomial likelihood
       !       sstrans = -999.0d0 
       !       updatesigma = 0
       !       ! sigma2 = 0.0 ! not allocated yet !!!
       !       S02 = 1.0
    elseif (sstype == 5) then ! t with fixed df
       if (sstrans < 1.0_dbl) then
          write(*,*) 'ERROR: please set sstrans = df, when sstype = 5'
          stop
       endif
       if (updatesigma /= 0) then
          write(*,*) 'no sigma2 update with t type ss'
          updatesigma = 0
       endif
       !       updatesigma = 0 !??
    elseif (sstype == 6) then ! Laplace aka L1
       sstrans = 1.0_dbl
       if (updatesigma /= 0) then
          write(*,*) 'no sigma2 update with Laplace ss'
          updatesigma = 0
       endif
       S02 = sqrt(S02)
    else
       write(*,*) 'unknown sstype in nml?'
       sstype = 1
       sstrans = 1.0_dbl
    end if

    !! some global variables

    !! scam
    if (method=='scam') then
       doscam = .true.
       if (verbosity>0) write(*,*) 'note: scam'
       if (condmax <= 0.0_dbl) condmax = 1.0e15_dbl
       doburnin = 0 !!!! XXX
       dodr = .false. !!! XXX
       drscale = 0.0_dbl
    else
       doscam = .false.
    end if

    !! ram
    if (method=='ram') then
       drscale = 0.0_dbl
    end if

    !! DR
    if (drscale <= 0.0_dbl) then
       dodr = .false.
    else
       dodr = .true.
       if (verbosity>0) write(*,*) 'note: using DR, with scale = ', real(drscale)
    end if

    if (condmax > 0.0_dbl) then
       usesvd = 1
       if (verbosity>0) write(*,*) 'note: using svd with condmax = ',real(condmax)
    else
       usesvd = 0
    end if

    !! s2dmpfile is s2file // .dmp
    if (index(s2file,'.')>0) then
       s2dmpfile = s2file(1:(index(s2file,'.')-1))//'.dmp'
    else
       s2dmpfile = trim(s2file)//'.dmp'
    end if
    !! chdmpfile is chainfile // .eva
    if (index(chainfile,'.')>0) then
       chdmpfile = chainfile(1:(index(chainfile,'.')-1))//'.eva'
    else
       chdmpfile = trim(chainfile)//'.eva'
    end if

  end subroutine check_mcmcinit_parameters

end module mcmcinit
