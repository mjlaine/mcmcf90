!!! $Id: MCMC_init.F90,v 1.29 2012/09/04 14:22:53 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mcmc library
!!! File: MCMC_init.F90
!!! Purpose: initialize the MCMC calculations, read namelists etc.
!!!
!!! Marko Laine 2008 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!!
!!! 

!!! contains: MCMC_init, MCMC_setpar0, MCMC_setcov0, MCMC_setsigma2
!!! if nsimu <= 0 then init is cancelled
!!!
subroutine MCMC_init()
  implicit none

  include 'external_inc.h'

  !! local variables
  integer :: allocstat, info, status

  if (inited == 1) then
     write(*,*) 'Warning(mcmcinit): allready inited'
     return
  end if

  vecss       = 1     ! vectorize SS and s2? (this must be 1 now)
!  npar        = 0     ! length of the parameter vector
!  nycol       = 0     ! number of y columns
  chainwsum   = 0.0e0_dbl

  call read_mcmcinit_namelist(status)
  if (status /= 0) call doerror('error in mcmcinit namelist')
  call check_mcmcinit_parameters()

  if (nsimu < 1) then
     return
  end if

  filepars = 1
  usrfunlen = 0

  if (nparok) then ! user has called MCMC_setpar0 or MCMC_setcov0
     if (verbosity>0) write(*,*) 'note: user init for par0'
     if (.not.par0ok) call doerror('user initialization error')
     if (.not.cmat0ok) then
        allocate(cmat0(npar,npar))
        call MCMC_initcmat0(cmat0)
        cmat0ok = .true.
     end if
     if (.not.sigma2ok) then
        nycol = 1
        allocate(sigma2(nycol))
        allocate(nobs(nycol))
        sigma2 = 1.0_dbl
        nobs = 1
        sigma2ok=.true.
     end if
  else
     !! Call a routine that returns initial values for
     !!  par0, npar, cmat0, initcmatn, sigma2, nobs, nycol.
     !! These are all global in mcmcmod or mcmcinit modules
     !! par0(:), cmat0(:), sigma2(:), nobs(:) need to be allocated
     call initialize(par0,npar,cmat0,initcmatn,sigma2,nobs,nycol)
     nparok = .true.
     par0ok = .true.
     cmat0ok = .true.
     sigma2ok = .true.
  end if

  !! debugging
  if (verbosity>5) then
     write(*,*) 'npar',npar,initcmatn,nobs
     call printmat(par0,header='par0',maxlines=10)
     call printmat(cmat0,header='cmat0',maxlines=10,maxcols=5)
  end if

  !! global R and R2 are used to save the Cholesky factor of cmat
  allocate(R(npar,npar), stat=allocstat)
  if (allocstat /= 0) stop 'allocate'
  if (dodr) then
     allocate(R2(npar,npar), stat=allocstat)
     if (allocstat /= 0) stop 'allocate'
     allocate(iC(npar,npar), stat=allocstat)
     if (allocstat /= 0) stop 'allocate'          
  end if
  if (doscam) then
     allocate(qcovstd(npar), stat=allocstat)
     if (allocstat /= 0) stop 'allocate'          
  end if

  allocate(chainmean(npar), stat=allocstat)
  if (allocstat /= 0) stop 'allocate'
  allocate(chaincmat(npar,npar), stat=allocstat)
  if (allocstat /= 0) stop 'allocate'

  !! save cmat to chaincmat and par0 to chainmean
  chaincmat = cmat0
  chainmean = par0
  chainwsum = dble(initcmatn)

  if (verbosity>0) write(*,*) 'note: sigma2 =',real(sigma2), &
                              ' nobs =',nobs,' nycol =',nycol

  !! Try to calculate Cholesky factor of the initial covariance
  if (verbosity>2) write(*,*) 'note: calling MCMC_calculate_R ...'
  call MCMC_calculate_R(cmat0,info)
  if (info .ne. 0) call doerror('could not factor the initial covariance')
  if (verbosity>2) write(*,*) 'note: ... called MCMC_calculate_R'

  !! prior for Gaussian error variance
  if (S02 <= 0.0_dbl) then
     S02 = sigma2(1)
  end if

  !! allocate chains
  if (updatesigma /= 0) then
     allocate(s2chain(nsimu,nycol), stat=allocstat)
     if (allocstat /= 0) stop 'allocate'
  end if
  if (vecss == 0) then ! not used
     nycol = 1
  end if
  ncolchain = npar+1 ! [par,stay]

  allocate(chain(nsimu,ncolchain), stat=allocstat)
  if (allocstat /= 0) stop 'allocate'

  allocate(sschain(nsimu,nycol+1), stat=allocstat)
  if (allocstat /= 0) stop 'allocate'

  if (usrfunlen>0) then
     allocate(usrfunmat(nsimu,usrfunlen), stat=allocstat)
     if (allocstat /= 0) stop 'allocate'
  end if

  !! initialize random number generator
  random_verbosity = verbosity
  call random_initialize()

  call signal_handler_init()

  !! MCMC_savechain saves the chain, there
  !! the initail value is accepted, and chainind is incremented
  chainind = 0                  ! the current row of the chain
  simuind  = 1                  ! global simuind

  stayed    = 0
  bndstayed = 0
  erstayed  = 0
  draccepted = 0
  drtries    = 0

  inited    = 1

end subroutine MCMC_init

!!!
!!! These try to be user friendly initialization routines
!!!
!!! MCMC_setpar0 for initial value
!!! MCMC_setcmat0 for initial proposal covariance
!!! MCMC_setsigma2nobs for sigma2 and nobs if doing Gibbs
!!!

subroutine MCMC_setpar0_vec(par)
  implicit none
  real(kind=dbl), intent(in) :: par(:)

  integer :: n

  n = size(par)
  if (nparok .and. npar.ne.n) call doerror('par0 and npar dont match')
  npar = n

  if (allocated(par0)) then
     if (size(par0) .ne. n) call doerror('par0 already set')
  else
     allocate(par0(npar))
  end if
  par0 = par

  nparok = .true.
  par0ok = .true.

end subroutine MCMC_setpar0_vec

subroutine MCMC_setpar0_n(n,par)
  implicit none
  integer, intent(in) :: n
  real(kind=dbl), intent(in) :: par

  if (nparok .and. npar.ne.n) call doerror('par0 and npar dont match')
  npar = n

  if (allocated(par0)) then
     if (size(par0) .ne. n) call doerror('par0 already set')
  else
     allocate(par0(npar))
  end if
  par0 = par

  nparok = .true.
  par0ok = .true.

end subroutine MCMC_setpar0_n

subroutine MCMC_setcmat0_mat(cmat)
  implicit none
  real(kind=dbl), intent(in) :: cmat(:,:)

  integer :: n

  n = size(cmat,1)
  if (nparok .and. npar.ne.n) call doerror('cmat0 and npar dont match')
  npar = n
  if (size(cmat,2).ne.n) call doerror('cmat already set')

  if (allocated(cmat0)) then
     if (size(cmat0,1) .ne. n .or. size(cmat0,2) .ne. n) &
          call doerror('initialization error with cmat0 size')
  else
     allocate(cmat0(npar,npar))
  end if

  cmat0 = cmat
  nparok = .true.
  cmat0ok = .true.

end subroutine MCMC_setcmat0_mat


subroutine MCMC_setcmat0_pct(pct)
  implicit none
  real(kind=dbl), intent(in) :: pct

  integer :: i

  if (.not.nparok.or..not.par0ok) call doerror('par0 not defined when calling setcmat0')

  if (allocated(cmat0)) then
     if (size(cmat0,1) .ne. npar .or. size(cmat0,2) .ne. npar) &
          call doerror('initialization error with cmat0 size')
  else
     allocate(cmat0(npar,npar))
  end if

  cmat0 = 0.0_dbl
  do i=1,npar
     if (par0(i) .ne. 0.0_dbl) then
        cmat0(i,i) = abs(pct*par0(i))**2
     else
        cmat0(i,i) = abs(pct)**2
     end if
  end do

  cmat0ok = .true.
end subroutine MCMC_setcmat0_pct

subroutine MCMC_setpar0_file(file)
  implicit none
  character(len=*), intent(in) :: file
  real(kind=dbl), allocatable :: par(:)
  
  call loaddata_veca(file,par)
  call MCMC_setpar0_vec(par)
  deallocate(par)
end subroutine MCMC_setpar0_file


subroutine MCMC_setcmat0_std(std)
  implicit none
  real(kind=dbl), intent(in) :: std(:)

  integer :: n, i

  n = size(std)

  if (.not.nparok) call doerror('npar not defined when calling setcmat0')
  if (n.ne.npar) call doerror('std size does not macth when calling setcmat0')

  if (allocated(cmat0)) then
     if (size(cmat0,1) .ne. npar .or. size(cmat0,2) .ne. npar) &
          call doerror('initialization error with cmat0 size')
  else
     allocate(cmat0(npar,npar))
  end if

  cmat0 = 0.0_dbl
  do i=1,npar
     cmat0(i,i) = std(i)**2
  end do

  cmat0ok = .true.
end subroutine MCMC_setcmat0_std

subroutine MCMC_setcmat0_file(file)
  implicit none
  character(len=*), intent(in) :: file
  real(kind=dbl), allocatable :: cmat(:,:)
  
  call loaddata_mata(file,cmat)
  call MCMC_setcmat0_mat(cmat)
  deallocate(cmat)
end subroutine MCMC_setcmat0_file

subroutine MCMC_setsigma2nobs_vec(sig2,n)
  implicit none

  real(kind=dbl), intent(in) :: sig2(:)
  integer, intent(in) :: n(:)

  integer :: l1, l2

  l1 = size(sig2)
  l2 = size(n)

  if (l1.ne.l2) call doerror('sig2 and nobs sizes')
  nycol = l1

  if (allocated(sigma2)) then
     if (size(sigma2).ne.nycol) call doerror('sigma2 size')
  else
     allocate(sigma2(nycol))
  end if

  if (allocated(nobs)) then
     if (size(nobs).ne.nycol) call doerror('nobs size')
  else
     allocate(nobs(nycol))
  end if

  sigma2 = sig2
  nobs = n

  sigma2ok = .true.
end subroutine MCMC_setsigma2nobs_vec


subroutine MCMC_setsigma2nobs_sca(sig2,n)
  implicit none
  real(kind=dbl), intent(in) :: sig2
  integer, intent(in) :: n
  call MCMC_setsigma2nobs_vec((/sig2/),(/n/))
end subroutine MCMC_setsigma2nobs_sca

!!!
subroutine MCMC_initcmat0(cmat)
  implicit none
  real(kind=dbl), intent(inout) :: cmat(:,:)
  integer :: i
  cmat = 0.0_dbl
  do i=1,size(cmat,1)
     cmat(i,i) = 1.0_dbl
  end do
end subroutine MCMC_initcmat0
