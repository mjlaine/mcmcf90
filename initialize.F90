!!! $Id: initialize.F90,v 1.4 2012/07/03 10:57:34 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mcmc library
!!! File: initialize.F90
!!! Purpose:
!!!
!!!  Initialize some variables for the MCMC run. This is called from
!!!  MCMC_init after the namelist mcmcinit.nml has bee read.
!!!  Allocation of the vectors shoud be done here, also.
!!! 
!!!  You can replace this by providing your own initialize subroutine.
!!!
!!! Marko Laine 2012 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!!

!!! Set initial values for:
!!!  par0, npar, cmat0, initcmatn, sigma2, nobs, nycol
!!!  allocate  par0(:), cmat0(:,:), sigma2(:) and nobs(:)
subroutine initialize(par0,npar,cmat0,initcmatn,sigma2,nobs,nycol)

  use mcmcprec
  use matutils, only : loaddata, readdata
  use mcmcinit, only : cov0file, covnfile, parfile, sigma2file
  use mcmcmod, only: verbosity

  implicit none
  
  integer, intent(inout) :: npar, initcmatn, nycol
  integer, intent(inout), allocatable :: nobs(:)
  real(kind=dbl), intent(inout), allocatable :: par0(:), cmat0(:,:)
  real(kind=dbl), intent(inout), allocatable :: sigma2(:)

  !! loadmat uses pointers
  real(kind=dbl), pointer :: cmat(:,:), xpar(:), s2n(:,:)
  real(kind=dbl), pointer :: xnycol(:,:)
  integer :: allocstat, fstat

  !! try to read nycol from file
  call loaddata('mcmcnycol.dat',xnycol,fstat)
  if (fstat /= 0) then
     if (verbosity>2) write(*,*) 'note: file mcmcnycol.dat not found, setting nycol = 1'
     nycol = 1
  else
     if (verbosity>0) write(*,*) 'note: read mcmcnycol.dat'
     if (verbosity>2) write(*,*) 'nycol=',nycol
     nycol = int(xnycol(1,1))
     deallocate(xnycol)
  end if

  allocate(sigma2(nycol), stat=allocstat)
  if (allocstat /= 0) stop 'ERROR: allocate sigma2'
  allocate(nobs(nycol), stat=allocstat)
  if (allocstat /= 0) stop 'ERROR: allocate nobs'

  call loaddata(parfile,xpar,fstat)
  if (fstat /= 0) then
     write(*,*) 'ERROR: Error reading file, ', trim(parfile)
     write(*,*) '       fstat = ', fstat
     write(*,*) 'STOPPING NOW'
     stop
  end if
  npar = size(xpar)
  if (npar < 1) then
     write(*,*) 'ERROR: npar = 0!!'
     stop
  end if
  if (verbosity>2) write(*,*) 'note: read ',trim(parfile),', npar =',npar
  !! allocate and copy xpar to par0
  allocate(par0(npar), stat=allocstat)
  if (allocstat /= 0) stop 'allocate'  
  par0 = xpar
  deallocate(xpar)


  call loaddata(cov0file,cmat,fstat,uselock=.true.)
  if (fstat /= 0 .or. size(cmat,1) /= npar .or. &
       size(cmat,2) /= npar ) then
     write(*,*) 'ERROR: Error reading file mcmccov.dat, status:',fstat
     write(*,*) '       npar = ',npar
     stop
  end if
  if (verbosity>2) write(*,*) 'note: read ', trim(cov0file)
  !! copy cmat to cmat0
  allocate(cmat0(npar,npar), stat=allocstat)
  if (allocstat /= 0) stop 'allocate'
  cmat0 = cmat
  deallocate(cmat)

  !! read mcmccovn.dat for initcmatn, this is also mcmcinit.nml variable!!???
  !! now also the name of the file is mcmcinit.nml variable
  if (len_trim(covnfile) > 0) then
     call readdata(covnfile,initcmatn,stat=fstat,uselock=.true.)
     if (fstat /=0) then
        if (verbosity>2) write(*,*) 'note: could not open: ', trim(covnfile)
        if (verbosity>2) write(*,*) 'note: initcmatn: ', initcmatn
     else
        if (verbosity>2) write(*,*) 'note: read initcmatn:', initcmatn
     end if
  end if

  ! mcmcsigma2.dat is 2*nycol, first row sigma2 vector, 2nd row nobs vector
  call loaddata(sigma2file,s2n,fstat)
  if (fstat /= 0) then
     sigma2 = 1.0e0_dbl
     nobs = 1
     if (verbosity>2) write(*,*) 'note: no mcmcsigma2.dat, using default values'
  else
     if ( size(s2n,1) /= 2 .and. size(s2n,2) /= nycol ) then
        write(*,*) 'sizes: ', size(s2n,1), size(s2n,2)
        stop 'ERROR: error in mcmcsigma2.dat (obs: new format 4.8.2006)'
     else
        sigma2 = s2n(1,1:nycol)
        nobs = int(s2n(2,1:nycol))
     end if
     deallocate(s2n)
     if (verbosity>2) write(*,*) 'note: ',trim(sigma2file),' read'
  end if
  
end subroutine initialize
