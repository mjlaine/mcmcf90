!!! $Id: MCMC_run1_er.F90,v 1.6 2012/07/03 10:57:34 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mcmc library
!!! File: MCMC_run1_er.F90
!!! Purpose: MCMC loop done once, with early rejection
!!!
!!! Marko Laine 2010 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!!
!!!
!!! here we assume that: mcmcparf.dat has the last accepted chain value
!!! (MCMC_init puts it in global variable par), 
!!! mcmcpar.dat has the new value to be tried
!!! chain(1,:) par, nrej 
!!!    par either thr old value from previous run or the new accepted
!!!    rej how many times we have tried to propose a new values before one inside the bounds is accepted
!!!
!!! this saves proposed new value in mcmcparnew.dat
!!!
!!! mcmcrun.nml has 
!!! nrej = how many times
!!! alpha12 = last acceptance probability
!!! drstage = what stage we are in
!!! isimu = how many actual rows in the chain so far
!!! ieval = how calls to ssfunction so far
!!!
subroutine MCMC_run1_er()

  use mcmcrun1

  implicit none


  real(kind=dbl) :: ss(nycol)
  real(kind=dbl) :: newpar(npar)
  real(kind=dbl) :: sspri
  real(kind=dbl) :: alpha = 0.0_dbl
  logical :: reject, inbounds

  real(kind=dbl) :: oldpar1(npar)
  real(kind=dbl) :: ssprev1(nycol)
  real(kind=dbl) :: sspri1

  integer :: fstat

!! these are in module mcmcrun1
!  real(kind=dbl) :: alpha12, sscrit
!  integer :: drstage=1, isimu = 1, nrej=0, ieval = 0
!  namelist /mcmcrun/ drstage, isimu, ieval, nrej, alpha12, sscrit

  if (inited /= 1) call doerror('we have not inited')

  dodr = .false.
  drstage = 1

  if (updatesigma==1) then
     write(*,*) 'Warning: do not use updatesigma in ER!'
     updatesigma=0
  end if

  if (any(abs(sigma2-1.0e0_dbl)>1.0e-6_dbl)) then
     write(*,*) 'Warning: in ER sigma2 should be exactly 1'
     write(*,*) sigma2,abs(sigma2-1.0e0_dbl)
     sigma2=1.0e0_dbl
  end if


  !! open namelist to get the values from previous run
  call read_mcmcrun_namelist()

  if (isimu == 1) then
     if (verbosity>=0) write(*,*) 'first simulation'
     if (verbosity>=0) write(*,*) 'Early rejection version'

     newpar  = par0
     oldpar1 = par0
     sspri   = MCMC_priorfun(newpar)
     ss      = MCMC_ssfunction(newpar)
     ieval = ieval + 1
     if (verbosity>=0) write(*,*) 'ss1 = ',real(ss), ' sspri1 = ', real(sspri)

     ssprev1 = ss

     isimu = isimu + 1
     nrej  = 1
     reject = .false.

     !! should we read meanfile for the mean if we are continuing from
     !! previous run
     call readdata(meanfile,par0,stat=fstat,uselock=.true.)
     if (fstat /=0) then
        if (verbosity>4) write(*,*) 'note: no file ',trim(meanfile),' using par0:',real(par0)
!        par0 = par
     else
        if (verbosity>0) write(*,*) 'note: read ',trim(meanfile),' for the chain mean'
     end if

  else

     !! must read previous accepted point from mcmcparf.dat !!!
     call readdata('mcmcparf.dat',oldpar1)
     call readdata('mcmcssprev1.dat',ssprev1)

     !! par0 has the intial value which should now have the previous chainmean
     !! need to fix MCMC_adapt???
     call readdata(meanfile,par0,uselock=.true.)

     simuind = isimu              ! global simuind ???

     newpar = par0 ! from 'mcmcpar.dat'

     !! evaluate the latest ss and prior
     sspri = MCMC_priorfun(newpar)
     ss  = MCMC_ssfunction(newpar)
     ieval = ieval + 1

     if (verbosity>=0) write(*,*) 'isimu:',isimu
     if (verbosity>=0) write(*,*) 'oldpar1:',real(oldpar1),' ssprev1:',real(ssprev1),' sscrit',real(sscrit)

     !! we will need to recalculate the old priors for now
     sspri1 = MCMC_priorfun(ssprev1)
     if (verbosity>=3) write(*,*) ' sspri1:',real(sspri1)

     !! do we accept it?
     !! now using early rejection and sscrit
     if (sum(ss/sigma2)>=sscrit) then
        reject = .true.
        if (verbosity>=1) write(*,*) 'ss>=sscrit, rejecting, ',real(sum(ss/sigma2)),real(sscrit)
     else
        reject = .false.
        if (verbosity>=1) write(*,*) 'ss<sscrit, accepting, ',real(sum(ss/sigma2)),real(sscrit)
     end if

!     alpha  = MCMC_alpha(oldpar1, ssprev1, sspri1, newpar, ss, sspri)
!     reject = MCMC_reject(alpha)

     if (reject) then
        isimu = isimu + 1
        drstage = 1
     else ! accept
        drstage = 1
        isimu = isimu + 1
        nrej = 1
        oldpar1 = newpar
        ssprev1 = ss
        alpha12 = alpha
        sspri1  = sspri

     end if

     if (verbosity>=0) write(*,*) 'newpar:',real(newpar),'ssnew:',real(ss),' ',.not.reject,' ',real(alpha)


  end if !!

  !! oldpar1 is the current good value !!

  !! propose the next value, do it until a value inside bounds is generated
  inbounds = .false.
  nrej = 1
  do while (.not.inbounds)

     newpar = MCMC_propose(oldpar1,R)

     inbounds = MCMC_checkbounds(newpar) 
     if (.not.inbounds) then ! reject immediately
        drstage = 1
        isimu = isimu+1  ! addvance mcmc index
        nrej  = nrej + 1  ! addvance number of rejects
        if (verbosity>=0) write(*,*) 'outbound, rej:',nrej, 'par:',real(newpar)

     else ! check if we reject by prior

        sspri = MCMC_priorfun(newpar)
        sscrit = MCMC_sscrit(ssprev1,sspri1)
        if (sspri >= sscrit) then
           inbounds = .false.
           isimu = isimu+1  ! addvance mcmc index
           nrej  = nrej + 1  ! addvance number of rejects
           if (verbosity>=0) write(*,*) 'er, rej by prior:',nrej, 'par:',real(newpar)
           if (verbosity>=0) write(*,*) 'sscrit',real(sscrit), 'sspri',real(sspri)
        else
           sscrit = sscrit - sspri
        end if

     end if
  end do
  !! now we have a value inside the box

  chain(1,1:npar)    = oldpar1
  sschain(1,1:nycol) = ssprev1

  chain(1,ncolchain) = dble(nrej)    ! last column has the number of stays
  sschain(1,nycol+1) = dble(nrej)

  chainind = 1 ! the length of the chain
  simuind  = 1

  call MCMC_updatesigma2(ssprev1)

  !! adapt
  call MCMC_adapt(simuind)

  !! calculate new sscrit
!  sspri = MCMC_priorfun(newpar)
!  sscrit = MCMC_sscrit(ssprev2,sspri1)-sspri

  !! save the control variables to a namelist
  call write_mcmcrun_namelist()

  !! save the proposed value to mcmcparnew.dat
  call writedata('mcmcparnew.dat', reshape(newpar,(/1, size(newpar,1)/)) ) ! row vector
  call writedata('mcmcssprev1.dat',ssprev1)
  call writedata('mcmcsscrit.dat',sscrit)

  !! we want to collect also the rejected points somewhere
  if (reject) then
     call delete_file('mcmc_accepted')
     call create_empty_file('mcmc_rejected')

  else
     call delete_file('mcmc_rejected')
     call create_empty_file('mcmc_accepted')
  end if

  !! is this the last round?
  if (ieval >= nsimu) then
     call create_empty_file('mcmc_run_done')
  end if

  return

end subroutine MCMC_run1_er
