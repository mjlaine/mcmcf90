!!! $Id: MCMC_run1.F90,v 1.17 2012/07/03 10:57:34 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mcmc library
!!! File: MCMC_run1.F90
!!! Purpose: MCMC loop done once
!!!
!!!
!!! Marko Laine 2010 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!!
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
subroutine MCMC_run1()

  use mcmcrun1

  implicit none


  real(kind=dbl) :: ss(nycol)
  real(kind=dbl) :: newpar(npar)
  real(kind=dbl) :: sspri
  real(kind=dbl) :: alpha
  logical :: reject, inbounds

  real(kind=dbl) :: oldpar1(npar), oldpar2(npar)
  real(kind=dbl) :: ssprev1(nycol), ssprev2(nycol)
  real(kind=dbl) :: sspri1, sspri2

  integer :: fstat

!! these are in module mcmcrun1
!  real(kind=dbl) :: alpha12
!  integer :: drstage=1, isimu = 1, nrej=0, ieval = 0
!  namelist /mcmcrun/ drstage, isimu, ieval, nrej, alpha12

  if (inited /= 1) call doerror('we have not inited')

  if (updatesigma==1) then
     write(*,*) 'Warning: do not use updatesigma in ER!'
     updatesigma=0
  end if

  !! open namelist to get the values from previous run
  call read_mcmcrun_namelist()

  if (isimu == 1) then
     if (verbosity>=0) write(*,*) 'first simulation'

     newpar  = par0
     oldpar1 = par0
     oldpar2 = par0
     sspri   = MCMC_priorfun(newpar)
     ss      = MCMC_ssfunction(newpar)
     ieval = ieval + 1
     if (verbosity>=0) write(*,*) 'ss1 = ',real(ss), ' sspri1 = ', real(sspri)

     ssprev1 = ss
     ssprev2 = ss

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
     if (drstage > 1 .and. dodr) then ! we return for DRstage==2
        call readdata('mcmcoldpar2.dat',oldpar2)
        call readdata('mcmcoldpar1.dat',oldpar1)
        call readdata('mcmcssprev2.dat',ssprev2)
        call readdata('mcmcssprev1.dat',ssprev1)
     else
        call readdata('mcmcparf.dat',oldpar1)
        call readdata('mcmcssprev1.dat',ssprev1)
     end if

     !! par0 has the intial value which should now have the previous chainmean
     !! need to fix MCMC_adapt???
     call readdata(meanfile,par0,uselock=.true.)

     simuind = isimu              ! global simuind ???

     newpar = par0 ! from 'mcmcpar.dat'

     !! evaluate the latest ss and prior
     sspri = MCMC_priorfun(newpar)
     ss  = MCMC_ssfunction(newpar)
     ieval = ieval + 1

     if (verbosity>=0) write(*,*) 'isimu:',isimu,' drstage:',drstage

     if (drstage>1 .and. dodr) then
        if (verbosity>=0) write(*,*) 'oldpar2:',real(oldpar2),' ssprev2:',real(ssprev2),' alpha12:',real(alpha12)
     end if
     if (verbosity>=0) write(*,*) 'oldpar1:',real(oldpar1),' ssprev1:',real(ssprev1)

     !! we will need to recalculate the old priors for now
     sspri1 = MCMC_priorfun(ssprev1)
     if (verbosity>=3) write(*,*) ' sspri1:',real(sspri1)


     !! do we accept it?
     if (drstage>1 .and. dodr) then
        sspri2 = MCMC_priorfun(ssprev2)
        if (verbosity>=3) write(*,*) 'sspri2:',real(sspri2)
        alpha = MCMC_DR_alpha13(oldpar2, ssprev2, sspri2, oldpar1, &
             ssprev1, sspri1, alpha12, &
             newpar, ss, sspri )
     else
        alpha = MCMC_alpha(oldpar1, ssprev1, sspri1, newpar, ss, sspri)
     end if
     reject  = MCMC_reject(alpha)

     if (reject) then
        if (dodr) then
           if (drstage == 1) then ! goto the next DR stage
              drstage = 2
              ssprev2 = ssprev1
              ssprev1 = ss
              oldpar2 = oldpar1
              oldpar1 = newpar
              alpha12 = alpha
           else ! (drstage == 2) then ! rejected the second stage also
              isimu = isimu + 1
              drstage = 1
              ssprev1 = ssprev2
              oldpar1 = oldpar2
           end if
        else
           isimu = isimu + 1
           drstage = 1
        end if
     else ! accept
        drstage = 1
        isimu = isimu + 1
        nrej = 1
        oldpar1 = newpar
        ssprev1 = ss
        alpha12 = alpha

        ssprev2 = ssprev1
        oldpar2 = oldpar1

     end if

     if (verbosity>=0) write(*,*) 'newpar:',real(newpar),'ssnew:',real(ss),' ',.not.reject,' ',real(alpha)


  end if

  !! oldpar2 is the current good value !!

  !! propose the next value, do it until a value inside bounds is generated
  inbounds = .false.
  nrej = 1
  do while (.not.inbounds)

     if (drstage > 1 .and. dodr) then
        newpar = MCMC_propose(oldpar2,R2)
     else
        newpar = MCMC_propose(oldpar2,R)
     end if

     inbounds = MCMC_checkbounds(newpar) 
     if (.not.inbounds) then ! reject immediately
        if (drstage == 1 .and. dodr) then
           drstage = 2
           oldpar1 = newpar
           ssprev1 = HUGE(ssprev1)
           alpha12 = 0.0d0
           if (verbosity>=0) write(*,*) 'outbound, stage 1, par:',real(newpar)
        else
           drstage = 1
           isimu = isimu+1  ! addvance mcmc index
           nrej  = nrej + 1  ! addvance number of rejects
           oldpar1 = oldpar2
           ssprev1 = ssprev2
           if (verbosity>=0) write(*,*) 'outbound, rej:',nrej, 'par:',real(newpar)
        end if
     end if
  end do
  !! now we have a value inside the box

  chain(1,1:npar)    = oldpar2
  sschain(1,1:nycol) = ssprev2

  chain(1,ncolchain) = dble(nrej)    ! last column has the number of stays
  sschain(1,nycol+1) = dble(nrej)

  chainind = 1 ! the length of the chain
  simuind  = 1

  call MCMC_updatesigma2(ssprev1)

  !! adapt, but not if
!  if (dodr == 0 .or. .not.reject .or. drstage == 1)) then
     call MCMC_adapt(simuind)
!  end if

  !! save the control variables to a namelist
  call write_mcmcrun_namelist()

  !! save the proposed value to mcmcparnew.dat
  call writedata('mcmcparnew.dat', reshape(newpar,(/1, size(newpar,1)/)) ) ! row vector
  call writedata('mcmcssprev1.dat',ssprev1)
  if (dodr) call writedata('mcmcoldpar2.dat',oldpar2)
  if (dodr) call writedata('mcmcoldpar1.dat',oldpar1)
  if (dodr) call writedata('mcmcssprev2.dat',ssprev2)
  !  if (dodr) call writedata('mcmcssprev1.dat',ssprev1)

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

end subroutine MCMC_run1
