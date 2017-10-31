!!! $Id: MCMC_run_er.F90,v 1.9 2012/07/03 10:57:34 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mcmc library
!!! File: MCMC_run_er.F90
!!! Purpose: the actual mcmc simulation loop, early rejection version
!!!
!!! Marko Laine 2010 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!!
!!!
subroutine MCMC_run_er()
  implicit none

  real(kind=dbl), dimension(nycol) :: ss1, ss2
  real(kind=dbl), dimension(npar) :: oldpar, newpar
  real(kind=dbl) :: sspri1, sspri2, sscrit
  integer :: i
  logical :: reject, inbounds, prirejected

  if (inited /= 1) call doerror('we have not inited')

  if (verbosity>0) write(*,*) 'note: Early rejection for prior'

  if (dodr.eqv..true.) then
     if (verbosity>0) write(*,*) 'note: no dr with er'
     dodr = .false.
  end if

  if (nycol > 1) then
     write(*,*) 'note: nycol should not be greater than 1 with er!'
  end if

  oldpar = par0
  if (verbosity>2) write(*,*) 'note: calling ssfunction for the first time'
  sspri1 = MCMC_priorfun(oldpar)
  ss1    = MCMC_ssfunction(oldpar)
  if (verbosity>0) write(*,*) 'note: ss1 = ',real(ss1),' sspri1 = ',real(sspri1)

  reject = .false.
  call MCMC_savechain(oldpar,ss1,sspri1,sigma2,reject)
  call MCMC_userfun(oldpar)
  call MCMC_dump_init()

  !! the simulation loop
  MCMC_LOOP: do i = 2, nsimu

     simuind = i              ! global simuind
     if (verbosity>2) write(*,*) 'isimu:', i

     prirejected = .false.

     newpar = MCMC_propose(oldpar,R)     
     inbounds = MCMC_checkbounds(newpar) ! constraints
     if (.not.inbounds) then
        bndstayed = bndstayed + 1
        reject  =.true.
        if (verbosity>2) write(*,*) 'outbound, newpar:', real(newpar)
     else
        sscrit = MCMC_sscrit(ss1,sspri1)
        sspri2 = MCMC_priorfun(newpar)
        if (verbosity>2) write(*,*) 'sspri2:', real(sspri2), ' sscrit:', real(sscrit)
        if (sspri2 >= sscrit) then
           reject  =.true.
           prirejected = .true.
           erstayed = erstayed + 1
           if (verbosity>2) write(*,*) 'note: rejected due sspri2 >= sscrit ', &
                real(sspri2),real(sscrit)
        else

           !! problem here if nycol > 1
           sscrit = sigma2(1)*(sscrit - sspri2)
           if (verbosity>3) write(*,*) 'note: sscrit is now:', real(sscrit)
           ss2 = MCMC_ssfunction_er(newpar,sscrit)
           if (sum(ss2) >= sscrit) then
              if (verbosity>2) write(*,*) 'note: rejected due ss2 >= sscrit ', &
                real(ss2),real(sscrit)
              reject = .true.
           else
              reject = .false.
           end if

           if (verbosity>3) write(*,*) 'oldpar:', real(oldpar)
           if (verbosity>3) write(*,*) 'newpar:', real(newpar)
           if (verbosity>2) write(*,*) 'ss1:', real(ss1),' ss2:', real(ss2), 'sscrit', real(sscrit),.not.reject

        end if
     end if

     if (reject) then
        stayed = stayed + 1
     else
        ss1    = ss2
        sspri1 = sspri2
        oldpar = newpar
        call MCMC_userfun(newpar)
     end if

     call MCMC_dump(oldpar)
     call MCMC_updatesigma2(ss1)
     call MCMC_savechain(oldpar,ss1,sspri1,sigma2,reject)
     call MCMC_adapt(i) ! adapt R if necessary, called on every simuloop

  end do MCMC_LOOP

  return
end subroutine MCMC_run_er
