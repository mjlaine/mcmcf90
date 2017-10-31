!!! $Id: MCMC_run.F90,v 1.14 2012/07/03 10:57:33 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mcmc library
!!! File: MCMC_run.F90
!!! Purpose: the actual mcmc simulation loop
!!!
!!!
!!! Marko Laine 2001 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!!
subroutine MCMC_run()
  implicit none

  real(kind=dbl), dimension(nycol) :: ss1, ss2, ss3
  real(kind=dbl), dimension(npar) :: oldpar, newpar, newpar2
  real(kind=dbl) :: sspri1, sspri2, sspri3
  real(kind=dbl) :: alpha12, alpha13
  integer :: i
  logical :: reject, inbounds

  if (inited /= 1) call doerror('we have not inited')

  if (verbosity>0.and.dodr) write(*,*) 'note: running DRAM with drscale',real(drscale)
  if (verbosity>0.and..not.dodr) write(*,*) 'note: running AM'

  oldpar = par0 ! par0 is global in mcmcmod and inited in MCMC_init
  if (verbosity>2) write(*,*) 'note: calling ssfunction for the first time'
  sspri1 = MCMC_priorfun(oldpar)
  ss1    = MCMC_ssfunction(oldpar)
  if (verbosity>0) write(*,*) 'note: ss1 = ',real(ss1), &
       ' sspri1 = ',real(sspri1)

  reject = .false.
  call MCMC_savechain(oldpar,ss1,sspri1,sigma2,reject)

  call MCMC_userfun(oldpar)
  call MCMC_dump_init()

  !! the simulation loop
  MCMC_LOOP: do i = 2, nsimu

     simuind = i              ! global simuind
     if (verbosity>2) write(*,*) 'isimu:', i

     newpar = MCMC_propose(oldpar,R)
     inbounds = MCMC_checkbounds(newpar) ! constraints
     if (.not.inbounds) then
        if (.not.dodr) bndstayed = bndstayed + 1
        ss2     = huge(0.0d0) ! realmax
        sspri2  = huge(0.0d0)
        alpha12 = 0.0d0
        reject  =.true.
     else
        sspri2  = MCMC_priorfun(newpar)
        ss2     = MCMC_ssfunction(newpar)
        alpha12 = MCMC_alpha(oldpar, ss1, sspri1, newpar, ss2, sspri2)
        reject  = MCMC_reject(alpha12)
     end if

     if (verbosity>3) write(*,*) 'oldpar:', real(oldpar)
     if (verbosity>3) write(*,*) 'newpar:', real(newpar)
     if (verbosity>2) write(*,*) 'ss1:', real(ss1),' ss2:', real(ss2), 'alpha12:', real(alpha12),.not.reject

     if (reject .and. dodr) then ! We reject but make one new DR try
        if (verbosity>2) write(*,*) 'dr try'
        drtries  = drtries+1
        newpar2  = MCMC_propose(oldpar,R2)
        inbounds = MCMC_checkbounds(newpar2)
        if (.not.inbounds) then
           bndstayed = bndstayed + 1
           reject    = .true.
        else
           sspri3  = MCMC_priorfun(newpar2)
           ss3     = MCMC_ssfunction(newpar2)
           alpha13 = MCMC_DR_alpha13(oldpar, ss1, sspri1, &
                newpar,ss2, sspri2, alpha12, &
                newpar2,ss3, sspri3)
           reject  = MCMC_reject(alpha13)

           if (.not.reject) then ! DR accepted
              draccepted = draccepted + 1
              newpar     = newpar2
              ss2        = ss3
              sspri2     = sspri3
           end if
        end if
        if (verbosity>3) write(*,*) 'newpar2:', real(newpar2)
        if (verbosity>2) write(*,*) 'ss1:', real(ss1),'ss2:', real(ss2),'ss3:', &
             real(ss3), 'alpha13:', real(alpha13),.not.reject
     end if ! DR

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

  if (dodr) then
     if (verbosity>0) write(*,*) 'note: draccepted=',draccepted,'drtries = ',drtries
  end if

  return
end subroutine MCMC_run
