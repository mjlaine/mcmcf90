!!! $Id: MCMC_run_ram.F90,v 1.6 2012/07/03 10:57:34 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mcmc library
!!! File: MCMC_run_ram.F90
!!! Purpose: the actual mcmc simulation loop,
!!! RAM version
!!!
!!! Marko Laine 2010 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!!
!!!
subroutine MCMC_run_ram()
  implicit none

  real(kind=dbl), dimension(nycol) :: ss1, ss2
  real(kind=dbl), dimension(npar) :: oldpar, newpar
  real(kind=dbl) :: sspri1, sspri2
  real(kind=dbl) :: alpha12
  real(kind=dbl), dimension(npar,npar) :: u
  integer :: i
  logical :: reject, inbounds

  if (inited /= 1) then
     write(*,*) "we have not inited"
     stop
  end if

  if (verbosity>0) write(*,*) 'note: running RAM'

  oldpar = par0
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

     newpar = MCMC_propose_ram(oldpar,R,u)
     inbounds = MCMC_checkbounds(newpar) ! constraints
     if (.not.inbounds) then
        bndstayed = bndstayed + 1
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
     call MCMC_adapt_ram(i,u,alpha12)

  end do MCMC_LOOP

  return
end subroutine MCMC_run_ram


!!! subroutines for RAM
function MCMC_propose_ram(oldpar,R,u) result(newpar)
  implicit none
  real(kind=dbl), intent(in) :: oldpar(:)
  real(kind=dbl), intent(in) :: R(:,:)
  real(kind=dbl)  :: newpar(size(oldpar))

  real(kind=dbl), intent(out)  :: u(size(oldpar))

  u = random_normal(npar)
  if (usesvd /= 0) then 
     newpar = oldpar + matmulx(R, u)
  else
     newpar = oldpar + matmulu(R, u,'t')
  end if
end function MCMC_propose_ram


subroutine MCMC_adapt_ram(simuind,u,alpha)
  implicit none
  integer :: simuind ! where we are
  real(kind=dbl), intent(in)  :: u(npar)
  real(kind=dbl), intent(in)  :: alpha
  integer, save :: stayed2 = 0
  !! local variables
  real(kind=dbl) :: staypc
  real(kind=dbl) :: a
  integer :: info
  integer :: i, doupdate=1

  real(kind=dbl) :: ram(npar,npar)

  if (mod(simuind,printint) == 0) then
     write(*,'(A,I10,A,F5.1,A,F5.1,A)') ' simu i =', simuind, &
          ', stayed % = ',real(stayed)/real(simuind)*100.0, &
          '  [',real(bndstayed)/real(simuind)*100.0,']'
  end if

  if (doadapt == 0) then
     if (verbosity>10) write(*,*) 'note: no adaptation'
     return ! we do not adapt
  end if

  if (simuind < burnintime .and. doburnin /= 0) then
     if (verbosity>10) write(*,*) 'note: in burnin'
     return ! still doing burnin
  end if

!!! do adapt
  if (verbosity>10) write(*,*) 'note: ',simuind,': adapting'

!!! Burn in is done,  adapt

  if (verbosity>3) write(*,*) 'R:'
  if (verbosity>3) call printmat(R)

  if (doupdate == 0) then
     ram = 0.0_dbl
     do i=1,npar
        ram(i,i) = 1.0_dbl
     end do
     if (verbosity>3) call printmat(ram)
     ram = ram + 1.0_dbl/real(simuind)**(0.7_dbl)*(alpha-0.234_dbl)* &
          matmul(reshape(u,(/npar,1/)),reshape(u,(/1,npar/)))/sum(u**2)
     chaincmat = matmulxy(R,matmulxy(R,ram,'l','t'),'r','n')

     if (verbosity>3) write(*,*) 'alpha:', alpha
     if (verbosity>3) write(*,*) 'u:',u
     if (verbosity>3) write(*,*) '|u|',sum(u**2)
     if (verbosity>3) call printmat(ram)
     if (verbosity>3) call printmat(matmul(reshape(u,(/npar,1/)),reshape(u,(/1,npar/)))/sum(u**2))

     call covtor(chaincmat,R,info=info)
     if (info .ne. 0) then
        write(*,*) 'Warning: error in Chol/SVD, not adapting, info = ', info
     end if

  else
     ! need to do Cholesky update

     a = 1.0_dbl/real(simuind)**nuparam*(alpha-alphatarget)
     if (a>=0.0_dbl) then
        call cholupdate(R, u/sum(u**2) * a)
     else
        call choldowndate(R, -u/sum(u**2) * a)
     end if

  end if

  if (verbosity>3) write(*,*) 'R: '
  if (verbosity>3) call printmat(R)
  
  return
end subroutine MCMC_adapt_ram

