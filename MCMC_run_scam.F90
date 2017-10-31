!!! $Id: MCMC_run_scam.F90,v 1.10 2012/11/21 11:37:39 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mcmc library
!!! File: MCMC_run_scam.F90
!!! Purpose: the actual mcmc simulation loop, scam version
!!!
!!! Marko Laine 2010 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!!
!!!
subroutine MCMC_run_scam()
  implicit none

  real(kind=dbl), dimension(nycol) :: ss1, ss2
  real(kind=dbl) :: sspri1, sspri2
  real(kind=dbl), dimension(npar) :: oldpar, newpar
  real(kind=dbl) :: alpha12
  integer :: i, j
  logical :: reject, inbounds, rejall

  if (inited /= 1) call doerror('we have not inited')

  if (verbosity>0) write(*,*) 'note: doing SCAM'

  oldpar = par0
  sspri1 = MCMC_priorfun(oldpar)
  ss1    = MCMC_ssfunction(oldpar)
  write(*,*) 'ss1 = ',ss1, ' sspri1 = ', sspri1

  rejall = .false.
  call MCMC_savechain(oldpar,ss1,sspri1,sigma2,rejall)

  call MCMC_userfun(oldpar)
  call MCMC_dump_init()

  !! the simulation loop
  MCMC_LOOP: do i = 2, nsimu

     simuind = i              ! global simuind

     rejall = .true. ! all rejected
     !! one component at a time loop
     PAR_LOOP: do j = 1,npar

        newpar   = MCMC_propose_sc(oldpar,j) !!!!!
        inbounds = MCMC_checkbounds(newpar) ! constraints
        if (.not.inbounds) then
           if (.not.dodr) bndstayed = bndstayed + 1
           ss2     = huge(0.0d0) ! realmax
           alpha12 = 0.0d0
           reject  =.true.
        else
           sspri2  = MCMC_priorfun(newpar)
           ss2     = MCMC_ssfunction(newpar)
           alpha12 = MCMC_alpha(oldpar, ss1, sspri1, newpar, ss2, sspri2)
           reject  = MCMC_reject(alpha12)
        end if
        if (verbosity > 2) then
           write(*,*) i,'/',j,' oldpar:', real(oldpar), ' newpar:',real(newpar), &
                'ss:',real(ss2),' ',.not.reject
        end if
        if (.not.reject) then
           ss1    = ss2
           sspri1 = sspri2
           oldpar = newpar
           rejall = .false.
        end if

        !! DR code removed

     end do PAR_LOOP

     if (rejall) then ! we rejected all the components
        stayed = stayed + 1
     else
        call MCMC_userfun(oldpar)
     end if

     sschain(chainind,nycol+1) = chain(chainind,ncolchain) ! copy the count

     call MCMC_dump(oldpar)
     call MCMC_updatesigma2(ss1)
     call MCMC_savechain(oldpar,ss1,sspri1,sigma2,rejall)
     !! must use svd in adapt
     call MCMC_adapt(i) ! adapt R if necessary, called on every simuloop
     
  end do MCMC_LOOP

  return
end subroutine MCMC_run_scam


function  MCMC_propose_sc(oldpar,j) result(newpar)
  implicit none
  real(kind=dbl),intent(in) :: oldpar(:)
  integer, intent(in) :: j
  real(kind=dbl) :: newpar(size(oldpar))

!   real(kind=dbl),intent(in) :: R(:,:) ! uses global R

  !! local variables
  real(kind=dbl) :: rotpar(size(oldpar)), z(1)

  ! rotation
  rotpar = MCMC_scam_rotate(oldpar,R,'f')
  z = random_normal(1)*qcovstd(j)
  rotpar(j) = rotpar(j) + z(1)

! testing without 1. rotation
!  rotpar = oldpar
!  z = random_normal(1)*qcovstd(j)
!  rotpar(j) = rotpar(j) + z(1)

  newpar = MCMC_scam_rotate(rotpar,R,'b')

end function MCMC_propose_sc

!!! rotate the parameter vector
!!! 'f' rotates oldpar from N(0,R'R) to N(0,I*s) variable
!!! 'b' rotates backwards
function MCMC_scam_rotate(oldpar,R,dire) result(newpar)
  implicit none
  real(kind=dbl),intent(in) :: oldpar(:)
  real(kind=dbl),intent(in) :: R(:,:)
  character(len=1),intent(in) :: dire
  real(kind=dbl) :: newpar(size(oldpar))

  if (usesvd /= 1) stop 'must use svd with scam'
  if ( dire.eq.'f'.or.dire.eq.'F' ) then
     newpar = matmulx(R, oldpar,'t')
  elseif ( dire.eq.'b'.or.dire.eq.'B' ) then
     newpar = matmulx(R, oldpar,'n')
  else
     stop 'wrong direction in rotate'
  end if

end function MCMC_scam_rotate

!!!
!!!
function MCMC_scam_update(oldpar,R,ipar,d) result(newpar)
  implicit none
  real(kind=dbl),intent(in) :: oldpar(:)
  real(kind=dbl),intent(in) :: R(:,:)
  integer, intent(in) :: ipar
  real(kind=dbl),intent(in) :: d
  real(kind=dbl) :: newpar(size(oldpar))

  if (usesvd /= 1) stop 'must use svd with scam'
  newpar = oldpar + R(:,ipar)*d

end function MCMC_scam_update
!!!
!!! SCAM version of ssfun.
!!! Calculate -2*log(likelihood(newpar)) aka 'ssfun' using
!!! user written ssfunction
!!!
function MCMC_ssfunction_scam(newpar, itask) result(ss)
  implicit none
  real(kind=dbl), intent(in) :: newpar(:)
  integer, intent(in) :: itask
  real(kind=dbl), dimension(nycol) :: ss

  ss = 0.0d0
  stop 'do not use MCMC_ssfunction_scam yet'

!!  ss = ssfunction(newpar,npar,nycol)
!  ss = ssfunction_scam(newpar,npar,nycol,ipar,p,d)
end function MCMC_ssfunction_scam
