!!! $Id: MCMC_DRAM.F90,v 1.8 2012/07/02 09:22:54 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mcmc library
!!! File: MCMC_DRAM.F90
!!! Purpose: Modular subroutines for the MCMC clculations
!!!
!!! Marko Laine 2008 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!!
!!! contains:
!!! MCMC_propose, MCMC_checkbounds, MCMC_ssfunction, MCMC_alpha, MCMC_reject, 
!!! MCMC_DR_alpha13, MCMC_updatesigma2

!!!
!!! propose a new value given current value and proposal Cholesky factor
!!! MH random walk proposal using multivariate Gaussian with covariance
!!! R'*R, and the current center at oldpar
!!!
function MCMC_propose(oldpar,R) result(newpar)
  implicit none
  real(kind=dbl), intent(in) :: oldpar(:)
  real(kind=dbl), intent(in) :: R(:,:)
  real(kind=dbl)  :: newpar(size(oldpar))

  if (usesvd /= 0) then 
     newpar = oldpar + matmulx(R, random_normal(npar))
  else
     newpar = oldpar + matmulu(R, random_normal(npar),'t')
  end if
end function MCMC_propose

!!!
!!! check if parameter newpar is inside bounds defined by
!!! user written subroutine checkbounds
!!!
function MCMC_checkbounds(newpar) result(inbounds)
  implicit none
  real(kind=dbl), intent(in) :: newpar(:)
  logical :: inbounds

  include 'external_inc.h'

  inbounds = checkbounds(newpar)

end function MCMC_checkbounds

!!!
!!! Calculate -2*log(likelihood(newpar)) aka 'ssfun' using user written
!!! ssfunction
!!!
function MCMC_ssfunction(newpar) result(ss)
  implicit none
  real(kind=dbl), intent(in) :: newpar(:)
  real(kind=dbl), dimension(nycol) :: ss

  include 'external_inc.h'

  ss = ssfunction(newpar,npar,nycol)
end function MCMC_ssfunction


!!!
!!! Calculate -2*log(likelihood(newpar)) aka 'ssfun' using user written
!!! ssfunction, early rejection version. 
!!! 
function MCMC_ssfunction_er(newpar,sscrit) result(ss)
  implicit none
  real(kind=dbl), intent(in) :: newpar(:), sscrit
  real(kind=dbl), dimension(nycol) :: ss

  include 'external_inc.h'

  ss = ssfunction_er(newpar,npar,nycol,sscrit)
end function MCMC_ssfunction_er


!!!
!!! Calculate -2*log(prior(newpar)) 
!!!
function MCMC_priorfun(newpar) result(sspri)
  implicit none
  real(kind=dbl), intent(in) :: newpar(:)
  real(kind=dbl)  :: sspri

  include 'external_inc.h'

  sspri = priorfun(newpar,npar)

end function MCMC_priorfun

!!!
!!! calculate the acceptance probability
!!! oldpar : previous velue, newpar: proposed value, 
!!! ss1: -2*log(lik(oldpar)),  ss2: -2*log(lik(newpar))
!!! 
!!! prior is calculated by priorfun, the default uses independent Gaussian 
!!! priors, whose parameters are read from priorsfile
!!! priors are in "ss" form, so -2*log(p(par))
function MCMC_alpha(oldpar, ss1, sspri1, newpar, ss2, sspri2) result(alpha12)
  implicit none
  real(kind=dbl), intent(in), dimension(:) :: oldpar, ss1, newpar, ss2
  real(kind=dbl), intent(in) :: sspri1, sspri2
  real(kind=dbl)  :: alpha12

  real(kind=dbl) :: tst

  tst = -0.5_dbl*(sum((ss2-ss1)/sigma2) + (sspri2 - sspri1))

  if (tst >= 0.0_dbl) then
     alpha12 = 1.0e0_dbl
  elseif (tst < log_realmin) then 
     alpha12 = 0.0e0_dbl
  else
     alpha12 = exp(tst)
  end if

end function MCMC_alpha

!!!
!!! critical value of -2*log(lik), if the actual value is larger, then we
!!! will reject
!!!
function MCMC_sscrit(ss1, priss1) result(sscrit)
  implicit none
  real(kind=dbl), intent(in), dimension(:) :: ss1
  real(kind=dbl), intent(in) :: priss1
  real(kind=dbl)  :: sscrit

  real(kind=dbl) :: u

  call random_number(u)
  sscrit = -2.0_dbl*log(u) + sum(ss1/sigma2) + priss1

end function MCMC_sscrit

!!!
!!! calculate wether we recject the point given the acceptance probability
!!!
function MCMC_reject(alpha) result(reject)
  implicit none
  real(kind=dbl), intent(in) :: alpha
  logical :: reject

  real(kind=dbl) :: u

  reject=.true.
  if (alpha >= 1.0e0_dbl) then
     reject=.false.
  elseif (alpha > 0.0e0_dbl) then
     call random_number(u)
     if (u<=alpha) reject = .false.
  end if

end function MCMC_reject

!!!
!!! second stage DR acceptance probability
!!! assumes Gaussian proposals, global variable iC = inv(C) is the inverse of
!!! first stage proposal of the point that was rejected
!!! 
function MCMC_DR_alpha13(oldpar,ss1, sspri1, newpar,ss2, sspri2, alpha12, &
     newpar2, ss3, sspri3) result(alpha13)

  implicit none
  real(kind=dbl), intent(in) :: ss1(:), ss2(:), ss3(:)
  real(kind=dbl), intent(in) :: sspri1, sspri2, sspri3
  real(kind=dbl), intent(in) :: oldpar(:), newpar(:), newpar2(:)
  real(kind=dbl), intent(in) :: alpha12
  real(kind=dbl) :: alpha13

  real(kind=dbl) :: tst32, l2, q1, alpha32

  if (alpha12 == 0.0e0_dbl) then
     alpha32 = 0.0e0_dbl
  else
     tst32 = -0.5_dbl*(sum((ss2-ss3)/sigma2) + (sspri2-sspri3))
     alpha32 = min(1.0e0_dbl,exp(tst32))
  end if
  l2 = -0.5_dbl*(sum((ss3-ss1)/sigma2) + (sspri3-sspri1) )
  q1 = -0.5_dbl*( &
       sum(matmuls(iC,newpar2-newpar)*(newpar2-newpar)) - &
       sum(matmuls(iC,oldpar-newpar)*(oldpar-newpar)))
  alpha13 = min(1.0e0_dbl, exp(l2+q1)*(1.0_dbl-alpha32)/(1.0_dbl-alpha12))

end function MCMC_DR_alpha13

!!!
!!! This routine updates the error variance sigma2 for Gaussian error models
!!! with conjugate invChisq(N0,S02) prior for the variance
!!!
subroutine MCMC_updatesigma2(ss)
  implicit none
  real(kind=dbl), intent(in) :: ss(:)

  real(kind=dbl) :: g(1)
  integer :: j

  if (updatesigma /= 0) then
     do j=1, nycol
        g = random_gamma(1,N0/2.0_dbl+dble(nobs(j))/2.0_dbl, 2.0_dbl/(N0*S02+ss(j)))
        sigma2(j) = 1.0_dbl / g(1)
     end do
  end if

end subroutine MCMC_updatesigma2

