!!! $Id: MCMC_userfun.F90,v 1.2 2012/06/27 10:10:37 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mcmc library
!!! File: MCMC_userfun.F90
!!! Purpose: call modest userfun inside MCMC loop
!!!
!!! Marko Laine 2008 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!!

subroutine MCMC_userfun(oldpar)
  implicit none

  real(kind=dbl), intent(in) :: oldpar(:)

  !! does not do anything at the moment

end subroutine MCMC_userfun
