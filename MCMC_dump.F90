!!! $Id: MCMC_dump.F90,v 1.3 2012/07/02 09:22:55 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mcmcrun library
!!! File: MCMC_dump.F90
!!! Purpose: handle the Modest dump files
!!!
!!! Marko Laine 2008 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!!

subroutine MCMC_dump_init()
  implicit none
  integer :: fstat 
  include 'external_inc.h'
  call dump_init()
end subroutine MCMC_dump_init

subroutine MCMC_dump_end()
  implicit none
  include 'external_inc.h'
  call dump_end()
end subroutine MCMC_dump_end

subroutine MCMC_dump(oldpar)
  implicit none
  real(kind=dbl), intent(in) :: oldpar(:)
  include 'external_inc.h'
  call dump(oldpar)
end subroutine MCMC_dump

