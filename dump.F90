!!! $Id: dump.F90,v 1.1 2012/07/02 06:25:00 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mcmc library
!!! File: dump.F90
!!! Purpose: dummy version of the dump routines
!!!
!!!
!!! Marko Laine 2012 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!!

subroutine dump_init()
  implicit none
end subroutine dump_init

subroutine dump_end()
  implicit none
end subroutine dump_end

subroutine dump(oldpar)
  use mcmcprec
  implicit none
  real(kind=dbl), intent(in) :: oldpar(:)
end subroutine dump

