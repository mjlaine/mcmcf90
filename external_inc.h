!!! external routines in mcmcrun library                           -*- mode: F90; -*-

!!! User specified prior, see the default prior in priorfun.f90.
  interface
     function priorfun(theta,len)
       real*8 priorfun
       integer*4 len
       real*8 theta(len)
     end function priorfun
  end interface

!!! User defined -2*log(likelihood) functions and bounds check function.
!!! The default versions are in ssfunction0.f90 and checkbounds0.f90.
  interface
     function ssfunction(theta,npar,ny)
       integer*4 npar, ny
       real*8 theta(npar)
       real*8 ssfunction(ny)
     end function ssfunction
     function ssfunction_er(theta,npar,ny,sscrit)
       integer*4 npar, ny
       real*8 theta(npar), sscrit
       real*8 ssfunction_er(ny)
     end function ssfunction_er
     function checkbounds(theta)
       real*8 theta(:)
       logical checkbounds
     end function checkbounds
  end interface

!!! The initialize function is used to set the initial 
!!! values for some variables. The default subroutine uses files.
  interface
     subroutine initialize(par0,npar,cmat0,initcmatn,sigma2,nobs,nycol)
       use mcmcprec
       implicit none
       integer, intent(inout) :: npar, initcmatn, nycol
       real(kind=dbl), intent(inout), allocatable :: par0(:), cmat0(:,:)
       real(kind=dbl), intent(inout), allocatable :: sigma2(:)
       integer, intent(inout), allocatable :: nobs(:)
     end subroutine initialize
  end interface

!!! Alternative to initialize is setparams
  interface
     subroutine setparams(par,qpct,sig2,n)
       use mcmcprec
       implicit none
       real(kind=dbl),intent(in) :: par(:)
       real(kind=dbl),intent(in), optional :: qpct
       real(kind=dbl),intent(in), optional :: sig2
       integer, intent(in), optional :: n
     end subroutine setparams
  end interface

!!! these dump routines are for advanced use only
  interface 
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
  end interface

