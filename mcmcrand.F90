!!! $Id: mcmcrand.F90,v 1.22 2012/06/27 10:10:37 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! this file of part of mcmc library
!!! File: mcmcrand.f90
!!! Purpose: random numbers module for the mcmc module
!!!          provides normal and gamma variates
!!!          uses the default random_number subroutine for univariate deviates
!!!          defines random_initialize and random_eoj for controlling the seed
!!!
!!! Marko Laine 2001 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!! 
module mcmcrand

  use mcmcprec
  implicit none

  integer, save :: random_verbosity = 0
  integer, save :: save_initial_seed = 0
  integer, save :: save_final_seed = 1
  integer, save :: load_initial_seed = 1

  logical, save :: inited = .false. ! have we inited the module

  character(len=*), parameter :: random_version="$Revision: 1.22 $ Compiled: "//__DATE__

#if defined(__IFC) || defined(__INTEL_COMPILER)
  character(len=*), parameter :: seedfile="ifc_seed.dat"
  character(len=*), parameter :: seedfile0="ifc_seed0.dat"
#elif defined(__PGI)
  character(len=*), parameter :: seedfile="pgi_seed.dat"
  character(len=*), parameter :: seedfile0="pgi_seed0.dat"
#elif defined(__GFORTRAN__)
  character(len=*), parameter :: seedfile="gfortran_seed.dat"
  character(len=*), parameter :: seedfile0="gfortran_seed0.dat"
#else
  character(len=*), parameter :: seedfile="random_seed.dat"
  character(len=*), parameter :: seedfile0="random_seed0.dat"
#endif

  integer(kind=ik4), private, save :: int32seed = 1234

  private gammar_mt, normal_bm, inited !, normal_zig

contains

!!! uniform
  function random_uniform(n) result(rn)
    implicit none
    integer, intent(in) :: n
    real(kind=dbl), dimension(n) :: rn

!    if (.not.inited) call random_initialize()
    call random_number(rn)

  end function random_uniform

!!! normal
  function random_normal(n,mu,sigma2) result(rn)
    implicit none
    integer, intent(in) :: n
    real(kind=dbl), dimension(n) :: rn
    real(kind=dbl), intent(in), optional :: mu, sigma2

    real(kind=dbl) xmu, xsigma
    integer i

    !! should we call randon_initialize?
    !! if yes, then we also should call random_eoj at the end
!    if (.not.inited) call random_initialize()

    xmu = 0.0e0_dbl
    xsigma = 1.0e0_dbl
    if (present(mu)) xmu=mu
    if (present(sigma2)) xsigma=sqrt(sigma2)
    do i=1,n
       ! use simple box-muller
       rn(i) = normal_bm()*xsigma+xmu
       ! ziggurat method (needs more testings)
 !      rn(i) = normal_zig()*xsigma+xmu
    end do
  end function random_normal

!!! gamma, check parameters
  function random_gamma(n,a,b) result(gr)
    implicit none
    integer, intent(in) :: n
    real(kind=dbl), dimension(n) :: gr
    real(kind=dbl), intent(in), optional :: a, b
    real(kind=dbl) :: xa, xb
    integer :: i
    real(kind=dbl) :: u

!    if (.not.inited) call random_initialize()

    xa=1.0e0_dbl
    xb=1.0e0_dbl
    if (present(a)) xa=a
    if (present(b)) xb=b
 
    do i=1,n
       if (xa<1.0_dbl) then
          call random_number(u)
          gr(i) = gammar_mt(1.0_dbl+xa,xb)*u**(1.0_dbl/xa)
       else
          gr(i)=gammar_mt(xa,xb)
       end if
    end do

  end function random_gamma
!!!
!!!  Gamma variates using method of Marsaglia and Tsang (2000)
!!!
!!! G. Marsaglia and W. W. Tsang:
!!! A Simple Method for Generating Gamma Variables,
!!! ACM Transactions on Mathematical Software, Vol. 26, No. 3,
!!! September 2000, 363-372.
!!!
  function gammar_mt(a,b) result(y)
    implicit none
    real(kind=dbl), intent(in) :: a
    real(kind=dbl), intent(in), optional :: b
    real(kind=dbl) y

    !! local variables
    real(kind=dbl) aa, bb, d, c, x(1), u, v

    aa = a
    if (present(b)) then
       bb = b
    else
       bb = 1.0_dbl
    end if
!! there is some problem when a < 1 but as this will not happen
!! in mcmc runs I will check it later ???
    if (aa<1.0_dbl) then
       call random_number(u)
       bb = bb*u**(1.0_dbl/aa)
       aa = aa+1.0_dbl
! should we use recursive call instead?
!       y = gammar_mt(1+aa,bb)*u**(1.0_dbl/aa)
!       return
       !! this should not happen now
       write(*,*) 'Warning: gamma rnd with a < 1 is not correct'
    end if
    d = aa-1.0_dbl/3.0_dbl
    c = 1.0_dbl/sqrt(9.0_dbl*d)
    do
       do
          x = random_normal(1)
          v = 1.0_dbl+c*x(1)
          if (v > 0.0_dbl) exit
       end do
       v = v**3
       call random_number(u)
       if (u < 1.0_dbl-0.0331_dbl*x(1)**4) exit
       if (log(u) < 0.5_dbl*x(1)**2+d*(1.0_dbl-v+log(v))) exit
    end do
    y = bb*d*v

  end function gammar_mt
!!
!! simple box-muller normal random number generator from csc f95 book
!!
  function normal_bm() result(y)
    implicit none
    real (kind=dbl) :: y
    real (kind=dbl), dimension(2) :: x
    real (kind=dbl) :: xx, z
    
    real (kind=dbl), save :: saved_y
    logical, save :: saved = .false.

    if (.not. saved) then
       do
          call random_number(x)
          x(:) = 2.0_dbl*x(:) - 1.0_dbl
          xx = x(1)**2 + x(2)**2
          if (xx < 1.0_dbl .and. xx /= 0.0_dbl)  exit
       end do
       z = sqrt(-2.0_dbl*log(xx)/xx)
       saved_y = z*x(1)
       saved = .true.
       y = z*x(2)
    else
       saved = .false.
       y = saved_y
    end if
  end function normal_bm
!!!
!!! initialize random seeds
!!!
  subroutine random_initialize()
    implicit none

    real(kind=dbl) u
    integer, allocatable :: seed(:)
    integer :: fstat, i
    integer, dimension(8) :: datetime
    logical :: exist, seedok

    if (inited) return ! already inited

    if (random_verbosity>1) &
         write(*,*) 'Random module: '//trim(random_version)

    !! do the default initialization first
    !!    (again some problems with the PGI compiler on Cray)
#ifndef __PGI
    call random_seed()
#endif
    !! try to read thee seed file
    seedok = .false.
    call random_seed(size=i)
    allocate(seed(i))

    if (load_initial_seed /= 0) then
       open(123,file=seedfile,status="old",iostat=fstat)
       if ( fstat .eq. 0 ) then
          read(123,*,iostat=fstat) seed(1:i)
          if ( fstat .eq. 0 ) then
             call random_seed(put=seed(1:i))
             seedok = .true.
             if (random_verbosity>0) &
                  write(*,*) 'note: read the seed from '//trim(seedfile)
          else
             if (random_verbosity>=0) &
                  write(*,*) 'warning: error reading seed file '//trim(seedfile)
          end if
          close(123)
       else
          if (random_verbosity>0) then
             write(*,*) 'note: error opening seed file '//trim(seedfile)
             write(*,*) 'note: will create it at the end of job'
          end if
       end if
    end if

    !! if no seedfile, then try to set the initial
    if (.not.seedok) then
       !! read the seed from /dev/urandom
       inquire(file="/dev/urandom",exist=exist)
       if (exist) then
#if defined(__GFORTRAN__)
          open(123,file="/dev/urandom",access="stream", & 
               form="unformatted",iostat=fstat)
#elif defined(__INTEL_COMPILER)
          open(123,file="/dev/urandom",access="sequential",&
               form="binary",status="old",action="read",& 
               recordtype="stream",iostat=fstat)
#elif defined(AIX) || defined(RS6K)
          open(123,file="/dev/urandom",iostat=fstat)
#else
          open(123,file="/dev/urandom",form="binary",iostat=fstat)
#endif
          if ( fstat .eq. 0 ) then
             read(123,iostat=fstat) seed
             if ( fstat .eq. 0 ) then
                call random_seed(put=seed)
                if (random_verbosity>0) &
                     write(*,*) 'note: used /dev/urandom to initialize the seed'
                if (random_verbosity>1) &
                     write(*,*) seed
                seedok = .true.
             else
                if (random_verbosity>0) then
                   write(*,*) 'note: failed to read /dev/urandom'
                   write(*,*) 'iostat: ',fstat
                end if
             end if
             close(123)
          else
             if (random_verbosity>0) then
                write(*,*) 'note: failed to open /dev/urandom'
                write(*,*) 'iostat: ',fstat
             end if
          end if
       end if
    end if
    if (.not.seedok) then
!!! does not work on the PGI compiler
#ifndef __PGI
       call date_and_time(values=datetime)
       seed = 100*datetime(7) + datetime(8)/10
       call random_seed(put=seed)
       if (random_verbosity>0) &
            write(*,*) 'note: used system date to initialize the seed'
       seedok = .true.
#endif
    end if

    if (.not.seedok) then
       if (random_verbosity>=0) &
            write(*,*) 'note: used the system default to set the seed'
    end if

    !! save the initial seed to seedfile0
    if (save_initial_seed /= 0) then
       open(123,file=seedfile0,status="replace",iostat=fstat)
       write(*,*) seed
       close(123)
    end if

    !! the internal int32seed 
    ! call random_number(u)
    ! int32seed = int((u-1.0_dbl)*2**16)
    int32seed = int(seed(1))
    deallocate(seed)

    inited=.true.

  end subroutine random_initialize


!!! random end-of-job
  subroutine random_eoj()
    implicit none
    integer, allocatable :: seed(:)
    integer :: n, fstat

    !!! save random seed
    if (save_final_seed /= 0) then
       call random_seed(size=n)
       allocate(seed(n))
       call random_seed(get=seed(1:n))
       open(123,file=seedfile,status="unknown",iostat=fstat)
       if ( fstat .eq. 0) then
          if (random_verbosity>0) &
               write(*,*) 'note: saving random seed to '//trim(seedfile)
          write(123,*) seed(1:n)
          close(123)
       else
          if (random_verbosity>=0) &
               write(*,*) 'warning: could not save seed to '//trim(seedfile)
       end if
       deallocate(seed)
    end if

  !  inited=.false. ! not really needed to set this to false

  end subroutine random_eoj

end module mcmcrand
