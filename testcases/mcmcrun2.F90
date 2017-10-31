!!! ------------------------------------------------------------------------
!!! Example that uses its own initialize 
!!!
!!! Marko Laine <marko.laine@fmi.fi>
!!! ------------------------------------------------------------------------

!!!
!!! Main program
!!! 
program mcmcmain
  implicit none
  !! call mcmc_main in the library libmcmcrun.a
  call mcmc_main()
end program mcmcmain


!!!
!!! -2*log(p(obs|params)) for
!!! general Gaussian likelihood, data are read from file
!!!
function ssfunction(theta,npar,ny) result(ss)

  use mcmcprec
  use matutils, only : loaddata

  implicit none
  integer(kind=ik4) :: npar, ny
  real(kind=dbl) :: theta(npar)
  real(kind=dbl) :: ss(ny)

  real(kind=dbl), save, pointer :: data(:,:)
  logical, save :: first = .true.

  interface
     function modelfunction(theta,x) result(y)
       use mcmcprec
       implicit none
       real(kind=dbl) :: theta(:), x(:), y(size(x))
     end function modelfunction
  end interface

  !! during the first call, data is loaded
  if (first) then

     call loaddata('data.dat',data)
     first = .false.

  endif
 
  !! calculate the actual sum of squares
  ss(1) = sum((data(:,2) - modelfunction(theta,data(:,1)))**2)

end function ssfunction

!!!
!!! model function used by ssfunction
!!!
function modelfunction(theta,x) result(y)

  use mcmcprec
  implicit none
  real(kind=dbl) :: theta(:), x(:), y(size(x))

  if (size(theta) .ne. 2) stop 'something wrong with sizes in modelfunction'
  !! simple exponential model
  y = theta(1)*exp(-theta(2)*x)

end function modelfunction

!!!
!!! this function returns false if any theta(i) is out of bounds
!!!
function checkbounds(theta)
  implicit none
  real*8 theta(:)
  logical checkbounds
  
!! example: all thetas must be positive
  checkbounds = .true.
  if (any(theta<=0.0)) checkbounds = .false.
  return

end function checkbounds

subroutine initialize(par0,npar,cmat0,initcmatn,sigma2,nobs,nycol)

  use mcmcprec
  use matutils, only : loaddata, readdata, sizecheck_mat

  implicit none
  
  integer, intent(inout) :: npar, initcmatn, nycol
  integer, intent(inout), allocatable :: nobs(:)
  real(kind=dbl), intent(inout), allocatable :: par0(:), cmat0(:,:)
  real(kind=dbl), intent(inout), allocatable :: sigma2(:)

  !! loaddata needs pointers
  real(kind=dbl), pointer :: cmat(:,:), xpar(:)

  interface 
     subroutine p2a(p,x)
       use mcmcprec
       implicit none
       real(kind=dbl), intent(inout), pointer :: p(:,:)
       real(kind=dbl), intent(out), allocatable :: x(:,:)
     end subroutine p2a
     subroutine p2av(p,x)
       use mcmcprec
       implicit none
       real(kind=dbl), intent(inout), pointer :: p(:)
       real(kind=dbl), intent(out), allocatable :: x(:)
     end subroutine p2av
     subroutine load(f,x, n1, n2)
       use mcmcprec
       use matutils
       implicit none
       character(len=*), intent(in) :: f
       real(kind=dbl), intent(out), allocatable :: x(:,:)
       integer, intent(in), optional :: n1, n2
     end subroutine load
     subroutine load2(f,x, n)
       use mcmcprec
       use matutils
       implicit none
       character(len=*), intent(in) :: f
       real(kind=dbl), intent(out), allocatable :: x(:)
       integer, intent(in), optional :: n
     end subroutine load2
  end interface

  nycol = 1

  !! par0
!  call load2('mcmcpar.dat',par0)
  call loaddata('mcmcpar.dat',xpar)
  npar = size(xpar)
  if (npar < 1) then
     write(*,*) 'ERROR: npar < 1!'
     stop
  end if  
  write(*,*) 'note: read mcmcpar.dat, npar =',npar
!  allocate(par0(npar)); par0 = xpar; deallocate(xpar)
  call p2av(xpar,par0)

  !! cmat0
!  call load('mcmccov.dat',cmat0, npar, npar)
  call loaddata('mcmccov.dat',cmat)
  call sizecheck_mat(cmat,npar,npar,'mcmccov.dat')
!  allocate(cmat0(npar,npar));  cmat0 = cmat; deallocate(cmat)
  call p2a(cmat,cmat0)

  !! sigma2 nobs
  allocate(sigma2(nycol))
  allocate(nobs(nycol))

  sigma2 = 0.5
  nobs = 11

  !! initcmatn
  if (initcmatn<0) initcmatn = 0

end subroutine initialize


!!! allocates x, copies pointer p to it and deallocates p
subroutine p2a(p,x)
  use mcmcprec
  implicit none
  real(kind=dbl), intent(inout), pointer :: p(:,:)
  real(kind=dbl), intent(out), allocatable :: x(:,:)

  allocate(x(size(p,1),size(p,2)))
  x = p ! copy p to x
  deallocate(p)
end subroutine p2a

subroutine p2av(p,x)
  use mcmcprec
  implicit none
  real(kind=dbl), intent(inout), pointer :: p(:)
  real(kind=dbl), intent(out), allocatable :: x(:)

  allocate(x(size(p)))
  x = p ! copy p to x
  deallocate(p)
end subroutine p2av

subroutine load(f,x, n1, n2)
  use mcmcprec
  use matutils
  implicit none
  character(len=*), intent(in) :: f
  real(kind=dbl), intent(out), allocatable :: x(:,:)
  integer, intent(in), optional :: n1, n2

  real(kind=dbl), pointer :: p(:,:)
  
  call loaddata(f,p)
!  call p2a(p,x)
  allocate(x(size(p,1),size(p,2)))
  x = p ! copy p to x
  deallocate(p)

  if (present(n1) .and. present(n2)) then
     call sizecheck_mat(x,n1,n2)
  end if
  
end subroutine load

subroutine load2(f,x, n1)
  use mcmcprec
  use matutils
  implicit none
  character(len=*), intent(in) :: f
  real(kind=dbl), intent(out), allocatable :: x(:)
  integer, intent(in), optional :: n1

  real(kind=dbl), pointer :: p(:)
  
  call loaddata(f,p)
!  call p2a(p,x)
  allocate(x(size(p)))
  x = p ! copy p to x
  deallocate(p)

  if (present(n1)) then
!     call sizecheck_vec(x,n1)
     if (size(x).ne.n1) then
        write(*,*) 'sizes do not match, while reading file ',trim(f)
     end if
  end if
  
end subroutine load2
