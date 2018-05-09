!!! $Id: matutils.F90,v 1.49 2017/08/10 07:18:33 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! this file is part of mcmc library
!!! File: matutils.F90
!!! Purpose: utility subroutines for the mcmc module
!!!     matrix utilities
!!!     error messages
!!!     allocate utils
!!!     file load and save
!!!
!!! Calls BLAS and LAPACK for matrix and vector operations, and
!!! auxiliary files dchud.f dchdd.f (from LINPACK) for Cholesky
!!! up/downdate
!!!
!!! Marko Laine 2001 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!!
module matutils

!#define EPPES 1  ! uncomment this if compiling the EPPES version
#define FILEUNIT 1234
  
  use mcmcprec    ! precision and machine constants

  implicit none

  character(len=*), parameter :: utils_version="$Revision: 1.49 $ Compiled: "//__DATE__

  include 'lapack_inc.h'

  interface writedata
     module procedure writedata_mat, writedata_vec, writedata_scal
  end interface
  interface readdata
     module procedure readdata_vec, readdata_mat, readdata_n, readdata_x
  end interface
  interface reallocate
     module procedure reallocate_rv, reallocate_rm
  end interface
  interface loaddata
     module procedure loaddata_mat, loaddata_vec, loaddata_x, &
          loaddata_int, loaddata_ints
  end interface
  interface loaddata2
     module procedure loaddata_mata, loaddata_veca
  end interface
  interface printmat
     module procedure printmat_mat, printmat_vec, printmat_ivec
  end interface


!!! A problem with output format here, the old G20.11 produces numbers like 0.12323200000-100
!!! with missing "E" is 3 numbers in the exponent, on th eother hand some older gfortran (4.3.4)
!!! does not allow for G0
#ifdef __GFORTRAN__
  character(len=*), parameter :: gfrmstr = 'G0'
#else
  character(len=*), parameter :: gfrmstr = 'G0'
#endif
  character(len=526), save :: msgbuff

  integer, save :: UtilInfoLevel = 0
  integer, save :: UtilErrorLevel = 0
  integer, parameter :: InfoLevelMinimal = 0
  integer, parameter :: InfoLevelDebug = 2
  integer, parameter :: InfoLevelInformational = 3

  private :: gfrmstr, UtilInfoLevel, UtilErrorLevel, msgbuff

contains

!!! we need these matrix multiplication utilities to use the Cholesky
!!! factor provided by LAPACK routine, where C=R'R, R upper diagonal
!!! and the lower diagonal entries are arbitrary.

!!! TODO: general matmulu(A,B, trans) for AB A'B and AX, also need A'BA
!!! with B matrix or vector

!!!
!!! vector matrix multiplication
!!! p = x'*v (default) or p = x*v , x upper diagonal
!!!
  function matmulu(x,v, trans) result(p)
    implicit none
    real(kind=dbl), intent(in) :: v(:)
    real(kind=dbl), intent(in) :: x(:,:)
    character(len=1), intent(in), optional :: trans
    real(kind=dbl), dimension(size(v)) :: p
    
    integer :: n,m
    character(len=1) t
    
    t = 'T'
    if (present(trans)) then
       if (trans == 'T' .or. trans == 't') then
          t=trans
       else
          t='n'
       end if
    end if

    n = size(v,1)
    m = size(x,1)
    if (n .ne. m) then
       stop 'matmulu: incorrect dimensions'
    end if
    call dcopy(n,v,1,p,1) ! copy v to p p = v
    call dtrmv('u',t,'n',n,x,n,p,1)
  end function matmulu
!!!
!!! same but with lower diagonal x
!!!
  function matmull(x,v) result(p)
    implicit none
    real(kind=dbl), intent(in) :: v(:)
    real(kind=dbl), intent(in) :: x(:,:)
    real(kind=dbl), dimension(size(v)) :: p

    integer :: n,m

    n = size(v,1)
    m = size(x,1)
    if (n .ne. m) then
       stop 'matmull: incorrect dimensions'
    end if
!   call dsymv('l',n,1.0d0,x,n,v,1,0.0d0,p,1)

    call dcopy(n,v,1,p,1) ! copy v to p p = v
    call dtrmv('l','n','n',n,x,n,p,1)

  end function matmull
!!!
!!! same but with general x (n * n)
!!! p = x*v (default) or x = x'*v
!!!
  function matmulx(x,v, trans) result(p)
    implicit none
    real(kind=dbl), intent(in) :: v(:)
    real(kind=dbl), intent(in) :: x(:,:)
    character(len=1), intent(in), optional :: trans
    real(kind=dbl), dimension(size(v)) :: p
    
    integer :: n,m
    character(len=1) t
    
    t = 'N'
    if (present(trans)) then
       if (trans == 'T' .or. trans == 't') then
          t=trans
       else
          t='n'
       end if
    end if

    n = size(v,1)
    m = size(x,1)
    if (n .ne. m) then
       stop 'matmulx: incorrect dimensions'
    end if
    call dgemv(t,n,n,1.0d0,x,n,v,1,0.0d0,p,1)
  end function matmulx
!!!
!!! vector matrix multiplication
!!! p = x*v, x symmetric, using upper part
!!!
  function matmuls(x,v) result(p)
    implicit none
    real(kind=dbl), intent(in) :: v(:)
    real(kind=dbl), intent(in) :: x(:,:)
    real(kind=dbl), dimension(size(v)) :: p
    
    integer :: n,m
    
    n = size(v,1)
    m = size(x,1)
    if (n .ne. m) then
       stop 'matmuls: incorrect dimensions'
    end if
    call dsymv('u',n,1.0d0,x,n,v,1,0.0d0,p,1)
  end function matmuls

!!!
!!! matrix matrix multiplication
!!! p = x*y, x upper diagonal, using upper part
!!!
!!!!! DO NOT USE YET !!!!
  function matmulxy(x,y,side0,trans0) result(p)
    implicit none
    real(kind=dbl), intent(in) :: y(:,:)
    real(kind=dbl), intent(in) :: x(:,:)
    character(len=1), intent(in), optional :: side0, trans0
    real(kind=dbl), dimension(size(x,1),size(y,2)) :: p
    
    integer :: mx,nx,my,ny
    character(len=1) side, trans

    side = 'L'
    if (present(side0)) then
       if (side0 == 'R' .or. side0 == 'r') then
          side=side0
       end if
    end if
    trans = 'n'
    if (present(trans0)) then
       if (trans0 == 'T' .or. trans0 == 't') then
          trans=trans0
       end if
    end if

    
    mx = size(x,1)
    nx = size(x,2)
    my = size(y,1)
    ny = size(y,2)
    if (nx .ne. my) then
       stop 'matmulxy: incorrect dimensions'
    end if
!    call dsymm(side,'u',mx,ny,1.0d0,x,mx,y,my,0.0d0,p,mx)
! ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
    p = y
    call dtrmm(side,'u',trans,'n',my,ny,1.0d0,x,mx,p,my)
!( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB )
  end function matmulxy

!!!
!!! calculates covariance matrix of x,
!!! with optional (frequency) weights in w,
!!! returns cmat and optionally mean and sum of weights in xmean and wsum
!!! updates cmat, cmean and wsum, if update=.true. and wsum>0
!!! 
  subroutine covmat(x, cmat, w, xmean, wsum, update)
    implicit none
    real(kind=dbl), intent(in) :: x(:,:)
    real(kind=dbl), intent(inout) :: cmat(:,:)
    real(kind=dbl), intent(in), optional :: w(:)
    real(kind=dbl), intent(inout), optional :: xmean(size(cmat,2))
    real(kind=dbl), intent(inout), optional :: wsum
    logical, intent(in), optional :: update

    ! local variables
    real(kind=dbl) :: xmean2(size(cmat,2))
    integer :: n,p,i,j
    real(kind=dbl) :: wsum2, w2, w3
    logical :: doupdate

    n=size(x,1)
    p=size(x,2)

    if (size(cmat,1) /= p .or. size(cmat,2) /= p) then
       stop 'covmat: invalid cmat'
    end if

    !! optional weight in w
    if (present(w) .and. size(w) == n) then
       w2 = -1.0e0_dbl
       wsum2 = sum(w)
    else if (present(w) .and. size(w) == 1) then
       w2 = w(1)
       wsum2 = dble(n)*w2
    else if (.not. present(w)) then
       w2 = 1.0e0_dbl
       wsum2 = dble(n)
    else
       stop 'covmat: invalid weight in w'
    end if

    !!  update
    doupdate=.false.
    if (present(update) .and. present(wsum)) then
       if (update .and. (wsum > 0.0_dbl)) then
          doupdate=.true.
          if (.not. present(xmean)) then
             stop 'invalid input in covmat'
          end if
       else
          doupdate=.false.
       end if
    end if
    if (doupdate) then
       !! covariance and mean update
       !! one line update at a time
       do i=1,n
          xmean2 = x(i,:)-xmean ! center x and store it in xmean2
          if (w2 == -1.0e0_dbl) then
             w3 = w(i)
          else
             w3 = w2
          end if
!          cmat = cmat + w3/(wsum+w3-1.0d0) &
!               * (wsum/(wsum+w3) &
!               * matmul(reshape(xmean2,(/p,1/)), &
!               reshape(xmean2,(/1,p/))) &
          !               - cmat)
          
          !! new version, as matmul(reshape(x,(/p,1/)),reshape(x,(/1,p/)))
          !! seems to generate run time error in gfortran 8.1.0
          cmat = cmat + w3/(wsum+w3-1.0d0) &
               * (wsum/(wsum+w3) &
               * outer_product(xmean2,xmean2) &
               - cmat)
          
          !! or simplify to this (needs test)
!          cmat = (wsum-1.0d0)/(w3+wsum-1.0d0)*cmat + &
!               w3*wsum/(w3+wsum-1.0d0)/(w3+wsum) &
!               * matmul(reshape(xmean2,(/p,1/)), &
!               reshape(xmean2,(/1,p/)))
          xmean = xmean + w3/(wsum+w3)*xmean2
          wsum = w3 + wsum
       end do
    else !!! no update
       do i=1,p
          if (w2 == -1.0e0_dbl) then
             xmean2(i) = sum(x(:,i)*w)/wsum2
          else
             xmean2(i) = sum(x(:,i)*w2)/wsum2
          end if
       end do
       do i=1,p
          do j=1,i
             if ( w2 == -1.0e0_dbl ) then
                cmat(i,j) = dot_product( (x(:,i)-xmean2(i)), &
                     (x(:,j)-xmean2(j))*w) /  (wsum2 - 1.0_dbl)
             else
                cmat(i,j) = dot_product( (x(:,i)-xmean2(i)), &
                     (x(:,j)-xmean2(j))*w2) /  (wsum2 - 1.0_dbl)
             end if
             if (i/=j) cmat(j,i) = cmat(i,j)
          end do
       end do

       if (present(xmean)) then
          xmean=xmean2
       end if
       if (present(wsum)) then
          wsum=wsum2
       end if
       
    end if ! update

  end subroutine covmat
!!!
!!! Calculates upper diagonal Cholesky transformation
!!! of symmetric matrix 'cmat'
  subroutine covtor(cmat,R,info)
    implicit none
    real(kind=dbl) :: cmat(:,:)
    real(kind=dbl), optional :: R(:,:)
    integer, intent(out), optional :: info

    integer :: info2, n

    n = size(cmat,1)
    if (size(cmat,2) /= n) then
       stop 'size error in covtor'
    end if

    if (present(R)) then
       if (size(R,1) /= n .or. size(R,1) /= size(R,2)) then
          stop 'size error in R in covtor'
       end if
       R = cmat
       call dpotrf('u',n,R,n,info2)
    else
       call dpotrf('u',n,cmat,n,info2)
    end if

    if (present(info)) then
       info=info2
    else  if (info2 .ne. 0) then
       stop 'error in Chol'
    end if

  end subroutine covtor
!!!
!!! test for svd sqrt
!!!
  subroutine covtor_svd(cmat,R,condmax,info)
    implicit none
    real(kind=dbl) :: cmat(:,:)
    real(kind=dbl) :: R(:,:)
    real(kind=dbl), intent(in) :: condmax
    integer, intent(out), optional :: info

    integer :: info2, n, i
    real(kind=dbl), allocatable :: work(:)
    real(kind=dbl), allocatable :: u(:,:), s(:)
    integer :: lwork

    real(kind=dbl) :: tol, cond

    n = size(cmat,1)
    if (size(cmat,2) /= n) then
       stop 'size error in covtor'
    end if

    lwork = 5*n
    allocate(work(lwork),u(n,n),s(n), stat=info2)
    if (info2 /= 0) then
       stop 'memory allocation in svd'
    end if

    if (size(R,1) /= n .or. size(R,1) /= size(R,2)) then
       stop 'size error in R in covtor'
    end if
    R = cmat

!!! call dgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info)
    call dgesvd('A','N',n,n,R,n,s,u,n,u,n,work,lwork,info2)
!   write(*,*) 'S       : ', real(s(max(1,n-3):n))
    if (s(1) == 0.0d0) then
       write(*,*) 'SVD: s(1)=0'
       if (present(info)) then
          info=n
       end if
       return
    end if
    if (info2 > 0) then
       write(*,*) 'svd info = ', info2
    end if
    tol = dble(n)*s(1)*2.2204d-16 !!!! check
    tol = dble(n)*s(1)*1.0d-6
    cond = condmax
    tol = s(1)/cond
    if (s(n) <= tol) then
       write(*,*) 'adjusting svd, tol:',real(tol)
       do i=1,n
          if (s(i) < tol) then
             s(i) = tol
          end if
       end do
       info2 = -1
    end if

!   R = u*sqrt(s)
    R = u

    !! scale R so that R'*R = cmat (or R*R' check!)
    !! for scam we would need to have the unscaled version!
    do i=1,n
!!!    R(1:n,i) = sqrt(s(i))*R(1:n,i)
       call dscal(n,sqrt(s(i)),R(1:n,i),1)
    end do

    if (present(info)) then
       info=info2
    else  if (info2 < 0 .or. info2 > n-2) then
       write(*,*) 'svd info = ', info2, ' Stopping'
       stop
    end if
    deallocate(work,u,s)

  end subroutine covtor_svd
!!!
!!! invert symmetric matrix using Cholesky decomposition
!!!
  function invertmat_sym(x,info0) result(y)
    implicit none
    real(kind=dbl), intent(in) :: x(:,:)
    integer, optional, intent(inout) :: info0
    real(kind=dbl) :: y(size(x,1),size(x,2))

    integer(kind=ik4) :: npar, info

    npar = size(x,1)
    if (npar .ne. size(x,2)) then
       write(*,*) 'ERROR: x must be square'
       stop
    end if
    y = x
    call dpotrf('u',npar,y,npar,info)
    if (info .ne. 0) then
       if (present(info0)) then
          info0 = info
          return
       else
          write(*,*) 'ERROR: cannot factor cmat, info = ', info
          stop
       end if
    end if
    call dpotri('u',npar,y,npar,info)
    if (info .ne. 0) then
       if (present(info0)) then
          info0 = info
          return
       else
          write(*,*) 'ERROR: cannot invert cmat, info = ', info
          stop
       end if
    end if

    !! fill the lower diagonal with the upper
    call copylower(y)
    if (present(info0)) info0=0

  end function invertmat_sym
!!!
!!! invert symmetric matrix using dsyrti
!!!
  function invertmat_sym2(x,info0) result(y)
    implicit none
    real(kind=dbl), intent(in) :: x(:,:)
    integer, optional, intent(inout) :: info0
    real(kind=dbl) :: y(size(x,1),size(x,2))
    
    real(kind=dbl) :: work0(1)
    real(kind=dbl),allocatable :: work(:)
    integer :: ipiv(size(x,1))
    integer(kind=ik4) :: npar, lwork, info

    npar = size(x,1)
    if (npar .ne. size(x,2)) then
       write(*,*) 'ERROR: x must be square'
       stop
    end if
    y = x
!    call dpotrf('u',npar,y,npar,info)
    call dsytrf('u', npar, y, npar, ipiv, work0, -1, info)
    lwork = max(npar,int(work0(1)))
    allocate(work(lwork))
    call dsytrf('u', npar, y, npar, ipiv, work, lwork, info)
    if (info .ne. 0) then
       if (present(info0)) then
          info0 = info
          deallocate(work)
          return
       else
          write(*,*) 'ERROR: cannot factor cmat, info = ', info
          stop
       end if
    end if
!    call dpotri('u',npar,y,npar,info)
    call dsytri('u', npar, y, npar, ipiv, work, info)
    deallocate(work)
    if (info .ne. 0) then
       if (present(info0)) then
          info0 = info
          return
       else
          write(*,*) 'ERROR: cannot invert cmat, info = ', info
          stop
       end if
    end if

    !! fill the lower diagonal with the upper
    call copylower(y)
    if (present(info0)) info0=0

  end function invertmat_sym2
!!!
!!! outer product of two vectors
!!!
  function outer_product(x1,x2) result(y)
    implicit none
    real(kind=dbl), intent(in) :: x1(:), x2(:)
    real(kind=dbl) :: y(size(x1),size(x2))

    y = spread(x1,2,size(x2))*spread(x2,1,size(x1))

  end function outer_product
!!!
!!! Mahalanobis distance of rows of x from mu given cmat
  function mahalanobis(x,mu,cmat) result(d)
    implicit none
    real(kind=dbl) :: x(:,:), mu(:), cmat(:,:), d(size(x,1))
    real(kind=dbl) :: cinv(size(cmat,1),size(cmat,2)), x0(size(x,1),size(x,2))
    integer :: i, n

    n = size(x,1)
    ! remove mean
    do i=1,n
       x0(i,:) = x(i,:)-mu
    end do
    ! invert cmat
    cinv=invertmat_sym(cmat)
    ! Mahalanobis distance for each row in x
    d = sum(matmul(x0,cinv)*x0,2)

  end function mahalanobis
!!!
!!! test for svd sqrt
!!!
  subroutine scam_svd(cmat,R,std,condmax,info)
    implicit none
    real(kind=dbl) :: cmat(:,:)
    real(kind=dbl) :: R(:,:)
    real(kind=dbl) :: std(size(R,1))
    real(kind=dbl), intent(in) :: condmax
    integer, intent(out), optional :: info

    integer :: info2, n, i
    real(kind=dbl), allocatable :: work(:)
    real(kind=dbl), allocatable :: u(:,:), s(:)
    integer :: lwork

    real(kind=dbl) :: tol, cond

    n = size(cmat,1)
    if (size(cmat,2) /= n) then
       stop 'size error in covtor'
    end if

    lwork = 5*n
    allocate(work(lwork),u(n,n),s(n), stat=info2)
    if (info2 /= 0) then
       stop 'memory allocation in svd'
    end if

    if (size(R,1) /= n .or. size(R,1) /= size(R,2)) then
       stop 'size error in R in covtor'
    end if
    R = cmat

!!! call dgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info)
    call dgesvd('A','N',n,n,R,n,s,u,n,u,n,work,lwork,info2)
!   write(*,*) 'S       : ', real(s(max(1,n-3):n))
    if (s(1) == 0.0d0) then
       write(*,*) 'SVD: s(1)=0'
       if (present(info)) then
          info=n
       end if
       return
    end if
    if (info2 > 0) then
       write(*,*) 'svd info = ', info2
    end if
    tol = dble(n)*s(1)*2.2204d-16 !!!! check
!   tol = dble(n)*s(1)*1.0d-6
    cond = condmax
    tol = s(1)/cond
    if (s(n) <= tol) then
       write(*,*) 'adjusting svd, tol:',real(tol)
       do i=1,n
          if (s(i) < tol) then
             s(i) = tol
          end if
       end do
       info2 = -1
    end if

!   R = u*sqrt(s)
    R = u
    std = sqrt(s)

    if (present(info)) then
       info=info2
    else  if (info2 < 0 .or. info2 > n-2) then
       write(*,*) 'svd info = ', info2, ' Stopping'
       stop
    end if
    deallocate(work,u,s)

  end subroutine scam_svd
!!!
!!! Cholesky update using dchud from Linpack
!!! 
!!! for the EPPES version we do not want to include dchud.f
#ifndef EPPES
subroutine cholupdate(r,x,info)
  implicit none
  real(kind=dbl), intent(inout) :: r(:,:)
  real(kind=dbl), intent(in) :: x(:)
  integer, intent(out), optional :: info

  integer ldr,p,ldz,nz
  real(kind=dbl) rho(1), c(size(r,2))
  real(kind=dbl) z(1,1),y(1),s(size(r,2))
  
  interface
     subroutine dchud(r,ldr,p,x,z,ldz,nz,y,rho,c,s)
       implicit none
       integer ldr,p,ldz,nz
       double precision rho(nz),c(p)
       double precision r(ldr,p),x(p),z(ldz,nz),y(nz),s(p)
     end subroutine dchud
  end interface

  ldr = size(r,1)
  p   = size(r,2)
  ldz = 0
  nz  = 0
  call dchud(r,ldr,p,x,z,ldz,nz,y,rho,c,s)
  if (present(info)) info = 0

end subroutine cholupdate
#endif
!!!
!!! Cholesky downdate using dchdd from Linpack
!!! 
#ifndef EPPES
subroutine choldowndate(r,x,info)
  implicit none
  real(kind=dbl), intent(inout) :: r(:,:)
  real(kind=dbl), intent(in) :: x(:)
  integer, intent(out), optional :: info

  integer ldr,p,ldz,nz
  real(kind=dbl) rho(1), c(size(r,2))
  real(kind=dbl) z(1,1),y(1),s(size(r,2))

  integer info0

  interface
     subroutine dchdd(r,ldr,p,x,z,ldz,nz,y,rho,c,s,info)
       implicit none
       integer ldr,p,ldz,nz,info
       double precision r(ldr,p),x(p),z(ldz,nz),y(nz),s(p)
       double precision rho(nz),c(p)
     end subroutine dchdd
  end interface

  ldr = size(r,1)
  p   = size(r,2)
  ldz = 0
  nz  = 0
  call dchdd(r,ldr,p,x,z,ldz,nz,y,rho,c,s,info0)
  if (present(info)) then 
     info = info0
  elseif (info0.ne.0) then
     write(*,*) 'error in coldowndate, info:', info
     stop
  end if

end subroutine choldowndate
#endif
!!!

!!!
!!! error routines NOT IN USE YET
!!!
!!! SetMessageLevel - set doerror and message levels
!!!
  subroutine SetMessageLevel(infolevel,errorlevel)
    implicit none
    integer, intent(in), optional :: infolevel
    integer, intent(in), optional :: errorlevel

    if (present(infolevel)) then
       UtilInfoLevel = max(0,infolevel)
    end if
    if (present(errorlevel)) then
       UtilErrorLevel = max(0,errorlevel)
    end if
    
  end subroutine SetMessageLevel
!!!
!!!GetInfoLevel - return current InfoLevel
!!!
  integer function getinfolevel()
    implicit none
    getinfolevel = UtilInfoLevel
  end function getinfolevel
!!!
!!!GetErrorLevel - return current ErrorLevel
!!!
  integer function geterrorlevel()
    implicit none
    geterrorlevel = UtilErrorLevel
  end function geterrorlevel
!!!
!!! DoError - error message and action
!!! action = 0, stop
!!! action = 1, 
  subroutine doerror(creason,elevel,action)
    implicit none
    character(len=*) , intent(in) :: creason ! character mesassage
    integer, intent(in), optional :: elevel ! error level
    integer, intent(in), optional :: action ! error action
    integer :: act, el
    
    if (present(elevel)) then
       el = elevel
    else
       el = 0
    end if
    if (present(action)) then
       act = action
    else
       act = 0
    end if

    if (el >= UtilErrorLevel) then
       write(*,*) 'ERROR: ',trim(creason)
       if (act .eq. 0) then
          stop 'error stop'
       end if
    end if
    return
  end subroutine doerror
!!!
!!! Message - send text to console
!!!
  subroutine message(text,ilevel)
    implicit none
    character(len=*), intent(in) :: text
    integer, intent(in), optional :: ilevel

    integer :: il

    if(.not.present(ilevel)) then
       il = 0
    else
       il = ilevel
    end if

    if (il >= UtilInfoLevel) then
       write(*,*) text
    end if

  end subroutine message
!!! note routines
  subroutine note_txt(verb,txt)
    implicit none
    integer, intent(in) :: verb
    character(len=*), intent(in) :: txt

    write(msgbuff,*) 'note: ',trim(txt)
    call message(msgbuff,verb)    
  end subroutine note_txt
  subroutine note_x(verb,txt,x)
    implicit none
    integer, intent(in) :: verb
    character(len=*), intent(in) :: txt
    real(kind=dbl), intent(in) :: x

    write(msgbuff,*) 'note: ',trim(txt),' ',real(x)
    call message(msgbuff,verb)    
  end subroutine note_x
  subroutine note_c(verb,txt,c)
    implicit none
    integer, intent(in) :: verb
    character(len=*), intent(in) :: txt
    character(len=*), intent(in) :: c

    write(msgbuff,*) 'note: ',trim(txt),' ',trim(c)
    call message(msgbuff,verb)    
  end subroutine note_c
!!!
!!! write data matrix to file, with optional locking
!!! 
  subroutine writedata_mat(file,xmat,stat,uselock)
    implicit none
    character(len=*), intent(in) :: file
    real(kind=dbl), intent(in) :: xmat(:,:)
    integer, intent(out), optional :: stat
    logical, intent(in), optional :: uselock
    
    integer :: m,n, fstat, i, j
    logical :: uselock2
    
    if (len_trim(file) == 0) then
       if (present(stat)) stat=-1
       return
    end if

    uselock2 = .false.
    if (present(uselock)) uselock2 = uselock
    
    m = size(xmat,1)
    n = size(xmat,2)
    
    if (n < 1 .or. m < 1) then
       write(*,*) 'Error in writedata, empty matrix, file:',trim(file)
       return
    end if

    if (uselock2) then
       call open_with_lock(file=file,unit=FILEUNIT,status='replace',iostat=fstat)
    else
#ifdef _WIN32
       open(unit=FILEUNIT, file=file, status='replace', &
            recordtype='stream_cr', iostat=fstat)
#else
       open(unit=FILEUNIT, file=file, status='replace',  iostat=fstat)
#endif
    end if
    if (fstat /= 0) then
       if (present(stat)) then
          stat = fstat
          return
       else
          write(*,*) 'Error opening file ', trim(file)
          stop
       end if
    end if
    do i=1,m
       do j=1,n
!         write(FILEUNIT,'('//gfrmstr//',x)', iostat=fstat, advance='NO') real(xmat(i,j))
          write(FILEUNIT,'('//gfrmstr//',x)', iostat=fstat, advance='NO') xmat(i,j)
          if (fstat /= 0) then 
             write(*,*) 'Write error on file ', trim(file)
             stop 
          end if
       end do
#ifdef __PGI
       !! PGI fortran adds extra newline at eof
       if (i<m) write(FILEUNIT,*)
#else
       write(FILEUNIT,*)
#endif
    end do
    call close_with_lock(unit=FILEUNIT,file=file,iostat=fstat,uselock=uselock2)

    if (present(stat)) stat=fstat

    return
  end subroutine writedata_mat
!!!
!!! write vector to file
!!!
  subroutine writedata_vec(file,xvec,stat)
    implicit none
    character(len=*), intent(in) :: file
    real(kind=dbl), intent(in) :: xvec(:)
    integer, intent(out), optional :: stat

!    real(kind=dbl) :: xmat(size(xvec),1)
    integer :: i, n, fstat

    if (len_trim(file) == 0) then
       if (present(stat)) stat=-1
       return
    end if

    n = size(xvec)
!   xmat = reshape(xvec,(/n,1/)) !! stack overflow to happen here    
!   call writedata(file,xmat,stat)

#ifdef _WIN32
    open(unit=FILEUNIT, file=file, status='replace', &
         recordtype='stream_cr', iostat=fstat)
#else
    open(unit=FILEUNIT, file=file, status='replace', iostat=fstat)
#endif
    if (fstat /= 0) then
       if (present(stat)) then
          stat = fstat
          return
       else
          write(*,*) 'Error opening file ', trim(file)
          stop
       end if
    end if

    do i=1,n
       write(FILEUNIT,'('//gfrmstr//')', iostat=fstat, advance='YES') xvec(i)
       if (fstat /= 0) then 
          write(*,*) 'Write error on file ', trim(file)
          stop 
       end if
    end do
    close(unit=FILEUNIT)

    if (present(stat)) stat=fstat

  end subroutine writedata_vec
!!!
!!! write dbl scalar to ascii file
!!! 
  subroutine writedata_scal(file,x,stat)
    implicit none
    character(len=*), intent(in) :: file
    real(kind=dbl), intent(in) :: x
    integer, intent(out), optional :: stat

    real(kind=dbl)  :: xvec(1)
    
    xvec(1) = x

    call writedata(file,xvec,stat)

  end subroutine writedata_scal
!!!
!!! write vector x to unit and try to avoid a newline
!!! this is not used anywhere?
  subroutine writevec(unit,x,stat)
    implicit none
    integer, intent(in) :: unit
    real(kind=dbl), intent(in) :: x(:)
    integer, intent(out), optional :: stat
    
    integer :: n, fstat, i

    n = size(x)
    do i=1,n
       write(unit,'('//gfrmstr//',x)', iostat=fstat, advance='NO') x(i)
       if (fstat /= 0) then 
          if (present(stat)) then 
             stat=fstat
             return
          else
             write(*,*) 'Write error on unit ', unit
             stop 
          end if
       end if
    end do
    write(unit,*) ! newline
    if (present(stat)) stat = 0    
    return
  end subroutine writevec
!!!
!!!
!!! Load ascii data from file to a newly allocated variable
!!! Optional file locking
!!!
!!! (this used to be readdata.c)
  subroutine loaddata_mat(file,xmat,status,uselock)
    implicit none
    character(len=*), intent(in) :: file
    real(kind=dbl), pointer :: xmat(:,:)
    integer, intent(out), optional :: status
    logical, intent(in), optional :: uselock
        

    integer  m, i, fstat, astat, nmem
    integer, parameter :: initmem=100, initmemincr=500
    integer, parameter :: MAXLINE=16384
    character(len=MAXLINE), target :: line
    character(len=*), parameter :: spaces = ' ,;'//achar(9) ! spaces
!    character(len=*), parameter :: spaces = achar(32)//achar(44)//achar(59)//achar(9) ! spaces
    character(len=*), parameter :: comments = '#%!Cc'       ! comments
    integer :: sprt, eprt
    integer :: ne, wl
    integer :: memincr
    logical :: ok, done
    logical :: uselock2

    uselock2 = .false.
    if (present(uselock)) uselock2 = uselock

    call open_with_lock(file=file,unit=FILEUNIT,status='old',iostat=fstat,uselock=uselock2)
    if (fstat /= 0) then
       if (present(status)) then
          status = fstat
          return
       else
          write(*,*) 'Error opening the file ', trim(file)
          stop
       end if
    end if

    !! read the first line
    ok = .false.
    do while (.not. ok)
       read(FILEUNIT,"(A)",iostat=fstat) line ! maybe use size=count here?
       if (fstat /= 0) then
          write(*,*) 'Error reading the first line of ', trim(file)
          if (present(status)) then
             status = fstat
             call close_with_lock(unit=FILEUNIT,file=file,uselock=uselock2)
             return
          else
             stop
          end if
       end if
       if (scan(line(1:1),comments) == 0) then ! scan for comment line
          ok=.true.
       end if
    end do

    !! scan the number of elements
    ne = 0
    wl = 0
    sprt = 0
    eprt = len_trim(line)
    do while (sprt < eprt)
       sprt=sprt+1
       if (scan(line(sprt:sprt),spaces)>0) then
          if (wl>0) ne=ne+1
          wl = 0
       else
          wl = wl+1
       end if
    end do
    if (wl>0) ne=ne+1

    !! allocate initial amount of memory
    memincr = initmemincr
    nmem = initmem
    nullify(xmat)
    xmat => reallocate(xmat,nmem,ne,astat)
    if (astat /= 0) then
       write(*,*) 'loaddata: error allocating memory'
       write(*,*) 'loaddata: nmem = ', nmem, 'ne = ',ne,' stat = ',astat
       if (present(status)) then
          status = -1
          call close_with_lock(unit=FILEUNIT,file=file,uselock=uselock2)
          return
       else
          stop
       end if
    end if    

    m = 1 ! number of data lines read
    done = .false.
    do                          ! loop over lines in the file
       read(line,*,iostat=fstat) (xmat(m,i), i = 1,ne)
       if (fstat < 0) then  !! eof
          write(*,*) 'eof while reading file ',trim(file)
          done = .true.
          call close_with_lock(unit=FILEUNIT,file=file,uselock=uselock2)
          if (present(status)) then
             status=fstat
             deallocate(xmat)
             return
          else
             stop
          end if
       else if (fstat > 0) then
          write(*,*) 'error reading file ',trim(file)
          call close_with_lock(unit=FILEUNIT,file=file,uselock=uselock2)
          if (present(status)) then
             status = fstat
             deallocate(xmat)
             return
          else
             stop
          end if
       end if
       !! read new line
       ok = .false.
       do while (.not. ok)
          read(FILEUNIT,"(A)",iostat=fstat) line
          if (fstat < 0) then
             done = .true.
             ok = .true.
          else if (fstat > 0) then
             write(*,*) 'Error reading line: ',m+1,' in file ', trim(file)
             if (present(status)) then
                status = fstat
                done = .true.
                ok  = .true.
             else
                stop
             end if
          else if (scan(line(1:1),comments) == 0) then
             ok=.true.
          end if
       end do
       !! are we done
       if (done) then
          call close_with_lock(unit=FILEUNIT,file=file,uselock=uselock2)
          xmat => reallocate(xmat,m,ne)
          if (present(status)) then
             status = 0
          end if
          return
       end if
       
       m = m+1
       !! do we need more memory
       if (m > nmem) then
          nmem = nmem+memincr
          xmat => reallocate(xmat,nmem,ne)
          memincr = memincr*2
       end if
    end do
  end subroutine loaddata_mat
!!! vector version
  subroutine loaddata_vec(file,xvec,status,uselock)
    implicit none
    character(len=*), intent(in) :: file
    real(kind=dbl), pointer :: xvec(:)
    integer, intent(out),optional :: status
    logical, intent(in), optional :: uselock

    real(kind=dbl), pointer :: xmat(:,:)
        
    call loaddata_mat(file,xmat,status,uselock)

    allocate(xvec(size(xmat,1)*size(xmat,2)))
    xvec = reshape(xmat,(/size(xvec)/))
    deallocate(xmat)
  end subroutine loaddata_vec
!!! scalar 
  subroutine loaddata_x(file,x,status,uselock)
    implicit none
    character(len=*), intent(in) :: file
    real(kind=dbl), intent(out) :: x
    integer, intent(out),optional :: status
    logical, intent(in), optional :: uselock

    real(kind=dbl), pointer :: xmat(:,:)
        
    call loaddata_mat(file,xmat,status,uselock)
    x = xmat(1,1)
    deallocate(xmat)
  end subroutine loaddata_x
!!! scalar integer
  subroutine loaddata_ints(file,n,status,uselock)
    implicit none
    character(len=*), intent(in) :: file
    integer, intent(out) :: n
    integer, intent(out),optional :: status
    logical, intent(in), optional :: uselock

    real(kind=dbl), pointer :: xmat(:,:)
        
    call loaddata_mat(file,xmat,status,uselock)
    n = int(xmat(1,1))
    deallocate(xmat)
  end subroutine loaddata_ints
!!! vector of integers
  subroutine loaddata_int(file,nvec,status,uselock)
    implicit none
    character(len=*), intent(in) :: file
    integer, pointer :: nvec(:)
    integer, intent(out),optional :: status
    logical, intent(in), optional :: uselock

    real(kind=dbl), pointer :: xmat(:,:)
        
    call loaddata_mat(file,xmat,status,uselock)
    allocate(nvec(size(xmat,1)*size(xmat,2)))
    nvec = reshape(int(xmat),(/size(nvec)/))
    deallocate(xmat)
  end subroutine loaddata_int

!!! xmat is allocatable, not a pointer, now we need to copy it.
!!! But due a bug/feature in gfortran, it does not distinguish this from 
!!! the pointer version.
  subroutine loaddata_mata(file,xmat,status,uselock)
    implicit none
    character(len=*), intent(in) :: file
    real(kind=dbl), intent(inout), allocatable :: xmat(:,:)
    integer, intent(out), optional :: status
    logical, intent(in), optional :: uselock

    real(kind=dbl), pointer :: xmatp(:,:)
    integer :: stat
    call loaddata_mat(file,xmatp,status,uselock)
    if (present(status)) then
       if (status /= 0) then
          return
       end if
    end if
    allocate(xmat(size(xmatp,1),size(xmatp,2)),stat=stat)
    if (stat /= 0) then
       if (present(status)) then
          status = stat
          deallocate(xmatp)
          return
       else
          call doerror('allocate error in loaddata_mata')
       end if
    end if
    xmat = xmatp
    deallocate(xmatp)
  end subroutine loaddata_mata
!!! same as mata but for a vector
  subroutine loaddata_veca(file,xvec,status,uselock)
    implicit none
    character(len=*), intent(in) :: file
    real(kind=dbl), intent(inout), allocatable :: xvec(:)
    integer, intent(out), optional :: status
    logical, intent(in), optional :: uselock

    real(kind=dbl), pointer :: xvecp(:)
    integer :: stat
    call loaddata_vec(file,xvecp,status,uselock)
    if (present(status)) then
       if (status /= 0) then
          return
       end if
    end if
    allocate(xvec(size(xvecp)),stat=stat)
    if (stat /= 0) then
       if (present(status)) then
          status = stat
          deallocate(xvecp)
          return
       else
          write(*,*) 'stat',stat
          call doerror('allocate error in loaddata_veca')
       end if
    end if
    xvec = xvecp
    deallocate(xvecp)
    if (present(status)) status = 0
  end subroutine loaddata_veca

!!!
!!!
!!! usage:  p => reallocate(p, 10000)
!!!
  function reallocate_rv(p, n, status)
    real(kind=dbl), pointer, dimension(:) :: p, reallocate_rv
    integer, intent(in) :: n
    integer, intent(out), optional :: status
    integer :: nold, ierr
    allocate(reallocate_rv(n), stat=ierr)
    if (present(status)) then
       status=ierr
    end if
    if(ierr /= 0) then 
       if (present(status)) then
          status=ierr
          return
       else
          call doerror("allocate error")
          return
        end if
    end if
    if(.not. associated(p)) return
    nold = min(size(p), n)
    reallocate_rv(1:nold) = p(1:nold)
    deallocate(p)
  end function reallocate_rv
!!!
  function reallocate_rm(p,m,n,status)
    real(kind=dbl), pointer, dimension(:,:) :: p, reallocate_rm
    integer, intent(in) :: m,n
    integer, intent(out), optional :: status
    integer :: mold,nold, ierr
    allocate(reallocate_rm(m,n), stat=ierr)
    if (present(status)) then
       status=ierr
    end if
    if(ierr /= 0) then 
       if (present(status)) then
          status=ierr
          return
        else
           call doerror("allocate error")
           return
       end if
    end if
    if(.not. associated(p)) return
    mold = min(size(p,1), m)
    nold = min(size(p,2), n)
    reallocate_rm(1:mold,1:nold) = p(1:mold,1:nold)
    deallocate(p)
  end function reallocate_rm
!!!
!!! read par from file
!!! like loaddata but par must be preallocated
!!! and data in the file must match the size of par
  subroutine readdata_vec(file,par,stat,uselock)
    implicit none
    real(kind=dbl), intent(out) :: par(:)
    character(*), intent(in) :: file
    integer, intent(out), optional :: stat
    logical, intent(in), optional :: uselock

    real(kind=dbl), pointer :: x(:,:)
    integer :: n, fstat

    n = size(par)
    fstat = 0

    call loaddata(file,x,fstat,uselock=uselock)
    if (fstat /= 0) then
       if (present(stat)) then
          stat = fstat
          return
       else         
          write(*,*) 'ERROR: Error reading file', trim(file)
          write(*,*) '       fstat = ', fstat
          write(*,*) 'Stopping now, sorry'
          stop
       end if
    end if

    if (n .ne. product(shape(x))) then
       if (present(stat)) then
          stat = -1
          deallocate(x)
          return
       else
          write(*,*) 'ERROR: size mismatch while reading file', trim(file)
          write(*,*) 'Stopping now, sorry'
          stop
       end if
    end if

    par = reshape(x,(/n/))

    deallocate(x)

  end subroutine readdata_vec
!!! subroutine readdata_mat
  subroutine readdata_mat(file,xmat,stat,uselock)
    implicit none
    real(kind=dbl), intent(out) :: xmat(:,:)
    character(*), intent(in) :: file
    integer, intent(out), optional :: stat
    logical, intent(in), optional :: uselock

    real(kind=dbl), pointer :: x(:,:)
    integer :: m,n, fstat

    m = size(xmat,1)
    n = size(xmat,2)
    fstat = 0

    call loaddata(file,x,fstat,uselock=uselock)
    if (fstat /= 0) then
       if (present(stat)) then
          stat = fstat
          return
       else         
          write(*,*) 'ERROR: Error reading file', trim(file)
          write(*,*) '       fstat = ', fstat
          write(*,*) 'Stopping now, sorry'
          stop
       end if
    end if

    if (m*n .ne.  product(shape(x))) then
       if (present(stat)) then
          stat = -1
          deallocate(x)
          return
       else
          write(*,*) 'ERROR: size mismatch while reading file', trim(file)
          write(*,*) 'want: ',shape(xmat),' is: ',shape(x)
          write(*,*) 'Stopping now, sorry'
          stop
       end if
    end if

    xmat = reshape(x,(/m,n/))

    deallocate(x)
  end subroutine readdata_mat
 !!! read a single integer
  subroutine readdata_n(file,n,stat,uselock)
    implicit none
    integer, intent(out) :: n
    character(*), intent(in) :: file
    integer, intent(out), optional :: stat
    logical, intent(in), optional :: uselock

    real(kind=dbl), pointer :: x(:,:)
    integer :: fstat

    fstat = 0
    call loaddata(file,x,fstat,uselock=uselock)
    if (fstat /= 0) then
       if (present(stat)) then
          stat = fstat
          return
       else         
          write(*,*) 'ERROR: Error reading file', trim(file)
          write(*,*) '       fstat = ', fstat
          write(*,*) 'Stopping now, sorry'
          stop
       end if
    end if

    if (1 .ne.  product(shape(x))) then
       if (present(stat)) then
          stat = -1
          deallocate(x)
          return
       else
          write(*,*) 'ERROR: size mismatch while reading file', trim(file)
          write(*,*) 'Stopping now, sorry'
          stop
       end if
    end if
    n = int(x(1,1))
    deallocate(x)
  end subroutine readdata_n

  !! read a single double
  subroutine readdata_x(file,x,stat,uselock)
    implicit none
    real(kind=dbl), intent(out) :: x
    character(*), intent(in) :: file
    integer, intent(out), optional :: stat
    logical, intent(in), optional :: uselock

    real(kind=dbl), pointer :: xx(:,:)
    integer :: fstat

    fstat = 0
    call loaddata(file,xx,fstat,uselock=uselock)
    if (fstat /= 0) then
       if (present(stat)) then
          stat = fstat
          return
       else         
          write(*,*) 'ERROR: Error reading file', trim(file)
          write(*,*) '       fstat = ', fstat
          write(*,*) 'Stopping now, sorry'
          stop
       end if
    end if

    if (1 .ne.  product(shape(xx))) then
       if (present(stat)) then
          stat = -1
          deallocate(xx)
          return
       else
          write(*,*) 'ERROR: size mismatch while reading file', trim(file)
          write(*,*) 'Stopping now, sorry'
          stop
       end if
    end if
    x = xx(1,1)
    deallocate(xx)
  end subroutine readdata_x
!!!
!!! create empty file
!!!
  subroutine create_empty_file(file,status)
    character(*), intent(in) :: file
    integer, intent(out), optional :: status
    integer :: fstat
    open(unit=10, file=file, status='replace', iostat=fstat)
    if (present(status)) status = fstat
    if (fstat /= 0 .and. (present(status) .eqv. .false.) ) then
       write(*,*) 'ERROR: in create_empty_file'
       write(*,*) 'ERROR: could no create file ', trim(file) 
       write(*,*) 'ERROR: iostat: ', fstat
       stop
    end if
    close(10)
  end subroutine create_empty_file
!!!
!!! delete file if it exists
!!!
  subroutine delete_file(file,status)
    character(*), intent(in) :: file
    integer, intent(out), optional :: status
    integer :: fstat
    open(unit=10, file=file, status='old', iostat=fstat)
    if (present(status)) status = fstat
    if (fstat == 0) then
       close(10, status='delete', iostat=fstat)
       if (present(status)) status = fstat
       if (fstat /= 0 .and. (present(status) .eqv. .false.) ) then
          write(*,*) 'ERROR: in delete_file'
          write(*,*) 'ERROR: could no delete file ', trim(file) 
          stop
       end if
    end if
  end subroutine delete_file
!!!
!!! open file with locking
!!!
  subroutine open_with_lock(file,unit,status,timeout,iostat,uselock)
    implicit none
    character(*), intent(in) :: file
    integer, intent(in)  :: unit
    character(*), intent(in), optional :: status
    real, intent(in), optional :: timeout
    integer, intent(out), optional :: iostat
    logical, intent(in), optional :: uselock

    character(len=len_trim(file)+5) :: lockfile
    character(len=9) :: status2
    integer :: iostat2
    real :: timeout2
    integer :: time1, time2, tunit
    logical :: exist, wait, uselock2

    iostat2 = 0

    if (present(uselock)) then
       uselock2 = uselock
    else
       uselock2 = .true.
    end if    

    lockfile = trim(file)//'.lock'

    if (present(timeout)) then
       timeout2 = timeout
    else
       timeout2 = 10.0
    end if

    if (present(status)) then
       status2 = status
    else
       status2 = 'unknown'
    end if

    
    !! wait for the lock file to disappear
    if (uselock2) then
       wait = .true.
       exist = .true.
       call system_clock(count_rate=tunit)
       call system_clock(count=time1)
       do
          inquire(file=lockfile,exist=exist)
          call system_clock(count=time2)
          if (real(time2-time1)/real(tunit) > timeout2) wait=.false.
          if (exist .and. wait) then
             write(*,*) 'waiting for lock file ...'// &
                  'time left:',real(timeout2)-real(time2-time1)/real(tunit)
#if defined(AIX)
!             call SYSTEM('sleep 1')
             call sleep_(1)
#elif defined(RS6K)
             call sleep(1)
#else
             call sleep(1)
#endif
          else
             exit
          end if
       end do

       if (exist) then
          if (present(iostat)) then
             iostat = -1
             return
          else
             stop 'timeout waiting for lockfile'
          end if
       end if
    end if

    !! create lock file
    if (uselock2) then
       call create_empty_file(lockfile,iostat2)
       if (iostat2 /= 0) then
          if (present(iostat)) then
             iostat = iostat2
             close(unit)
             return
          else
             write(*,*) 'Could not create lockfile ',trim(lockfile) 
             stop
          end if
       end if
    end if

    !! open the file
    open(file=file,unit=unit,iostat=iostat2,status=status2)
    if (iostat2 /= 0) then
       if (present(iostat)) then
          iostat = iostat2
          if (uselock2) call delete_file(lockfile,iostat2)
          return
       else
          write(*,*) 'Could not open file ',trim(file)
          write(*,*) 'stopping'
          if (uselock2) call delete_file(lockfile)
          stop
       end if
    end if

    if (present(iostat)) iostat = iostat2
    
  end subroutine open_with_lock
!!!
!!! close file with locking
!!!
  subroutine close_with_lock(unit,file,iostat,uselock)
    implicit none
    integer, intent(in)  :: unit
    character(*), intent(in) :: file
    integer, intent(out), optional :: iostat
    logical, intent(in), optional :: uselock

    character(len=len_trim(file)+5) :: lockfile
    integer :: iostat2
    logical :: uselock2

    uselock2 = .true.
    if (present(uselock)) uselock2 = uselock

    iostat2 = 0

    close(unit=unit,iostat=iostat2)

    if (uselock2) then
       lockfile = trim(file)//'.lock'
       call delete_file(lockfile,iostat2)
    end if

    if (present(iostat)) iostat = iostat2

  end subroutine close_with_lock

!!!
!!! text waitbar
!!!
  subroutine waitbar(pctot,msg)
    implicit none
    real, intent(in) :: pctot
    character(len=*), intent(in), optional :: msg
    !! local variables
    character(len=1), parameter :: bar='=', back=char(8), sp=' '
    real :: pctotloc
    integer :: k, loc
    integer, parameter :: numbars = 60, msglen=128
    character(len=msglen) :: msgloc
    integer :: linelen = numbars + msglen + 9

    if (present(msg)) then
       msgloc = msg(1:min(len_trim(msg),msglen))
    else
       msgloc = ''
    end if

    pctotloc = pctot
    if (pctotloc<0.0) pctotloc=0.0
    if (pctotloc>1.0) pctotloc=1.0
    loc = int(pctotloc*real(numbars))

    ! uses ideas from comp.lang.fortran
    ! first delete to the beginning of line
    write(6,'(256a1)', advance='no') (back, k=1,linelen)
    ! print the percentage and the bar
    write(6,'(2x,1i3,1a1,2x,1a1,256a1)', advance='no') &
         int(pctotloc*100.0),'%','|', (bar, k=1,loc)
    ! print spaces at the end of the waitbar
    write(6,'(256a1,1a1)', advance='no') &
         (sp, k=loc+1,numbars),'|'
    ! add message at the end of line
    write(6,'(1x,a)', advance='no') trim(msgloc)
    ! flush stdout (is this a standard function?)
    call flush(6)

  end subroutine waitbar
!!!
!!! generate unit matrix
!!! x => genmat_eye(n)
!!!
  function genmat_eye(n) result(x)
    implicit none
    integer, intent(in) :: n
    real(kind=dbl), pointer :: x(:,:)
    !!
    integer :: i
    allocate(x(n,n))
    x = 0.0_dbl
    do i=1,n
       x(i,i) = 1.0_dbl
    end do
  end function genmat_eye
!!!
!!! zero the lower diagonal of symmetric x
!!!
  subroutine zerolower(x)
    implicit none
    real(kind=dbl), intent(inout) :: x(:,:)
    integer i,m,n
    m = size(x,1)
    n = size(x,2)
    do i = 1,min(n,m)-1
       x((i+1):m,i) = 0.0_dbl
    end do
  end subroutine zerolower
!!!
!!! copy the lower diagonal from upper diagonal for symmetric x
!!!
  subroutine copylower(x)
    implicit none
    real(kind=dbl), intent(inout) :: x(:,:)
    integer i,m,n
    m = size(x,1)
    n = size(x,2)
    do i=1,min(n,m)-1
       x((i+1):n,i) = x(i,(i+1):n)
    end do
  end subroutine copylower
!!!
!!! simple data printing
!!!
  subroutine printmat_mat(x,maxlines,maxcols,header)
    implicit none
    real(kind=dbl), intent(in) :: x(:,:)
    integer, optional, intent(in) :: maxlines
    integer, optional, intent(in) :: maxcols
    character(len=*), optional, intent(in) :: header

    integer :: i,n,m
    if (present(maxlines).and.maxlines.gt.0) then
       n = min(size(x,1),maxlines)
    else
       n = size(x,1)
    end if
    if (present(maxcols).and.maxcols.gt.0) then
       m = min(size(x,2),maxcols)
    else
       m = size(x,2)
    end if
    if (present(header)) write(*,*) header
    do i=1,n
       write(*,*) real(x(i,1:m))
    end do
  end subroutine printmat_mat
!!!
  subroutine printmat_vec(x,maxlines,header)
    implicit none
    real(kind=dbl), intent(in) :: x(:)
    integer, optional, intent(in) :: maxlines
    character(len=*), optional, intent(in) :: header

    integer :: n
    if (present(maxlines).and.maxlines.gt.0) then
       n = min(size(x),maxlines)
    else
       n = size(x)
    end if

    if (present(header)) write(*,*) header
    write(*,*) real(x(1:n))
  end subroutine printmat_vec
!!!
  subroutine printmat_ivec(x,maxlines,header)
    implicit none
    integer, intent(in) :: x(:)
    integer, optional, intent(in) :: maxlines
    character(len=*), optional, intent(in) :: header
    integer :: i,n

    if (present(maxlines).and.maxlines.gt.0) then
       n = min(size(x),maxlines)
    else
       n = size(x)
    end if

    if (present(header)) write(*,*) header
    do i=1,n
       write(*,*) x(i)
    end do
  end subroutine printmat_ivec

!!! size check utility
  subroutine sizecheck_mat(xmat,n1,n2,txt)
    implicit none
    real(kind=dbl), intent(in) :: xmat(:,:)
    integer, intent(in) :: n1, n2
    character(len=*), intent(in), optional :: txt

    if (size(xmat,1) .ne. n1 .or. size(xmat,2) .ne. n2) then
       if (present(txt)) then
          write(*,*) 'Size mismatch:',txt
       else
          write(*,*) 'Size mismatch'
       end if
       write(*,*) 'rows:',n1,' columns:',n2
       write(*,*) ' Will stop now, sorry'
       stop
    end if
  end subroutine sizecheck_mat

end module matutils
