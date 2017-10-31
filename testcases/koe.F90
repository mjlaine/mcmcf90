!!! random tests

program koe

  use matutils
  use mcmcprec
  use matfiles
  implicit none


  real(kind=dbl), allocatable, target :: x(:,:)
  real(kind=dbl), pointer :: p(:,:)

  call readmat('data.mat',p)

  write(*,*) shape(p)

  call printmat(p, header='data')

  return

  call loaddata('data.dat',p)

  write(*,*) associated(p)
  write(*,*) allocated(x)
  allocate(x(size(p,1),size(p,2)))
  x = p
  call printmat(x)
  deallocate(p)

  write(*,*) associated(p)
  write(*,*) allocated(x)

  p => x

  deallocate(x)
  allocate(x(size(p,1),size(p,2)))
  x=0

!  deallocate(p)

  write(*,*) associated(p)
  write(*,*) allocated(x)

  call printmat(p)

end program koe
