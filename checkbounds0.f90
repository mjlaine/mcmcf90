!!! default/dummy checkbunds for libmcmcrun

function checkbounds(theta)
  implicit none
  real*8 theta(:)
  logical checkbounds

  logical, save :: warn = .true.

  if (warn) then
     write(*,*) 'note: using the default dummy checkbound routine'
     warn = .false.
  end if
  
  checkbounds = .true.
  return
  
end function checkbounds

