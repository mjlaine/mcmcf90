!!! $Id: priorfun.f90,v 1.4 2012/06/27 10:10:37 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mcmc library
!!! File: priorfun.f90
!!! Purpose:
!!!
!!!  default version of the prior for mcmc code
!!!  this calculates -2*log(p(theta)
!!! 
!!!  you can replace this by providing your own priorfun function
!!!
!!! Marko Laine 2001 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!!

!!! 2007-03-26: new version with Gaussian priors read from a dat file
!!! priorsfile contains two rows as:
!!!
!!!    mu1  mu2  mu3 ...
!!!    sig1 sig2 sig3 ...
!!!
!!! (OBS: sigs are stds)
!!! sig<=0 means sig = infinity
!!! priorsfile is given in the mcmc namelist
!!! the default value is to use no priors = flat/uninformative priors

!!!
!!! This returns -2*log(p(theta))
!!!
function priorfun(theta,len) result(priss)
  use mcmcmod, only:priorsfile
  use mcmcprec, only:dbl,ik4

  implicit none

  real(dbl) priss
  integer(ik4) len
  real(dbl) theta(len)
  integer i, fstat
  logical, save,allocatable :: inds(:)
  integer, save :: plen, np, first_run=1
  real(dbl), save, allocatable, dimension(:,:) :: pmus

  priss = 0.0e0_dbl

  !! if the string priorsfile is empty, no priors set
  if (len_trim(priorsfile) <= 0) then
     if (first_run == 1) then
!        write(*,*) 'note: no priors file'
        first_run = 0
     end if
     priss = 0.0e0_dbl
     return
  end if

  !! during the first call load prior mu and sigma from 'priorsfile'
  if (first_run == 1) then
     allocate(inds(len),pmus(2,len))
     open(333,file=priorsfile,status='old',iostat=fstat) 
     if (fstat /= 0) then
        write(*,*) 'could not open file ',trim(priorsfile)
        write(*,*) 'create it or remove priorsfile from namelist'
        stop
     end if
     read(333,*,iostat=fstat) (pmus(1,i), i=1,len)
     if (fstat /= 0) then
        write(*,*) 'priors.dat should have  2*npar elements'
        write(*,*) 'npar=',len
        stop
     end if
     read(333,*,iostat=fstat) (pmus(2,i), i=1,len)
     if (fstat /= 0) then
        write(*,*) 'priors.dat should have  2*npar elements'
        write(*,*) 'npar=',len
        stop
     end if
     close(333)
     inds = (pmus(2,:) > 0.0e0_dbl)
     np = count(pmus(2,:) > 0.0e0_dbl)
     plen = len
     first_run = 0

     write(*,*) 'read ',trim(priorsfile),': npriors=', np
!     write(*,*) 'mu   :', real(pmus(1,:))
!     write(*,*) 'sigma:', real(pmus(2,:))

  end if

  if ( len /= plen ) then
     write(*,*) 'priors.f90: something wrong with prior length plen'
     write(*,*) 'plen =',plen,' length(theta) =',len
     stop
  end if

  !! Gaussian prior distribution
  priss = 0.0e0_dbl
  if (np > 0) then
     priss = sum(((theta(:)-pmus(1,:))/pmus(2,:))**2, mask=inds)
  end if

  return
end function priorfun
