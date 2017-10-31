libmcmcrun.a - MCMC library

inputfiles for libmcmrun

mcmcpar.dat  : initial first guess for theta
mcmccov.dat  : initial proposal covariance, a first guess for covariance of theta
mcmcsigma2.dat : initial error variance and number of observations

MCMC run time parameters are read from mcmcinit.nml:

!
! mcmcinit.nml -- initialization for mcmc F90 code
!
&mcmc
 nsimu       = 10000   ! length of the chain
 updatesigma = 0       ! update error variance
 N0          = 1       ! prior for error variance,
 S02         = 0       !    1/s^2 ~ Gamma(N0/2,2/N/S02), or fixed sigma2
 priorsfile  = ''      ! file from which prior parameters are read
 chainfile   = 'chain.mat'    ! file to save the chain
 s2file      = 's2chain.mat'  ! file to save sigma2 chain
 ssfile      = 'sschain.mat'  ! file to save sum-of-squares chain

 doburnin    = 0       ! do we have 'burn-in'
 burnintime  = 1000    ! burn-in time
 scalelimit  = 0.05    ! when to scale
 scalefactor = 2.5     ! scale factor
 drscale     = 0.0     ! DR scale
 greedy      = 1       ! "greedy" burn in adaptation 
 badaptint   = -1      ! burn-in adaptation interval, if < 0 use adaptint
 doadapt     = 1       ! do we adapt
 adaptint    = 100     ! interval for adaptation
 adaptend    = 0       ! end adaptation at this time (if >0)
 initcmatn   = 0       ! "imaginary chain size" for initial proposal cmat 
 printint    = 1000    ! interval to print statistics
/


NOTES:

chainfile can be file.mat or file.dat, in the first case a binary
Matlab mat-file is written, in the second case an ascii file. The last
column of the chain file tells the number of times the corresponding
line repeats itself. So the actual chain (if needed) must be expanded
from this information (e.g. by loadchain matlab function).

