!!! $Id: MCMC_aux.F90,v 1.19 2012/07/03 10:57:33 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mcmc library
!!! File: MCMC_aux.F90
!!! Purpose: auxiliary MCMC routines
!!!
!!! Marko Laine 2008 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!!
!!! some auxiliary MCMC routines for the mcmc modulde
!!! contains: MCMC_writechains, MCMC_cleanup, MCMC_savechain

!!!
!!! write chains to file
!!!
subroutine MCMC_writechains(status)
  implicit none
  integer, intent(out) :: status
  integer :: fstat

  status = 0
  fstat = 0

  if (chainfile(len_trim(chainfile)-3:len_trim(chainfile)) == '.mat') then
     call writemat(chainfile, chain(1:chainind,:),'chain',fstat)
  else
     call writedata(chainfile, chain(1:chainind,:))
  end if
  if (ssfile(len_trim(ssfile)-3:len_trim(ssfile)) == '.mat') then
     call writemat(ssfile, sschain(1:chainind,:),'sschain',fstat)
  else
     call writedata(ssfile, sschain(1:chainind,:))
  end if
  if (updatesigma /= 0) then
     if (s2file(len_trim(s2file)-3:len_trim(s2file)) == '.mat') then
        call writemat(s2file,s2chain(1:simuind,:),'s2chain',fstat)
     else
        call writedata(s2file,s2chain(1:simuind,:))
     end if
  end if
  if (usrfunlen > 0) then
     !      call writedata('usrfunmat.dat',usrfunmat(1:chainind,:))
     call writemat('usrfunmat.mat',usrfunmat(1:chainind,:),'usrfunmat',fstat)
  end if
  !! save the final covariance matrix.
  call writedata(covffile,chaincmat,uselock=.true.)
  !! save also mcmccovn.dat if it was used
  if (len_trim(covnfile) > 0) then
     initcmatn = int(chainwsum) !! check this
     call writedata(covnfile,reshape((/dble(initcmatn)/),(/1,1/)),uselock=.true.)
     if (verbosity>2) write(*,*) 'note: wrote initcmatn:', initcmatn
  end if

  call writedata(meanfile, &
       reshape((/chainmean/),(/size(chainmean,1),1/)) &
       ,uselock=.true.)
  call writedata(parffile,chain(chainind,1:npar))
  if (updatesigma /= 0) then
     call writedata(sigma2ffile, &
          transpose(reshape( (/s2chain(simuind,1:nycol), dble(nobs)/), &
          (/nycol,2/) )))
  end if

  status = fstat

  if (status /= 0) then
     write(*,*) 'Something wrong when writing files'
  else
     if (updatesigma.ne.0) then
        if (verbosity>0) write(*,*) 'note: saved results in ', trim(chainfile), ' and ', &
             trim(s2file), ' and ', trim(ssfile),'.'
     else
        if (verbosity>0) write(*,*) 'note: saved results in ', trim(chainfile), ' and ', &
             trim(ssfile),'.'
     end if
  end if

  !! lastly write mcmcinit.nml
  initcmatn = initcmatn + simuind
  burnintime = 0
  if (len_trim(nmlffile)>0) & 
       call write_mcmcinit_namelist(fstat,file=nmlffile)

end subroutine MCMC_writechains
!!!
!!! clear memory and other cleanup stuff
!!! This is still not complete, so only one MCMC run per program invocation 
!!! is actually supported. 
subroutine MCMC_cleanup(status)
  implicit none
  integer, intent(out) :: status

  status = 0
  deallocate(par0, cmat0, sigma2, nobs)
  deallocate(chain, sschain, R, chainmean, chaincmat)
  if (updatesigma /= 0) then
     deallocate(s2chain)
  end if
  if (dodr) then
     deallocate(R2,iC)
  end if
  if (doscam) then
     deallocate(qcovstd)
  end if

  nparok = .false.
  par0ok = .false.
  cmat0ok = .false.
  sigma2ok = .false.

  npar = 0
  nsimu = 0
  inited = 0
  call random_eoj()           ! random number end of job
  call MCMC_dump_end()        ! close possible dump files

end subroutine MCMC_cleanup

!!!
!!! Save the current chain values
!!! 
subroutine MCMC_savechain(par,ss,sspri,sigma2,reject)

  implicit none
  real(kind=dbl), intent(in) :: par(npar),ss(nycol),sspri,sigma2(nycol)
  logical, intent(in) :: reject

! we need to initialize the writing somewhere
!  character(len=100), intent(in) :: task

  character(len=100) :: save_method = 'memory'

  !! chainind is global, lchainind local
  integer, save :: lchainind = 1, nrep = 1

  !! save methods: memory, disk, interval

  !! disk save: binary dump with mat v4 header, saved row wise (transposed). 
  !! Either flush (close) the file after each write or only at the end.
  if (save_method == 'disk') then

     if (reject) then
        nrep = nrep+1
     else
        nrep = 1
     end if

     chain(chainind,:) = par
     sschain(chainind,:) = ss
     if (updatesigma /= 0) s2chain(chainind,:) = sigma2

     !! if we have filled up the chain then append it to a binary matlab v4 file
     if (chainind >= size(chain,1)) then
        call addtomat(chainfile,transpose(chain))
        call addtomat(ssfile,transpose(sschain))
        if (updatesigma /= 0) call addtomat(s2file,transpose(s2chain))
        chainind = 1
     else
        chainind = chainind + 1
     end if

     return
  end if

  !! memory save: allocate space for the whole chain at once
  if (save_method == 'memory') then

     if (reject) then
        chain(chainind,ncolchain) = chain(chainind,ncolchain) + 1.0_dbl
     else
        chainind = chainind + 1
        chain(chainind,1:npar)    = par
        chain(chainind,ncolchain) = 1.0_dbl
        sschain(chainind,1:nycol) = ss
     end if
     sschain(chainind,nycol+1) = chain(chainind,ncolchain) ! copy the count

     !! save the sigma2 chain
     if (updatesigma /= 0) then
        s2chain(simuind,1:nycol) = sigma2
     end if

     return
  end if

  !! We should not get here
  write(*,*) 'ERROR: unknown chain save method'
  stop
  
end subroutine MCMC_savechain
