!!! $Id: MCMC_adapt.F90,v 1.12 2012/06/27 10:10:36 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mcmc library
!!! File: MCMC_adapt.F90
!!! Purpose: proposal covariance matrix adaptation
!!!
!!! Marko Laine 2008 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!!
!!!
subroutine MCMC_adapt(simuind)
  implicit none
  integer :: simuind ! where we are
  integer, save :: stayed2 = 0
  !! local variables
  real(kind=dbl) :: staypc
  integer :: info
  integer, save :: istart=1, istartind=1, lastind=1, lastfreq=0, newfreq
  integer :: histsum

  if (mod(simuind,printint) == 0) then
     if (dodr) then
        write(*,'(A,I10,A,F5.1,A,F5.1,A,F5.1,A)') ' simu i =', simuind, &
             ', stayed % = ',real(stayed)/real(simuind)*100.0, &
             ' /',(1.0-real(draccepted)/real(drtries))*100.0, &
             '  [',real(bndstayed)/real(simuind)*100.0,']'
     else if (erstayed>0) then
        write(*,'(A,I10,A,F5.1,A,F5.1,A)') ' simu i =', simuind, &
             ', stayed % = ',real(stayed)/real(simuind)*100.0, &
             '  [',real(erstayed)/real(simuind)*100.0,'](er)'
     else
        write(*,'(A,I10,A,F5.1,A,F5.1,A)') ' simu i =', simuind, &
             ', stayed % = ',real(stayed)/real(simuind)*100.0, &
             '  [',real(bndstayed)/real(simuind)*100.0,']'
     end if
  end if
  if (simuind == adaptend .and. doadapt > 0) then
     write(*,*) 'note: last adapt'
  end if

  if (doadapt == 0 .and. doburnin == 0) return ! we do not adapt
  if (adaptend > 0 .and. simuind > adaptend) return ! not any more
  !   if (simuind == nsimu) return ! no need to adapt any more
  if (mod(simuind,adaptint) /= 0 .and. &
       mod(simuind,badaptint) /= 0 ) return
  !! but do we want to still calculate covmat?


  !    write(*,*) 'size of chain=',size(chain,1),size(chain,2)
  !    write(*,*) 'size of sschain=',size(sschain,1),size(sschain,2)
  !    write(*,*) 'size of s2chain=',size(s2chain,1),size(s2chain,2)
  !    write(*,*) 'npar =',npar,'nycol=',nycol

!!! do adapt
  if (verbosity>1) write(*,*) 'note: ',simuind,': adapting'

!!! not very many points accepted, scale down    
!!! only during burnin
  if (simuind < burnintime .and. doburnin /= 0 .and. &
       mod(simuind,badaptint) == 0) then
     staypc = dble(stayed)/dble(simuind)
     istartind = chainind ! save current index to chain =====XXXX
     if (staypc>1.0_dbl-scalelimit) then
        R = R/scalefactor
        stayed2 = stayed2 + stayed ! this is not used?
        if (dodr) then
           R2 = R2/scalefactor
           iC = iC*scalefactor*scalefactor
        end if
        if (verbosity>0) write(*,*) 'note: scaling down'
        return
     else if (staypc<scalelimit) then
!!! scale up
        R = R*scalefactor
        stayed2 = stayed2 + stayed
        if (dodr) then
           R2 = R2*scalefactor
           iC = iC/scalefactor/scalefactor
        end if
        if (verbosity>0) write(*,*) 'note: scaling up'
        return
     else if (greedy /= 0) then !!! greedy adaptation
        if (verbosity>0) write(*,*) 'note: greedy adaptation'
        !! in burnin allways start with cmat0
        chainwsum = dble(initcmatn)
        chaincmat = cmat0
        chainmean = par0
        !          call covmat(chain(1:chainind,1:npar),chaincmat,(/1.0d0/), &
        !              chainmean, chainwsum, .false.)
        call covmat(chain(1:chainind,1:npar),chaincmat,(/1.0e0_dbl/), &
             chainmean, chainwsum, .true.)

        if (verbosity>9) write(*,*) 'Greedy: called covmat'
        if (verbosity>9) write(*,*) istart, chainind, chainwsum ! debug
        if (verbosity>9) write(*,*) lastfreq, newfreq

        !! save frequency info for update in adapt below
        lastfreq = int(chain(chainind,ncolchain))

     end if
     lastind  = chainind

!!! Burn in is done,  adapt
  else if (simuind >= burnintime+adaptint+adapthist .and. doadapt /= 0) then 

     !! this is the first time we adapt
     if (simuind == burnintime+adaptint+adapthist) then
        if (verbosity>1) write(*,*) 'note: first adapt'
        chainwsum = dble(initcmatn)
        chaincmat = cmat0
        chainmean = par0
!!! how about AP??
     end if

     if (adapthist > 1) then ! AP adaptation at simuind intervals

        !! need to calculate istart
        istart = chainind
        histsum = int(chain(istart,ncolchain))
        do while (histsum < adapthist .and. istart > 1) 
           istart = istart - 1
           histsum = histsum + int(chain(istart,ncolchain))
        end do
        if (istart == chainind) then
           if (verbosity>0) write(*,*) 'Warning: no move in last adaptint steps'
           if (verbosity>0) write(*,*) '         Adaptint = ',adaptint
        end if
        newfreq = int(chain(istart,ncolchain))
        chain(istart,ncolchain) = dble(newfreq-histsum+adapthist)
        call covmat(chain(istart:chainind,1:npar), &
             chaincmat,chain(istart:chainind,ncolchain), &
             chainmean, chainwsum, .false.)
        chain(istart,ncolchain) = dble(newfreq) !! resume freq
        !         write(*,*) 'AP adaptation: ', istart, chainind, histsum
        !         write(*,*) newfreq, newfreq-histsum+adapthist

     else ! adapt from the end of burn in

        !! check frequency data
        newfreq = int(chain(lastind,ncolchain))
        chain(lastind,ncolchain) = dble(newfreq-lastfreq)

        istart = lastind  !!!  +1
        !          if (istart == chainind) then
        !             write(*,*) 'Warning: no move in last adaptint steps'
        !             write(*,*) '         Adaptint = ',adaptint
        !          end if
        call covmat(chain(istart:chainind,1:npar),chaincmat, &
             chain(istart:chainind,ncolchain), chainmean, chainwsum, .true.)
        
        if (verbosity>9) write(*,*) 'D: called covmat'
        if (verbosity>9) write(*,*) 'D:',istart,chainind,real(chainwsum),lastfreq,newfreq

        chain(lastind,ncolchain) = dble(newfreq) !! resume freq
        lastfreq = int(chain(chainind,ncolchain))
        lastind  = chainind

     end if

  else if (doadapt /= 0) then
     if (verbosity>3)  write(*,*) 'Not enough good chain to do real adaptation yet'
     return
  else !!! nothing to do
     return
  end if

  call MCMC_calculate_R(chaincmat,info)
  if (info .ne. 0) then
     write(*,*) 'Warning: error in Chol/SVD, not adapting, info = ', info
  end if

  return
end subroutine MCMC_adapt

!!! 
!!! calculate Cholesky factorization of covmat and store it in R, scale R
!!! also calculate other global variables needed
!!! scaling 2.4/sqrt(d): Gelman, et al. Bayesian Statistics 5,
!!! Efficient Metropolis Jumping Rules, 1996.
subroutine MCMC_calculate_R(cmat,info)
  implicit none
  real(kind=dbl), intent(inout) :: cmat(:,:)
  integer, intent(out) :: info

  real(kind=dbl) :: R0(npar,npar)

  info = 0
  if (doscam) then
     !! scam adaptation = rotation
     call scam_svd(cmat,R0,qcovstd,condmax,info=info)
     if (info == -1) then ! cmat adjusted so we change cmat
        !! SCAM?    cmat = matmul(R0,transpose(R0))
        info = 0
     end if
     !! scale in SCAM?
     !! qcovstd = qcovstd*2.4      ! ???
     R = R0
     if (verbosity>1) write(*,*) 'SCAM adaptation'
     !! no DR code yet

  else
     !! calculate Cholesky factorization of covmat and store it in R, scale R
     if (usesvd /= 0) then
        call covtor_svd(cmat,R0,condmax,info=info)
        if (info == -1) then ! cmat adjusted so we change cmat
           cmat = matmul(R0,transpose(R0))
           info = 0
        end if
     else
        call covtor(cmat,R0,info=info)
     end if
     if (info .ne. 0) then
        write(*,*) 'warning: error in Chol/svd, covmat singular?, info = ', info
     else
        R = R0*2.4_dbl/sqrt(dble(npar)) !!
        if (dodr) then ! invert covmat using factorization in R        
           iC = R           
           call dpotri('u',npar,iC,npar,info)
           if (info .ne. 0) then
              write(*,*) 'ERROR: cannot invert cmat, info = ', info
              stop
           end if
           R2 = R/drscale
        end if
     end if
  end if


end subroutine MCMC_calculate_R

!!!
!!!
!!! update the covariance matrix and chain mean stored in global values
!!! chaincmeat, chainmean, ...
!!!
!!!
subroutine MCMC_updatecov(newpar)
  implicit none
  real(kind=dbl), intent(in), dimension(npar) :: newpar

!  call covmat(newpar,chaincmat,1,chainmean,chainwsum, .true.)
!
end subroutine MCMC_updatecov
