subroutine bounds( kL, lambda, nl, sigma_L0MP, beta, ebar1, elower, elowernew, euppernew, &
                   emiddle, qboundMP, rho, qfixed, test, eref, fmprec )
  use fmzm
  implicit none
  
  ! This subroutine computes IMPROVED lower bounds using arbirtary precision arithmetic
  ! Please refer to the "Notes for Eli", and section 8 in particular
  ! On entry
  !        nl        is the largest size of the Lanczos matrices (for dimensioning only)
  !        kL        is the "reference" Lanczos dimension (i.e. it is 'L+1')
  !        lambda    is the matrix of eigenvalues for dimensions 1,2,..k
  !        sigma_L0  is the ground-state variance of the reference problem 
  !        beta      is the first diagonal of the Lanczos matrix ( beta(1), beta(2),etc. for beta_1, beta_2,etc.)
  !        elower    is a lower bound estimate
  !        ebar1     is the lower bound Ebar_1
  !        qboundMP  is an estimate of the Q-ratio
  !        qfixed    is a flag to use an energy-independent Q-bound (qfixed = T)
  !        test      is used only when qfixed = T, and instruct the code to use the provided value of qboundMP
  ! On exit
  !        elowernew is the improved estimate for the lower bound
  !        euppernew is the improved estimate for the upper bound
  !        emiddle   is the improved mean energy
  !        qboundMP  is an upper bound to the Q-ratio effectively used
  !        rho       is (half the squared) gap between the upper and lower energy estimates

  logical(1) :: test, qfixed, quasiexact
  integer(4) :: fmprec
  integer(4) :: kL, L, nl
  real(8) :: lambda(nl,kL), beta(nl-1) 
  real(8) :: sigma_L0, ebar1, elower, elowernew, euppernew, qbound, rho, eref, emiddle
  real(8) :: zero, one, two
  integer(4) :: i, j, n
  
  TYPE(FM) :: lambdaMP(nl,kL), betaMP(nl-1)
  TYPE(FM) :: gammaMP(kL-1), pMP(kL-1), plessMP
  TYPE(FM) :: sigma_L0MP, ebar1MP, elowerMP, elowernewMP, eupperMP, qboundMP
  TYPE(FM) :: zeroMP, oneMP, twoMP

  TYPE(FM) :: lamrefMP, erefMP
  TYPE(FM) :: ebestMP, deltaMP
  TYPE(FM) :: anMP, bnMP, cnMP, banMP, canMP
  TYPE(FM) :: tmpMP, rhoMP
  
  L = kL-1                     ! to simplify reading, L is the Lanczos order of the reference problem
  quasiexact = .false.         ! instruct the code to evaluate the L-polynomial at the exact energy eref
                               ! and, if qfixed=.true., to use the exact energy for computing the energy-independent Q0
  ! Set precision
  call FMSET(fmprec)
  
  ! Convert from double precision to multiple precision
  zero = 0.d0
  one = 1.d0
  two = 2.d0
  CALL FMDP2M(zero, zeroMP)
  CALL FMDP2M(one, oneMP)
  CALL FMDP2M(two, twoMP)
  do i = 1, kL
     do j = 1, i
        CALL FMDP2M( lambda(j,i),lambdaMP(j,i) )
     enddo
  enddo
  lamrefMP = lambdaMP(1,kL)      ! this is the approximate ground-state energy from the Lth order problem 
  do i = 1, nl-1
     CALL FMDP2M( beta(i),betaMP(i) )
  enddo
  call FMDP2M( ebar1, ebar1MP)
  call FMDP2M( elower, elowerMP)
  call FMDP2M( eref, erefMP)

  ! compute the relevant coefficients, for n = 1,2, .., L
  gammaMP(1) = betaMP(2)**2
  do n = 2, L
     gammaMP(n) = betaMP(n+1)**2 * gammaMP(n-1)
  enddo
  
  do n = 1, L
     pMP(n) = oneMP / gammaMP(n)
     do j = 1, n
        pMP(n) = ( lamrefMP-lambdaMP(j+1,n+1) )**2 * pMP(n)    ! this is the squared \tilde(p)_n evaluated at lambdaref, i.e.
                                                               ! the product of all eigenvalues but the ground-state, normalized by gamma_n
     enddo
  enddo
  
  plessMP = oneMP / gammaMP(L)

  ebestMP = (elowerMP + lamrefMP)/twoMP
  deltaMP = (lamrefMP-elowerMP)/twoMP

  if(quasiexact) then
     do j = 1, L
        plessMP = (erefMP - lambdaMP(j+1,L+1) )**2 * plessMP      ! this is the L polynomial evaluated at the exact energy
     enddo
  else
     do j = 1, L
        plessMP = (elowerMP - lambdaMP(j+1,L+1) )**2 * plessMP    ! this is the L polynomial evaluated at the lower bound energy
        ! test alternative bounds
        ! plessMP = (elowerMP+deltaMP-lambdaMP(j+1,L+1) )*(elowerMP-deltaMP-lambdaMP(j+1,L+1) )*plessMP
        ! 
     enddo
  endif
  
  ! qboundMP = plessMP / betaMP(1)**2*sigma_L0MP**2 - oneMP      ! this is the used Q-bound at the eigenvalue energy
  ! call FMM2DP(qboundMP, qbound)
  
  ! build the coefficients of the quadratic form:
  ! this only works for L+1 < nlmax, because we need the variance associated to lambda_0^L,
  ! i.e. the coupling with the L+1 Lanczos state
  
  anMP = oneMP
  bnMP = lambdaMP(1,1)      ! this is alpha_0
  cnMP = bnMP**2            ! this is alpha_0^2

  do n = 1, L
     anMP = anMP + pMP(n)
     bnMP = bnMP + pMP(n) * lambdaMP(1,n+1)
     cnMP = cnMP + pMP(n) * (lambdaMP(1,n+1))**2
  enddo

  if (qfixed) then
     if( .not.test ) then
        if(quasiexact) then
           qboundMP = - oneMP + plessMP * sigma_L0MP**2/betaMP(1)**2 * (ebar1MP - erefMP)/(ebar1MP -lamrefMP)
        else
           qboundMP = - oneMP + plessMP * sigma_L0MP**2/betaMP(1)**2 * (ebar1MP - elowerMP)/(ebar1MP -lamrefMP)
        endif
     endif
     cnMP = cnMP - betaMP(1)**2 * qboundMP  ! use energy-independent Q
  else
     bnMP = bnMP - oneMP/twoMP * sigma_L0MP**2 * plessMP / (ebar1MP-lamrefMP)
     cnMP = cnMP + betaMP(1)**2 - ebar1MP* sigma_L0MP**2 * plessMP / (ebar1MP-lamrefMP)
  endif

  ! ratio
  banMP = bnMP / anMP
  canMP = cnMP / anMP
  
  rhoMP = banMP**2 - canMP                   ! this should be larger than zero
  call FMM2DP(rhoMP, rho)                    ! output rho before taking the square root

  call FMMAX(zeroMP, rhoMP, tmpMP )
  call FMEQ(tmpMP, rhoMP)
  call FMSQRT_R1(rhoMP)

  call FMM2DP( banMP, emiddle )
  elowerMP = banMP -rhoMP
  eupperMP = banMP +rhoMP

  ! convert final result in double precision
  call FMM2DP(elowerMP, elowernew)
  call FMM2DP(eupperMP, euppernew)
  
  
end subroutine bounds


subroutine tbounds( kL, lambda, nl, sigma_L0MP, beta, ebar1, elower, elowernew, euppernew, &
                    qboundMP, rho, qfixed, fmprec )
  use fmzm
  implicit none
  
  ! This subroutine computes TIGHT lower bounds using arbirtary precision arithmetic
  ! assuming that the Ritz error is less than that of the input elower
  
  ! Please refer to the "Notes for Eli", and section 8 in particular
  ! On entry
  !        nl        is the largest size of the Lanczos matrices (for dimensioning only)
  !        kL        is the "reference" Lanczos dimension (i.e. it is 'L+1')
  !        lambda    is the matrix of eigenvalues for dimensions 1,2,..k
  !        sigma_L0  is the ground-state variance of the reference problem 
  !        beta      is the first diagonal of the Lanczos matrix ( beta(1), beta(2),etc. for beta_1, beta_2,etc.)
  !        elower    is a lower bound estimate
  !        ebar1     is the lower bound Ebar_1
  !        qfixed    is a flag to use an energy-independent Q-bound (qfixed = T)
  ! On exit
  !        elowernew is the improved estimate for the lower bound
  !        euppernew is the improved estimate for the upper bound
  !        qboundMP  is an upper bound to the Q-ratio effectively used
  !        rho       is (half the squared) gap between the upper and lower energy estimates

  logical(1) :: qfixed
  integer(4) :: fmprec
  integer(4) :: kL, L, nl
  real(8) :: lambda(nl,kL), beta(nl-1) 
  real(8) :: sigma_L0, ebar1, elower, elowernew, euppernew, qbound, rho
  real(8) :: zero, one, two
  integer(4) :: i, j, n
  
  TYPE(FM) :: lambdaMP(nl,kL), betaMP(nl-1)
  TYPE(FM) :: gammaMP(kL-1), pMP(kL-1), plessMP
  TYPE(FM) :: sigma_L0MP, ebar1MP, elowerMP, elowernewMP, eupperMP, qboundMP
  TYPE(FM) :: zeroMP, oneMP, twoMP

  TYPE(FM) :: lamrefMP, erefMP
  TYPE(FM) :: ebestMP, deltaMP
  TYPE(FM) :: anMP, bnMP, cnMP, banMP, canMP
  TYPE(FM) :: tmpMP, rhoMP
  
  L = kL-1                     ! to simplify reading, L is the Lanczos order of the reference problem

  ! Set precision
  call FMSET(fmprec)
  
  ! Convert from double precision to multiple precision
  zero = 0.d0
  one = 1.d0
  two = 2.d0
  CALL FMDP2M(zero, zeroMP)
  CALL FMDP2M(one, oneMP)
  CALL FMDP2M(two, twoMP)
  do i = 1, kL
     do j = 1, i
        CALL FMDP2M( lambda(j,i),lambdaMP(j,i) )
     enddo
  enddo
  lamrefMP = lambdaMP(1,kL)      ! this is the approximate ground-state energy from the Lth order problem 
  do i = 1, nl-1
     CALL FMDP2M( beta(i),betaMP(i) )
  enddo
  call FMDP2M( ebar1, ebar1MP)
  call FMDP2M( elower, elowerMP)

  erefMP = (lamrefMP + elowerMP ) /twoMP
  
  ! compute the relevant coefficients, for n = 1,2, .., L
  gammaMP(1) = betaMP(2)**2
  do n = 2, L
     gammaMP(n) = betaMP(n+1)**2 * gammaMP(n-1)
  enddo
  
  do n = 1, L
     pMP(n) = oneMP / gammaMP(n)
     do j = 1, n
        pMP(n) = ( lamrefMP-lambdaMP(j+1,n+1) )**2 * pMP(n)    ! this is the squared \tilde(p)_n evaluated at lambdaref, i.e.
                                                               ! the product of all eigenvalues but the ground-state, normalized by gamma_n
     enddo
  enddo
  
  plessMP = oneMP / gammaMP(L)

  do j = 1, L
     ! this is the tight bound that one can use when the Ritz estimate is better than elower
     plessMP = (erefMP-lambdaMP(j+1,L+1) )**2 *plessMP
  enddo
  
  ! build the coefficients of the quadratic form:
  ! this only works for L+1 < nlmax, because we need the variance associated to lambda_0^L,
  ! i.e. the coupling with the L+1 Lanczos state
  
  anMP = oneMP
  bnMP = lambdaMP(1,1)      ! this is alpha_0
  cnMP = bnMP**2            ! this is alpha_0^2

  do n = 1, L
     anMP = anMP + pMP(n)
     bnMP = bnMP + pMP(n) * lambdaMP(1,n+1)
     cnMP = cnMP + pMP(n) * (lambdaMP(1,n+1))**2
  enddo

  if (qfixed) then
     qboundMP = - oneMP + plessMP * sigma_L0MP**2/betaMP(1)**2 * (ebar1MP - elowerMP)/(ebar1MP -lamrefMP)
     cnMP = cnMP - betaMP(1)**2 * qboundMP  ! use energy-independent Q
  else
     bnMP = bnMP - oneMP/twoMP * sigma_L0MP**2 * plessMP / (ebar1MP-lamrefMP)
     cnMP = cnMP + betaMP(1)**2 - ebar1MP* sigma_L0MP**2 * plessMP / (ebar1MP-lamrefMP)
  endif

  ! ratio
  banMP = bnMP / anMP
  canMP = cnMP / anMP
  
  rhoMP = banMP**2 - canMP                   ! this should be larger than zero
  call FMM2DP(rhoMP, rho)                    ! output rho before taking the square root

  call FMMAX(zeroMP, rhoMP, tmpMP )
  call FMEQ(tmpMP, rhoMP)
  call FMSQRT_R1(rhoMP)

  elowerMP = banMP -rhoMP
  eupperMP = banMP +rhoMP

  ! convert final result in double precision
  call FMM2DP(elowerMP, elowernew)
  call FMM2DP(eupperMP, euppernew)
  
  
end subroutine tbounds


subroutine qbounds( kL, lambda, nl, sigma_L0MP, beta, ebar1, elower, elowernew, euppernew, &
                   emiddle, qboundMP, rho, qfixed, test, eref, fmprec )
  use fmzm
  implicit none
  
  ! This subroutine computes IMPROVED lower bounds using arbirtary precision arithmetic
  ! On entry
  !        nl        is the largest size of the Lanczos matrices (for dimensioning only)
  !        kL        is the "reference" Lanczos dimension (i.e. it is 'L+1')
  !        lambda    is the matrix of eigenvalues for dimensions 1,2,..k
  !        sigma_L0  is the ground-state variance of the reference problem 
  !        beta      is the first diagonal of the Lanczos matrix ( beta(1), beta(2),etc. for beta_1, beta_2,etc.)
  !        elower    is a lower bound estimate
  !        ebar1     is the lower bound Ebar_1
  !        qboundMP  is an estimate of the Q-ratio
  !        qfixed    is a flag to use an energy-independent Q-bound (qfixed = T)
  !        test      is used only when qfixed = T, and instruct the code to use the provided value of qboundMP
  ! On exit
  !        elowernew is the improved estimate for the lower bound
  !        euppernew is the improved estimate for the upper bound
  !        emiddle   is the improved mean energy
  !        qboundMP  is an upper bound to the Q-ratio effectively used
  !        rho       is (half the squared) gap between the upper and lower energy estimates

  logical(1) :: test, qfixed, quasiexact
  integer(4) :: fmprec
  integer(4) :: kL, L, nl
  real(8) :: lambda(nl,kL), beta(nl-1) 
  real(8) :: sigma_L0, ebar1, elower, elowernew, euppernew, qbound, rho, eref, emiddle
  real(8) :: zero, one, two
  integer(4) :: i, j, n
  
  TYPE(FM) :: lambdaMP(nl,kL), betaMP(nl-1)
  TYPE(FM) :: gammaMP(kL-1), pMP(kL-1), plessMP
  TYPE(FM) :: sigma_L0MP, ebar1MP, elowerMP, elowernewMP, eupperMP, qboundMP
  TYPE(FM) :: zeroMP, oneMP, twoMP

  TYPE(FM) :: lamrefMP, erefMP, alpha0MP
  TYPE(FM) :: ebestMP, deltaMP
  TYPE(FM) :: anMP, bnMP, cnMP, banMP, canMP
  TYPE(FM) :: tmpMP, rhoMP
  
  L = kL-1                     ! to simplify reading, L is the Lanczos order of the reference problem

  ! Set precision
  call FMSET(fmprec)
  
  ! Convert from double precision to multiple precision
  zero = 0.d0
  one = 1.d0
  two = 2.d0
  CALL FMDP2M(zero, zeroMP)
  CALL FMDP2M(one, oneMP)
  CALL FMDP2M(two, twoMP)
  do i = 1, kL
     do j = 1, i
        CALL FMDP2M( lambda(j,i),lambdaMP(j,i) )
     enddo
  enddo

  lamrefMP = lambdaMP(1,kL)      ! this is the approximate ground-state energy from the Lth order problem 
  alpha0MP = lambdaMP(1,1)       ! this is alpha_0

  do i = 1, nl-1
     CALL FMDP2M( beta(i),betaMP(i) )
  enddo
  call FMDP2M( ebar1, ebar1MP)
  call FMDP2M( elower, elowerMP)
  call FMDP2M( eref, erefMP)

  
  ! compute the relevant coefficients, for n = 1,2, .., L
  gammaMP(1) = betaMP(2)**2
  do n = 2, L
     gammaMP(n) = betaMP(n+1)**2 * gammaMP(n-1)
  enddo
  
  do n = 1, L
     pMP(n) = oneMP / gammaMP(n)
     do j = 1, n
        pMP(n) = ( lamrefMP-lambdaMP(j+1,n+1) )**2 * pMP(n)    ! this is the squared \tilde(p)_n evaluated at lambdaref, i.e.
                                                               ! the product of all eigenvalues but the ground-state, normalized by gamma_n
     enddo
  enddo
  
  ! build the coefficients of the quadratic form:
  ! this only works for L+1 < nlmax, because we need the variance associated to lambda_0^L,
  ! i.e. the coupling with the L+1 Lanczos state
  
  anMP = oneMP
  bnMP = lambdaMP(1,1)      ! this is alpha_0
  cnMP = bnMP**2            ! this is alpha_0^2

  do n = 1, L
     anMP = anMP + pMP(n)
     bnMP = bnMP + pMP(n) * lambdaMP(1,n+1)
     cnMP = cnMP + pMP(n) * (lambdaMP(1,n+1))**2
  enddo

  if (qfixed) then
     qboundMP = (alpha0MP - elowerMP)/(ebar1MP -alpha0MP)
     cnMP = cnMP - betaMP(1)**2 * qboundMP  ! use energy-independent Q
  else
     bnMP = bnMP - oneMP/twoMP * betaMP(1)**2  / (ebar1MP-alpha0MP)
     cnMP = cnMP - betaMP(1)**2 * alpha0MP / (ebar1MP-alpha0MP)
  endif

  ! ratio
  banMP = bnMP / anMP
  canMP = cnMP / anMP
  
  rhoMP = banMP**2 - canMP                   ! this should be larger than zero
  call FMM2DP(rhoMP, rho)                    ! output rho before taking the square root

  call FMMAX(zeroMP, rhoMP, tmpMP )
  call FMEQ(tmpMP, rhoMP)
  call FMSQRT_R1(rhoMP)

  call FMM2DP( banMP, emiddle )
  elowerMP = banMP -rhoMP
  eupperMP = banMP +rhoMP

  ! convert final result in double precision
  call FMM2DP(elowerMP, elowernew)
  call FMM2DP(eupperMP, euppernew)
  
  
end subroutine qbounds

subroutine bounds_new( kL, lambda, nl, sigma_L0MP, beta, ebar1, elower, elowernew, &
                   qboundMP, test, qflag, eref, fmprec, etemple )
  use fmzm
  implicit none
  
  ! This subroutine computes the NEW lower bounds using arbirtary precision arithmetics
  ! On entry
  !        nl        is the largest size of the Lanczos matrices (for dimensioning only)
  !        kL        is the "reference" Lanczos dimension (i.e. it is 'L+1')
  !        lambda    is the matrix of eigenvalues for dimensions 1,2,..k
  !        sigma_L0  is the variance array of the reference problem 
  !        beta      is the first diagonal of the Lanczos matrix ( beta(1), beta(2),etc. for beta_1, beta_2,etc.)
  !        elower    is a lower bound estimate
  !        ebar1     is the lower bound Ebar_1, not used when qflag = .true.
  !        qboundMP  is an estimate of the Q-ratio
  !        test      instruct the code to use the provided value of qboundMP
  !        qflag     instruct the code to use the new bound to Q0, exploiting the Lanczos construct
  ! On exit
  !        elowernew is the improved estimate for the lower bound
  !        euppernew is the improved estimate for the upper bound
  !        etemple   is the bare Temple lower bound 

  logical(1) :: test, qflag, qfixed, quasiexact
  integer(4) :: fmprec
  integer(4) :: kL, L, nl
  real(8) :: lambda(nl,kL), beta(nl-1) 
  real(8) :: sigma_L0, ebar1, elower, elowernew, euppernew, qbound, rho, eref, etemple
  real(8) :: an, eg
  real(8) :: zero, one, two
  integer(4) :: i, j, n
  real(8), parameter :: sth=1.d-8
  
  TYPE(FM) :: lambdaMP(nl,kL), betaMP(nl-1)
  TYPE(FM) :: gammaMP(kL-1), pMP(kL-1), plessMP
  TYPE(FM) :: sigma_L0MP(nl), ebar1MP, elowerMP, etempleMP, eupperMP, qboundMP
  TYPE(FM) :: zeroMP, oneMP, twoMP

  TYPE(FM) :: lamrefMP, erefMP, alpha0MP
  TYPE(FM) :: ebestMP, deltaMP
  TYPE(FM) :: anMP, bnMP, cnMP, banMP, canMP, aLMP
  TYPE(FM) :: rhoMP, e0MP, egMP, eg1MP, sthMP
  
  L = kL-1                     ! to simplify reading, L is the Lanczos order of the reference problem

  ! Set precision
  call FMSET(fmprec)
  
  ! Convert from double precision to multiple precision
  zero = 0.d0
  one = 1.d0
  two = 2.d0
  CALL FMDP2M(zero, zeroMP)
  CALL FMDP2M(one, oneMP)
  CALL FMDP2M(two, twoMP)
  CALL FMDP2M(sth, sthMP)
  do i = 1, kL
     do j = 1, i
        CALL FMDP2M( lambda(j,i),lambdaMP(j,i) )
     enddo
  enddo

  lamrefMP = lambdaMP(1,kL)      ! this is the approximate ground-state energy from the Lth order problem 

  do i = 1, nl-1
     CALL FMDP2M( beta(i),betaMP(i) )
  enddo
  call FMDP2M( ebar1, ebar1MP)
  call FMDP2M( elower, elowerMP)

  anMP = oneMP
  do n = 1, L
     anMP = anMP + sigma_L0MP(n+1)**2/(elowerMP-lambdaMP(n+1,L+1))**2
  enddo
  if(test)then
     bnMP = sigma_L0MP(1)**2*qboundMP           
     rhoMP = bnMP / anMP                     ! test Eq. 22
     call FMSQRT_R1(rhoMP)
     elowerMP = lamrefMP - rhoMP  
     ! here for Temple's bound
     rhoMP = bnMP
     call FMSQRT_R1(rhoMP)
     etempleMP = lamrefMP - rhoMP
  else
     if(.not.qflag)then
        bnMP = sigma_L0MP(1)**2/(ebar1MP-lamrefMP)
        elowerMP = lamrefMP - bnMP/anMP
        etempleMP =lamrefMP - bnMP
     else
        rhoMP = lamrefMP-elowerMP
        aLMP = rhoMP**2/sigma_L0MP(1)**2*(anMP-oneMP)+oneMP
        if(rhoMP.lt.zeroMP) return
        egMP = lambdaMP(2,kL)-lamrefMP
        eg1MP = egMP                ! this is the most reasonable assumption                 
        if(rhoMP.lt.(aLMP-oneMP)/aLMP*egMP)write(*,*)' TEST!!'
        if(eg1MP.gt.rhoMP)then
           qboundMP = ( (aLMP-oneMP)*(eg1MP-egMP)+ rhoMP )/(eg1MP-rhoMP )        
           bnMP = sigma_L0MP(1)**2*qboundMP 
           rhoMP = bnMP / anMP
           call FMSQRT_R1(rhoMP)
           elowerMP = lamrefMP - rhoMP
           ! here for Temple's bound
           rhoMP = bnMP
           call FMSQRT_R1(rhoMP)
           etempleMP = lamrefMP - rhoMP
        else
           bnMP = sigma_L0MP(1)**2/(ebar1MP-lamrefMP)
           elowerMP = lamrefMP - bnMP/anMP
           etempleMP = lamrefMP - bnMP
        endif
     endif
  endif
  ! convert final result in double precision
  call FMM2DP(elowerMP, elowernew)
  call FMM2DP(etempleMP, etemple)
  
end subroutine bounds_new

subroutine bounds_newex( kL, lambda, nl, variance, lbar1, elower, elowernew, fmprec, etemple )
  use fmzm
  implicit none
  
  ! This subroutine computes the NEW lower bound for the first excited state using arbirtary precision arithmetic
  ! On entry
  !        nl        is the largest size of the Lanczos matrices (for dimensioning only)
  !        kL        is the "reference" Lanczos dimension (i.e. it is 'L+1')
  !        lambda    is the matrix of eigenvalues for dimensions 1,2,..k
  !        variance  is the variance array of the reference problem 
  !        elower    is a lower bound estimate
  !        lbar1     is the lower bound lbar_1
  ! On exit
  !        elowernew is the improved estimate for the lower bound
  !        etemple   is the bare Temple lower bound 

  integer(4) :: fmprec
  integer(4) :: kL, L, nl
  real(8) :: lambda(nl,kL), variance(nl,kL), sigma 
  real(8) :: sigma_L0, lbar1, elower, elowernew, euppernew, qbound, rho, eref, etemple
  real(8) :: an, eg
  real(8) :: zero, one, two
  integer(4) :: i, j, n
  real(8), parameter :: sth=1.d-8
  
  TYPE(FM) :: lambdaMP(nl,kL)
  TYPE(FM) :: sigma_L0MP(nl), lbar1MP, elowerMP, etempleMP, eupperMP, qboundMP
  TYPE(FM) :: zeroMP, oneMP, twoMP

  TYPE(FM) :: lamrefMP, erefMP, alpha0MP
  TYPE(FM) :: ebestMP, deltaMP
  TYPE(FM) :: anMP, bnMP, cnMP, banMP, canMP, aLMP
  TYPE(FM) :: rhoMP, e0MP, egMP, eg1MP, sthMP
  
  L = kL-1                     ! to simplify reading, L is the Lanczos order of the reference problem

  ! Set precision
  call FMSET(fmprec)
  
  ! Convert from double precision to multiple precision
  zero = 0.d0
  one = 1.d0
  two = 2.d0
  CALL FMDP2M(zero, zeroMP)
  CALL FMDP2M(one, oneMP)
  CALL FMDP2M(two, twoMP)
  CALL FMDP2M(sth, sthMP)
  do i = 1, kL
     do j = 1, i
        CALL FMDP2M( lambda(j,i),lambdaMP(j,i) )
     enddo
  enddo
  do j = 1, kL 
     sigma = sqrt(variance(j,kL))              ! for easy-reading, compared to ground-state theory   
     call FMDP2M(sigma, sigma_L0MP(j) )
  enddo
  
  lamrefMP = lambdaMP(2,kL)      ! this is the Ritz first excited-state energy from the Lth order problem 

  call FMDP2M( lbar1, lbar1MP)
  call FMDP2M( elower, elowerMP)

  anMP = sigma_L0MP(1)**2/(lambdaMP(2,L+1)-lambdaMP(1,L+1))**2
  do n = 2, L
     anMP = anMP + sigma_L0MP(n+1)**2/(elowerMP-lambdaMP(n+1,L+1))**2
  enddo

  bnMP = sigma_L0MP(2)**2/(lbar1MP-lamrefMP)
  elowerMP = lamrefMP - bnMP/anMP
  etempleMP =lamrefMP - bnMP 

  ! convert final result in double precision
  call FMM2DP(elowerMP, elowernew)
  call FMM2DP(etempleMP, etemple)
  
end subroutine bounds_newex

