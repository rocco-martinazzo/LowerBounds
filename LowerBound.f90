Program lowerb
  use fmzm

  implicit none

  ! this program reads the matrix elements of a Lanczos matrix of dimensions nl as produced by HPhi
  ! and diagonalize each upperleft subblock of dimension 1,2,..nl
  ! Please refer to the "Notes for Eli", and section 8 in particular for the theory
  
  integer(4), parameter :: maxd=150
  character(1) :: ans 
  character(10) :: answer
  logical(1) :: test, qfixed, scf, qflag, qflag2, qtmp, eread, newbound, tscf
  integer(4) :: nlmax, ni
  integer(4) :: idiag                              ! select diagonalization routine (for testing purposes)
  integer(4) :: icyc, jcyc, icycx, jcycx, jscf
  real(8), allocatable :: alpha0(:), beta0(:)      ! "0" is appended here to not conflict with high precision module
  real(8), allocatable :: alphanew(:), betanew(:), psi0(:)
  real(8), allocatable :: ulanc(:,:), usave(:), usavex(:)
  real(8), allocatable :: vlanc(:,:)
  real(8), allocatable :: lambda(:,:)              ! the columns of this matrix are the Lanczos eigenvalues at 1st, 2st, etc.  order
  real(8), allocatable :: variance(:,:)            ! and this contains the corresponding variances
  real(8), allocatable :: lambdarest(:)
  real(8), parameter :: sth = 1.d-16               ! parameter controlling the self-consistent cycle

  ! scratch & temporary variables and arrays
  integer(4) :: i, j, k, n, lwork, liwork, info, io, lz, worki, lunit
  real(8) :: workl, csi, am, eref, erefx, small, ebarold
  real(8) :: an, bn, cn, anold, bnold, cnold, ban, can, banold, canold, eba, etmp
  real(8) :: sigma_L0
  real(8) :: var1, var2, ebar1, ebar2, etemple1, etemple2, elower, elowernew, euppernew, aratio, rad, u(2), v
  real(8) :: elsave0, lbar1, wb1, wb2, wb3, den, dl
  real(8) :: elowerx, elowernewx
  real(8) :: euppernew2, elowernew2, rad2, elower2
  real(8) :: qbound, aex, qex, qcheck, ebar1ex, ebar2ex, ebar2save, factor, ebar1app, ebar2app, gap
  
  real(8), allocatable :: epsu(:), epsl(:), epslx(:), epsm(:), qb(:), al(:), rho(:), ebar(:), etemple(:), etemplex(:)
  real(8), allocatable :: epsu2(:), epsl2(:), qb2(:), rho2(:)
  real(8), allocatable :: epstight(:,:)
  integer(4), parameter :: jscfmax = 5
  real(8), parameter :: eps=1.d-8
  integer(4), parameter :: icycmax=200
  real(8), allocatable :: work(:), betaprod(:)
  integer(4), allocatable :: iwork(:)
  
  ! this is for high precision arithmetics (MUST be the same precision as in bounds)
  real(8) :: zero, one, two
  integer(4), parameter :: fmprec=64
  TYPE(FM) :: zeroMP, oneMP, twoMP
  TYPE(FM) :: qboundMP, qboundMP2, uMP, betaMP
  TYPE(FM) :: qcheckMP, erefMP, lambdaMP
  TYPE(FM) :: factorMP
  TYPE(FM), allocatable :: sigma_L0MP(:)
     
  ! Set precision for higher precision arithmetics
  call FMSET(fmprec)
  ! Convert from double precision to multiple precision
  zero = 0.d0
  one = 1.d0
  two = 2.d0
  CALL FMDP2M(zero, zeroMP)
  CALL FMDP2M(one, oneMP)
  CALL FMDP2M(two, twoMP)
  ! set flag for NEW bound
  newbound=.false.
  ! set flag for Q bound
  qflag=.false.
  ! set flag for Q bound
  qflag2=.false.
  ! set flag for SCF
  scf = .false.
  ! set flag for Tight SCF
  tscf = .false.
  ! set flag for testing
  test = .false. 

  write(*,*)' *********************************** '
  write(*,*)' *                                 * '
  write(*,*)' *   LLB (Lanczos Lower bounds)    * '
  write(*,*)' *                                 * '
  write(*,*)' *********************************** '
  
  open(file='alpha.dat',unit=11)
  write(*,*) ' Opening file alpha.dat for reading max dim '
  nlmax =0 
  do
     read(11,*,iostat=io) 
     if(io/=0) exit
     nlmax = nlmax + 1
  enddo
  close(11)
  write(*,*) ' Full size of the Lanczos matrix = ', nlmax
  
  open(file='alpha.dat',unit=11)
  open(file='beta.dat',unit=12)
  write(*,*) ' Opening file alpha.dat for reading '
  write(*,*) ' Opening file  beta.dat for reading '

  open(file='ebar.dat',unit=16)
  open(file='ebarx.dat',unit=116)
  open(file='residual.dat',unit=19)
  open(file='residual_energies.dat',unit=119)
  open(file='upper.dat',unit=20)
  open(file='lower.dat',unit=21)
  open(file='lower_weinstein.dat',unit=22)
  open(file='errors.dat',unit=24)
  open(file='ratios.dat',unit=124)
  open(file='errorsx.dat',unit=224)
  open(file='ratiosx.dat',unit=324)
  open(file='gap.dat',unit=424)
  write(*,*) ' Opening file ebar.dat                for writing' 
  write(*,*) ' Opening file upper.dat               for writing'    ! this contains the lowest eigenvalues 
  write(*,*) ' Opening file lower.dat               for writing'    ! this contains the computed lower bound to the ground-state energy 
  write(*,*) ' Opening file lower_weinstein.dat     for writing'    ! this contains the simple Weinstein bounds for the lowest eigenvalues
  write(*,*) ' Opening file errors.dat              for writing'    ! these are for the (better) upper bound given by the approximate eigenstate 
  write(*,*) ' Opening file ratios.dat              for writing'
  allocate(alpha0(nlmax))
  allocate(beta0(nlmax-1))
  allocate(epsu(nlmax))
  allocate(epsl(nlmax))
  allocate(epslx(nlmax))
  allocate(epsm(nlmax))
  allocate(qb(nlmax))
  allocate(rho(nlmax))
  allocate(al(nlmax))
  allocate(ebar(nlmax))
  allocate(etemple(nlmax))
  allocate(etemplex(nlmax))
  !
  allocate(epsu2(nlmax))
  allocate(epsl2(nlmax))
  allocate(rho2(nlmax))
  allocate(qb2(nlmax))
  allocate(epstight(nlmax,jscfmax))
  !
  ! if(test)allocate(usave(nlmax))
  allocate(usave(nlmax))                ! save ground- and first-excited state vectors
  allocate(usavex(nlmax))               ! in any case
  
  ! notice that in the notes we introduce the Lanczos order L = n-1 and Lmax = nlmax - 1
  ! and consider alpha0(0),..alpha0(L) (differently from below)
  !              beta0(1),..beta0(L)   (same as below)
  
  read(11,*) alpha0(1) 
  do i = 1, nlmax-1
     read(11,*) alpha0(i+1)
     read(12,*) beta0(i)
  enddo
  
  allocate(ulanc(nlmax,nlmax))
  allocate(lambda(nlmax,nlmax))
  allocate(variance(nlmax,nlmax))
  allocate(betaprod(nlmax))
  allocate(vlanc(nlmax,nlmax))
  allocate(lambdarest(nlmax))
  
  ! set dimensions of work arrays at its optimal for nl

  idiag = 2                                 ! idiag = 1, 2 for using dsyev or dsyevd throughout the code
                                            !       There is probably no real difference between the two
                                            !       Preliminary results showing some effects were probably
                                            !       affected by an erroneous query call to dsyev
                                            !       Set it to 2 for a safe use
  liwork = 1
  if(idiag.eq.1) then
     call dsyev('V', 'U', nlmax, ulanc, nlmax, lambda(1,1), workl, -1, info)
  elseif(idiag.eq.2)then
     call dsyevd('V', 'U', nlmax, ulanc, nlmax, lambda(1,1), workl, -1, worki, -1, info)
  else
     write(*,*) ' Error! Select a valid diagonalization routine! '
     stop
  endif
  if(info.eq.0) then
     lwork = int(workl)
     liwork = worki
  else
     write(*,*)' *** Check your code/data !! ***'
     write(*,*)' *** Something went wrong in dimensioning !'
  end if
  allocate(work(lwork))
  allocate(iwork(liwork))
  
  ! Diagonalize subblocks
  
  write(*,*) ''
  write(*,*) ' ******  Diagonalize subblocks '
  write(*,*) ' '
  lambda(1,1) = alpha0(1)

  ebarold = 1.d+10
  elowernew = -1.0d+2
  elowernewx = elowernew
  elower = 0.0d0
  elowerx = elower
  
  ! **************
  ! this is for testing purpose: perform diagonalization of the larger matrix available
  ulanc = 0.0d0
  do i = 1, nlmax
     ulanc(i, i) = alpha0(i)
     if(i+1.le.nlmax) ulanc(i, i+1) = beta0(i)
  enddo
  if(idiag.eq.1) then
     call dsyev('V', 'U', nlmax, ulanc, nlmax, lambda(1,nlmax), work, lwork, info)
  elseif(idiag.eq.2)then
     call dsyevd('V', 'U', nlmax, ulanc, nlmax, lambda(1,nlmax), work, lwork, iwork, liwork,  info)
  endif
  write(*,*) ''
  write(*,*) ' ***** Reference results ***** ' 
  write(*,*) ''
  if(info == 0) then
     write(*,*) ' SubL = ', nlmax-1, ' E(i) = ', (lambda(i,nlmax), i = 1,min(6,nlmax))
  elseif(info < 0 ) then
     write(*,*) ' *** Error in dsyev !! **** '
     write(*,*) ' '
     write(*,*) ' *** The ',-info,'-th argument had an illegal value'
  else
     write(*,*) ' *** Error in dsyev !! **** '
     write(*,*) ' '
     write(*,*) ' *** The algorithm failed to converge;', info
     write(*,*) ' *** off-diagonal elements of an intermediate tridiagonal'
     write(*,*) ' *** form did not converge to zero.'
  endif

  eref = lambda(1,nlmax)   ! this is the best estimate of the ground-state energy we have
  erefx = lambda(2,nlmax)  ! this is the best estimate of the first excited state energy we have
  
  ! compute Q0 using high-precision arithmetics
  
  call FMDP2M( ulanc(1,1), uMP) 
  qboundMP = (oneMP-uMP**2)/uMP**2                  ! this is the exact Q_0 in quadruple precision
  call FMM2DP(qboundMP,qex)                         ! this is the exact Q_0 in double precision
  
  ! this is the dp value of the exact Ebar1 (in the new notation this is "\bar(lambda)_0")
  ebar1ex  = ( alpha0(1)-eref ) / ( 1.d0-ulanc(1,1)**2 ) + eref
  write(*,*) ' Exact Energy, Ebar1, Q0 = ', eref, ebar1ex, qex
  write(*,*) ''
  write(*,*) ' ******************** ' 
  write(*,*) ''

  ! compute spectrum of residual Hamiltonian for L=0
  vlanc = 0.0d0
  do i = 2, nlmax
     vlanc(i, i) = alpha0(i)
     if(i+1.le.nlmax) vlanc(i, i+1) = beta0(i)
  enddo
  if(idiag.eq.1) then
     call dsyev('N', 'U', nlmax, vlanc, nlmax, lambdarest, work, lwork, info)
  elseif(idiag.eq.2)then
     call dsyevd('N', 'U', nlmax, vlanc, nlmax, lambdarest, work, lwork, iwork, liwork,  info)
  endif
  if( info < 0 ) then
     write(*,*) ' *** Error in dsyev !! **** '
     write(*,*) ' '
     write(*,*) ' *** The ',-info,'-th argument had an illegal value'
  elseif( info > 0) then
     write(*,*) ' *** Error in dsyev !! **** '
     write(*,*) ' '
     write(*,*) ' *** The algorithm failed to converge;', info
     write(*,*) ' *** off-diagonal elements of an intermediate tridiagonal'
     write(*,*) ' *** form did not converge to zero.'
  endif
  
  write(*,*) ' We are almost ready to compute lower bound estimates..'
  write(*,*) '   Selected newbound estimate (Y/N)?'
  read(*,*) answer
  if (trim(answer)=="Y" )  newbound = .true.
  if(.not. newbound) then
     ni = 1
     write(*,*) '   Select Q-bound estimate' 
     write(*,*) '   Type <bare|improved> for bare or improved overlap  '
     read(*,*) answer
     if (trim(answer)=="bare" )  qflag = .true.
          
     ! if qfixed = .true.  uses energy-independent qbound  
     !           = .false. uses energy-dependent   qbound
     !     test  = .true. for using exact Q0 (enforces qfixed = .true. )
     !      scf  = .true. instruct the code to reach self-consistency
     qfixed=.false.
     write(*,*) '   Select energy dependent Q-bound ' 
     write(*,*) '   Type <T|F> for energy dependent or independent Q bound'
     read(*,*) qfixed
     open(file='errors_Lanczos.dat',unit=23)
     open(file='qbounds.dat',unit=25)
     write(*,*) ''
     write(*,*) ' Opening file errors_Lanczos.dat      for writing'    ! these are for the both the upper and the lower bounds provided by Lanczos
     write(*,*) ' Opening file qbounds.dat             for writing'    ! this contains the bound to the effective Q-ratio used in the calculation
     write(*,*) ''
  else
     ni = nlmax
     write(*,*) '   Select Q-bound estimate' 
     write(*,*) '   Type <new|old> for new or old estimate  '
     read(*,*) answer
     if (trim(answer)=="new" )  qflag2 = .true.          ! we use qflag2 for the new qbound estimate that exploit Lanczos construct in e1bar estimate
     ! if(newbound.and.test) usave(:) = ulanc(:,1)       ! when test=.true. save the exact ground-state vector
     if(newbound) usave(:)  = ulanc(:,1)                 ! save ground- and first-excited vectors 
     if(newbound) usavex(:) = ulanc(:,2)                 ! anyway when newbound=.true.
  endif
  write(*,*) '   Ask for self-consistency ' 
  write(*,*) '   Type <T|F> for self-consistent or bare estimate'
  read(*,*) scf
  write(*,*) '   Ask for tight self-consistency'                       ! this assumes that the Ritz error is the smallest 
  write(*,*) '   Type <T|F> for tight self-consistent estimate'
  read(*,*) tscf
  if(tscf)then
     open(file='errors_tbound.dat',unit=26)
     open(file='terrors_scf.dat',unit=28)
     write(*,*) ' Opening file errors_tbound.dat       for writing'    ! these are for the (better) upper bound given by the approximate eigenstate 
     write(*,*) ' Opening file terrors_scf.dat         for writing'    ! save tbound iterative results for the lower bound
     if(.not.newbound)then
        open(file='qbounds_tbound.dat',unit=27)
        write(*,*) ' Opening file qbounds_tbound.dat      for writing'    ! this contains the bound to the effective Q-ratio used in the calculation
     endif
  endif

  allocate(sigma_L0MP(ni))

  write(*,*) ' Current status of flags newbound | qflag | qflag2 | QE | SC =', newbound, qflag, qflag2, qfixed, scf
  write(*,*) ' '
  if( test ) qfixed = .true.
  if( test ) scf = .false.
  if( test ) qflag = .false.
  if( test ) qflag2 = .false.
  if(test) then
     write(*,*) ''
     write(*,*) ' !!!!! This is a test calculation using exact Q0 !!!!!!'
     write(*,*) ''
  endif

  ! ************* add this for printing purposes only
  n=1
  ulanc(1,1) = 1.d0
  icyc = 1
  jscf = 0
  ebar1 = lambda(1,1) + eps
  if(newbound)then   
     v = usave(1)
     v=min(abs(v),1.0d0)
     call FMDP2M( v, uMP) 
     qboundMP = (oneMP-uMP**2)/uMP**2                  ! this is the exact Q_0 in quadruple precision
     call FMM2DP( qboundMP, qbound)
     ebar1ex  = lambda(1,n) + (lambda(1,n)-eref)/ qbound
  endif
  ebar(n) = ebar1
  write(*,*) ' SubL = ', n-1, ' E(i) = ', (lambda(i,n), i = 1,min(6,n))
  do i = 1, n
     variance(i, n) = abs( beta0(n)*ulanc(n,i) )**2
  enddo
  write(*,*) ''
  write(*,*) '                    StE(i) = ', ( abs( beta0(n)*ulanc(n,i) ), i = 1,min(6,n) )
  write(*,*) ' ResL = ', n-1, ' E(i) = ', (lambdarest(i), i = 1,min(6,nlmax-n))
  write(*,*) ' Lambda_0 (bar)        = ', ebar1
  if(newbound) write(*,*) ' Lambda_0 (bar) exact     = ', ebar1ex  
  write(*,*) ' E_lb (in)          = ', elower
  write(*,*) ' E_lb (out)         = ', elowernew
  if(scf)then
     write(*,*) ' # of self-cycles   = ', icyc-1
  endif
  if(tscf)then
     write(*,*) ' # of tself-cycles  = ', jscf-1
     write(*,*) ' E_lb (tight)       = ', epsl2(n)
  endif
  write(*,*) ' '
  write(16,*) n-1, ebar1
  write(19,*) n-1, (lambdarest(i),i=1,6)
  write(20,*) n-1, (lambda(i,n), i=1,6)
  write(21,*) n-1, elowernew
  write(22,*) n-1, (lambda(i,n)-abs(beta0(n)*ulanc(n,i)), i=1,6 )
  ! *****************

  do n = 2, min(nlmax,maxd)   ! n is the problem dimension, the Lanczos order L is given by L = n-1
     icyc = 1
     icycx = 1
     jscf = 0
     ulanc = 0.0d0
     vlanc = 0.0d0
     do i = 1, n
        ulanc(i, i) = alpha0(i)
        if(i+1.le.nlmax) ulanc(i, i+1) = beta0(i)
     enddo
     do i = n+1, nlmax
        vlanc(i, i) = alpha0(i)
        if(i+1.le.nlmax) vlanc(i, i+1) = beta0(i)
     enddo
     ! compute the spectrum of the projected Hamiltonian
     if(idiag.eq.1)then
        call dsyev('V', 'U', n, ulanc, nlmax, lambda(1,n), work, lwork, info)
     elseif(idiag.eq.2)then
        call dsyevd('V', 'U', n, ulanc, nlmax, lambda(1,n), work, lwork, iwork, liwork, info)
     endif
     if(info == 0) then
        write(*,*) ' ------------------------------------------------------------------------ ', n-1, &
                   ' ------------------------------------------------------------------------ '
        write(*,*) ' '
        write(*,*) ' SubL = ', n-1, ' E(i) = ', (lambda(i,n), i = 1,min(6,n))
     elseif(info < 0 ) then
        write(*,*) ' *** Error in dsyev !! **** '
        write(*,*) ' '
        write(*,*) ' *** The ',-info,'-th argument had an illegal value'
     else
        write(*,*) ' *** Error in dsyev !! **** '
        write(*,*) ' '
        write(*,*) ' *** The algorithm failed to converge;', info
        write(*,*) ' *** off-diagonal elements of an intermediate tridiagonal'
        write(*,*) ' *** form did not converge to zero.'
     endif
     ! compute the spectrum of the residual Hamiltonian (notice we keep it of the same max dimension)
     if(idiag.eq.1)then
        call dsyev('N', 'U', nlmax, vlanc, nlmax, lambdarest, work, lwork, info)
     elseif(idiag.eq.2)then
        call dsyevd('N', 'U', nlmax, vlanc, nlmax, lambdarest, work, lwork, iwork, liwork, info)
     endif
     if(info == 0) then
        write(*,*) ' ResL = ', n-1, ' E(i) = ', (lambdarest(i), i = 1,min(6,nlmax-n))
     elseif(info < 0 ) then
        write(*,*) ' *** Error in dsyev !! **** '
        write(*,*) ' '
        write(*,*) ' *** The ',-info,'-th argument had an illegal value'
     else
        write(*,*) ' *** Error in dsyev !! **** '
        write(*,*) ' '
        write(*,*) ' *** The algorithm failed to converge;', info
        write(*,*) ' *** off-diagonal elements of an intermediate tridiagonal'
        write(*,*) ' *** form did not converge to zero.'
     endif
     ! compute eigenstate variance, using high precision arithmetics
     do i = 1, ni
        call FMDP2M( ulanc(n,i), uMP )
        call FMDP2M( beta0(n), betaMP)         
        sigma_L0MP(i) = abs(betaMP*uMP)
     enddo
     ! when test = .true. compute exact Q_0 for the new bound estimate
     ! if(newbound.and.test)then
     ! compute exact Q_0 and residual energies always when newbound=.true. 
     if(newbound) then
        ! this is for the first excited state
        v = 0.0d0
        do i = 1, n
           v = v + usavex(i)*ulanc(i,2)
        enddo
        v=min(abs(v),1.0d0)
        call FMDP2M( v, uMP) 
        qboundMP = (oneMP-uMP**2)/uMP**2                  ! this is the exact Q_0 in quadruple precision
        call FMM2DP( qboundMP, qbound)
        qbound = max(qbound,1.d-16)
        ebar2ex  = lambda(2,n) + (lambda(2,n)-erefx)/ qbound
        ! write(*,*) 'Excited', qbound, v**2
        ! this is for the ground-state
        v = 0.0d0
        do i = 1, n
           v = v + usave(i)*ulanc(i,1)
        enddo
        v=min(abs(v),1.0d0)
        call FMDP2M( v, uMP) 
        qboundMP = (oneMP-uMP**2)/uMP**2                  ! this is the exact Q_0 in quadruple precision
        call FMM2DP( qboundMP, qbound)
        qbound = max(qbound,1.d-16)
        ebar1ex  = lambda(1,n) + (lambda(1,n)-eref)/ qbound
        ! write(*,*) 'Ground', qbound, v**2
        write(119,*) n-1, ebar1ex, ebar2ex 
     endif
     if( n.lt. nlmax) then
        ! lower bound estimate for the first excited state: it can be larger than E1!
        ebar2 = lambda(2,n)-abs(beta0(n)*ulanc(n,2))
        if(ebar2.lt.lambda(1,n)) then   !!! To be carefully checked for crossing  !!
           if ( qflag ) then
              ebar1 = lambda(1,1) + eps
           else
              ebar1 = lambda(1,n) + eps
           endif
           ebarold = ebar2
        endif  
        if(ebar2.lt.ebarold .and. ebar1.lt.ebarold) then
           ! if qflag take ebar1 infinitesimally larger than alpha0, otherwsise larger than E0
           if ( qflag ) then
              ebar1 = lambda(1,1) + eps
           else
              ebar1 = lambda(1,n) + eps
           endif
           ebarold = ebar2
        else
        ! take best Weinstein lower bound estimate for the first excited state  
           if( ebar2.gt.lambda(1,n) ) ebar1 = max(ebar2, ebar1)
        endif
        
        
10      continue
        if(elowernew.gt.elower) elower = elowernew
        if(elower.gt.lambda(1,n)) elower = lambda(1,n) -abs(beta0(n)*ulanc(n,1))
        ! compute bound with high precision (see precision set in routine bounds)
        if(newbound)then
           qtmp = .false.
           if(n.gt.15) qtmp = qflag2

           ! !! test exact
           ! ebar1 = ebar1ex
           
           call bounds_new( n, lambda, nlmax, sigma_L0MP, beta0, ebar1, elower, elowernew, & 
                qboundMP, test, qtmp, eref, fmprec, etemple(n) )
           icyc = icyc + 1
           if( scf .and. (elowernew-elower.gt.sth) .and. icyc.lt.icycmax ) goto 10
        else
           if ( qflag )then
              if( .not. qfixed ) scf = .false.        ! nothig to make self-consistent if qfixed is false
              call qbounds( n, lambda, nlmax, sigma_L0MP, beta0, ebar1, elower, elowernew, euppernew, & 
                   epsm(n), qboundMP, rad, qfixed, test, eref, fmprec )
           else
              call bounds( n, lambda, nlmax, sigma_L0MP, beta0, ebar1, elower, elowernew, euppernew, & 
                   epsm(n), qboundMP, rad, qfixed, test, eref, fmprec )
           endif
           icyc = icyc + 1
           if(rad.lt.0.d0) then
              write(*,*) ' Warning!! Rho =', rad,' < 0.d0 !! &
                   Resetting lower bound estimate for self-consistency'
              ! this ensures that elower at the next step is Temple's lower bound
              elower = lambda(1,n) - 2.0d0* abs(beta0(n)*ulanc(n,1))**2/( ebar1-lambda(1,n) ) 
              elowernew = elower +  abs(beta0(n)*ulanc(n,1))**2/( ebar1-lambda(1,n) ) 
           endif
           if( scf .and. (elowernew-elower.gt.sth) .and. icyc.lt.icycmax .and. rad.gt.0.0d0 ) goto 10
        endif

        elower2 = elowernew
        ! this mught be a very bad estimate for the first n
        if( elowernew.lt.elower ) elower2 = elower 

        if(tscf)then
           if(newbound) then
              elowernew2 = elowernew
              jscf = 0
120           continue
              elowernew2 = (elowernew2 + lambda(1,n))/2.d0
110           continue
              jcyc = 1
              elower2 = elowernew2
              call bounds_new( n, lambda, nlmax, sigma_L0MP, beta0, ebar1, elower2, elowernew2, & 
                   qboundMP, test, qtmp, eref, fmprec, etemple(n) )
              jcyc = jcyc + 1
              if( scf .and. (elowernew2-elower2.gt.sth) .and. jcyc.lt.icycmax ) goto 110
              elowernew2 = max(elower2, elowernew2)
              !write(*,*) jscf, elowernew, elowernew2, lambda(1,n) 
              jscf = jscf + 1
              if( jscf.le.jscfmax ) then
                 epstight(n, jscf) = elowernew2
                 ! save first tbound results
                 if(jscf.eq.1) then
                    ! this is for the tbound
                    epsl2(n) = elowernew2
                 endif
                 goto 120
              endif
           else
              if( .not. qflag) then           
20               call tbounds( n, lambda, nlmax, sigma_L0MP, beta0, ebar1, elower2, elowernew2, euppernew2, & 
                      qboundMP2, rad2, qfixed, fmprec )
                 jscf = jscf + 1
                 if (rad2.lt.0.0d0) write(*,*) ' TBound warning!! Rho =', rad2,' < 0.0d0 !!'
                 if( jscf.le.jscfmax ) then
                    epstight(n, jscf) = elowernew2
                    elower2 = max(elowernew2,elowernew)
                    ! save first tbound results
                    if(jscf.eq.1) then
                       ! this is for the tbound
                       epsu2(n) = euppernew2
                       epsl2(n) = elowernew2
                       call FMM2DP(qboundMP2,qb2(n))    ! this is an upper bound to the effective Q-ratio used 
                       rho2(n) = rad2
                    endif
                    goto 20
                 endif
              endif
           endif
        endif
        ebar(n) = ebar1
        do i = 1, n
           variance(i, n) = abs( beta0(n)*ulanc(n,i) )**2 ! save variances
        enddo
        write(*,*) '                    StE(i) = ', ( abs( beta0(n)*ulanc(n,i) ), i = 1,min(6,n) )
        write(*,*) ''
        write(*,*) ' Lambda_0 (bar)     = ', ebar1
        if(newbound) write(*,*) ' Lambda_0 (bar) exact = ', ebar1ex
        write(*,*) ' E_lb (in)          = ', elower
        write(*,*) ' E_lb (out)         = ', elowernew
        if(scf) then
           write(*,*) ' # of self-cycles   = ', icyc-1
        endif
        if(tscf) then
           write(*,*) ' # of tself-cycles  = ', jscf-1
           write(*,*) ' E_lb (tight)       = ', epsl2(n)
        endif
        write(16,*) n-1, ebar1
        write(19,*) n-1, (lambdarest(i),i=1,6)
        write(20,*) n-1, (lambda(i,n), i=1,6)
        write(21,*) n-1, elowernew
        write(22,*) n-1, (lambda(i,n)-abs(beta0(n)*ulanc(n,i)), i=1,6 )
        epsu(n) = euppernew
        epsl(n) = elowernew
        call FMM2DP(qboundMP,qb(n))    ! this is an upper bound to the effective Q-ratio used 
        rho(n) = rad                   ! output "rad" in qbounds.dat

        ! let's now compute the lower bound to the first excited-state energy for L-1 (only for newbound=.true.)
        if(newbound.and.n.gt.3) then
           ! compute \bar(lambda)_1 for L-1. This can be improved when a better estimate for the lower bound
           ! is available (just replace wb1 with the new bound)
           wb1 = lambda(2,n-1) - sqrt(variance(2,n-1))
           wb2 = lambda(3,n-1) - sqrt(variance(3,n-1))
           wb3 = lambda(4,n-1) - sqrt(variance(4,n-1))

           ! "standard" M = 2
           den = variance(2, n-1)/variance(1,n-1)*(lambda(1,n-1)-elsave0)**2/(lambda(2,n-1)-lambda(1,n-1))**2 * &
                (lambda(3, n-1) - elsave0)
           dl = max( lambda(2,n-1)-lambda(2,n), 1.0d-15)
           ! test
           ! dl = max( lambda(2,n-1)-erefx, 1.0d-15)
           den  =  1.d0 + den / dl             
           lbar1 = wb1 + ( wb2 - lambda(2,n-1) ) / den
           ! write(*,*) ' ********** 1', den
           lbar1 = max(lbar1, lambda(2,n-1))
           
           ! test "improved expression" M = 2
           ! den = variance(2, n-1)/variance(1,n-1)*(lambda(1,n-1)-elsave0)**2/(lambda(2,n-1)-lambda(1,n-1))**2 * &
           !     (lambda(3, n-1) - elsave0)
           ! den = den + &
           !     variance(2, n-1)/variance(3,n-1) * (lambda(3,n-1)-wb2)**2/(lambda(2,n-1)-wb2)**2 * &
           !     (lambda(3, n-1) - wb2)
           ! dl = max( lambda(2,n-1)-lambda(2,n), 1.0d-15)
           ! den = 1.0d0 + den / dl
           ! lbar1 = wb1 + ( lambda(3,n-1) - lambda(2,n-1) ) / den
           ! write(*,*) ' ********** 2', den
           ! lbar1 = max(lbar1, lambda(2,n-1))
           !
           ! test "improved expression" M = 3
           ! den = variance(2, n-1)/variance(1,n-1)*(lambda(1,n-1)-elsave0)**2/(lambda(2,n-1)-lambda(1,n-1))**2 * &
           !     (lambda(4, n-1) - elsave0)
           ! den = den + &
           !     variance(2, n-1)/variance(3,n-1) * (lambda(3,n-1)-wb2)**2/(lambda(2,n-1)-wb2)**2 * &
           !     (lambda(4, n-1) - wb2)
           ! den = den + &
           !     variance(2, n-1)/variance(4,n-1) * (lambda(4,n-1)-wb3)**2/(lambda(2,n-1)-wb3)**2 * &
           !     (lambda(4, n-1) - wb3)
           ! dl = max( lambda(2,n-1)-lambda(2,n), 1.0d-15)
           ! den = 1.0d0 + den / dl
           ! lbar1 = wb1 + ( lambda(4,n-1) - lambda(2,n-1) ) / den
           ! write(*,*) ' ********** 3', den
           ! lbar1 = max(lbar1, lambda(2,n-1))
           !
           ! test "standard" M = 3
           ! den = variance(2, n-1)/variance(1,n-1)*(lambda(1,n-1)-elsave0)**2/(lambda(2,n-1)-lambda(1,n-1))**2 * &
           !     (lambda(4, n-1) - elsave0)
           ! den = den + &
           !     variance(2, n-1)/variance(3,n-1) * (lambda(3,n-1)-wb2)**2/(lambda(2,n-1)-wb2)**2 * &
           !     (lambda(4, n-1) - wb2)
           ! dl = max( lambda(2,n-1)-lambda(2,n), 1.0d-15)
           ! den = 1.0d0 + den / dl
           ! lbar1 = wb1 + ( wb3 - lambda(2,n-1) ) / den
           ! write(*,*) ' ********** 4', den
           ! lbar1 = max(lbar1, lambda(2,n-1))
           
           write(116,*) n-2, lbar1
           ! we are now ready to compute the bound..!
100        continue
           if(elowernewx.gt.elowerx) elowerx = elowernewx
           if(elowerx.gt.lambda(2,n)) elowerx = lambda(2,n) -abs(beta0(n)*ulanc(n,2))
           ! compute bound with high precision (see precision set in routine bounds)
           ! !! test exact !!
           ! lbar1 = ebar2save
           
           call bounds_newex( n-1, lambda, nlmax, variance, lbar1, elowerx, elowernewx, fmprec, etemplex(n-1) )
           icycx = icycx + 1
           if( scf .and. (elowernewx-elowerx.gt.sth) .and. icycx.lt.icycmax ) goto 100
           epslx(n-1) = elowernewx
           write(*,*) ' --- First excited state L =  ', n-2
           write(*,*) ' Lambda_1 (bar)     = ', lbar1
           write(*,*) ' E_lb (in)          = ', elowerx
           write(*,*) ' E_lb (out)         = ', elowernewx
           if(scf) then
              write(*,*) ' # of self-cycles   = ', icycx-1
           endif
           ! wb1 = elowernewx
        endif
        ! save best ground-state lower bound to compute the first excited-state lower bound
        ! at the next step. For the first and second excited state we use Wienstein bounds to start with
        if(newbound.and.n.ge.3) elsave0 = elowernew

        !!! test
        if(newbound.and.n.ge.3) ebar2save = ebar2ex
        
        write(*,*) ' '
     endif
  enddo
  
  do n = 2, min(nlmax,maxd) -1
     write(24,*) n-1, eref-epsl(n),  lambda(1,n)-eref, 0.5d0*( lambda(1,n)+epsl(n) )-eref         ! errors.dat
     write(124,*) lambda(1,n)-eref, abs(eref-epsl(n))/ max(abs((lambda(1,n)-eref)),sth*1.d-1),  & ! ratios.dat
                                    abs(eref-etemple(n))/ max(abs((lambda(1,n)-eref)),sth*1.d-1)
  enddo
  if(newbound) then
     gap =  erefx-eref
     do n = 3, min(nlmax-1,maxd) -1
        write(224,*) n-1, erefx-epslx(n),  lambda(2,n)-erefx,  erefx-etemplex(n)                           ! errorsx.dat
        write(324,*) lambda(2,n)-erefx, abs(erefx-epslx(n))/ max(abs((lambda(2,n)-erefx)),sth*1.d-1),  &   ! ratiosx.dat
                                   abs(erefx-etemplex(n))/ max(abs((lambda(2,n)-erefx)),sth*1.d-1)
        write(424,*) n-1, lambda(2,n)-epsl(n) - gap, epslx(n)-lambda(1,n)-gap, lambda(2,n)-lambda(1,n)-gap ! gap.dat
     enddo
  endif
  close(16)
  close(116)
  close(19)
  close(20)
  close(21)
  close(22)
  close(24)
  if(tscf)then
     do n = 2, nlmax -1
        write(26,*) n-1, eref-epsl2(n), lambda(1,n)-eref, 0.5d0*( lambda(1,n)+epsl2(n) )-eref  ! errors-TBOUND.dat
        write(28,*) n-1, (eref-epstight(n,jscf), jscf = 1,jscfmax)                             ! terrors_scf.dat
     enddo
     close(26)
     close(28)
     if(.not.newbound)then
        do n = 2, nlmax-1 
           write(27,*) n-1, qb2(n), qex, qb2(n)-qex, rho2(n)                                   ! qbounds-TBOUND.dat
        enddo
        close(27)
     endif
  endif
  if(.not.newbound)then
     do n = 2, nlmax-1
        write(23,*) n-1, eref-epsl(n), epsu(n)-eref, epsm(n)-eref                              ! errors_lanczos.dat
        write(25,*) n-1, qb(n), qex, qb(n)-qex, rho(n)                                         ! qbounds.dat
     enddo
     close(23)
     close(25)
  endif


  ! restart Lanczos iteration from a given intermediate eigenstate

  ! set TEST flag
  test = .false.

  write(*,*) ''
  write(*,*) ''
  write(*,*) ' Target step (L value)? [A negative value stops the code]' 
  read(*,*) n
  if( n.le.0 ) stop
  n = n + 1
  write(*,*)
  write(*,*) ' Restart Lanczos iterations using the intermediate eigenstate of the', n-1, ' step '
  write(*,*) ''
  open(file='ebar-phi0.dat',unit=16)
  lunit=17
  open(file='spectrum-phi0.dat',unit=lunit)       ! this is passed to Lanczos below
  open(file='coeffs-phi0.dat',unit=18)
  open(file='residual-phi0.dat',unit=19)
  open(file='upper-phi0.dat',unit=20)
  open(file='lower-phi0.dat',unit=21)
  open(file='lower_weinstein-phi0.dat',unit=22)
  open(file='errors-phi0.dat',unit=24)
  if(.not.newbound) then
     open(file='errors_Lanczos-phi0.dat',unit=23)
     open(file='qbounds-phi0.dat',unit=25)
  endif

  write(*,*) ' Opening file ebar-phi0.dat               for writing'
  write(*,*) ' Opening file coeffs-phi0.dat             for writing'
  write(*,*) ' Opening file spectrum-phi0.dat           for writing'
  write(*,*) ' Opening file upper-phi0.dat              for writing'
  write(*,*) ' Opening file lower-phi0.dat              for writing'
  write(*,*) ' Opening file lower_weinstein-phi0.dat    for writing'    ! this contains the simple Weinstein bounds for the lowest eigenvalues
  write(*,*) ' Opening file errors-phi0.dat             for writing'    ! these are for the (better) upper bound given by the approximate eigenstate 
  if(.not.newbound) then
     write(*,*) ' Opening file errors_Lanczos-phi0.dat     for writing'    ! these are for the both the upper and the lower bounds provided by Lanczos
     write(*,*) ' Opening file qbounds-phi0.dat            for writing'    ! this contains the bound to the effective Q-ratio used in the calculation
  endif
  write(*,*)''
  
  n = min(n,nlmax)
  ulanc = 0.0d0
  do i = 1, n-1
     ulanc(i, i) = alpha0(i)
     ulanc(i, i+1) = beta0(i)
  enddo
  ulanc(n,n) = alpha0(n)
  if(idiag.eq.1)then
     call dsyev('V', 'U', n, ulanc, nlmax, lambda(1,n), work, lwork, info)
  elseif(idiag.eq.2)then
     call dsyevd('V', 'U', n, ulanc, nlmax, lambda(1,n), work, lwork, iwork, liwork, info)
  endif
  if(info == 0) then
     write(*,*) ' Ground-state energy at this step ', lambda(1,n)
  elseif(info < 0 ) then
     write(*,*) ' *** Error in dsyev !! **** '
     write(*,*) ' '
     write(*,*) ' *** The ',-info,'-th argument had an illegal value'
  else
     write(*,*) ' *** Error in dsyev !! **** '
     write(*,*) ' '
     write(*,*) ' *** The algorithm failed to converge;', info
     write(*,*) ' *** off-diagonal elements of an intermediate tridiagonal'
     write(*,*) ' *** form did not converge to zero.'
  endif
  
  lz = nlmax - n
  allocate(psi0(nlmax))
  allocate(alphanew(lz))
  allocate(betanew(lz))
  
  psi0 = 0.0d0
  do i = 1, n
     psi0(i) = ulanc(i,1)
  enddo
  
  ! build the full Hamiltonian matrix

  ulanc = 0.0d0
  do i = 1, nlmax-1
     ulanc(i, i) = alpha0(i)
     ulanc(i, i+1) = beta0(i)
     ulanc(i+1, i) = beta0(i)
  enddo
  ulanc(nlmax,nlmax) = alpha0(nlmax)
  call lanczos( ulanc, nlmax, psi0, lz, alphanew, betanew, lunit, n )

  ! now alphanew and betanew contain the Lanczos matrix in the desired representation

  ebarold = 1.d+10
  elower = 0.0d0
  elowernew = -1.0d+2

  ! **************
  ! this is for testing purpose: perform diagonalization of the larger matrix available
  ! recompute eref with the nwe set of parameters
  
  ulanc = 0.0d0
  do i = 1, lz
     ulanc(i, i) = alphanew(i)
     if(i+1.le.lz) ulanc(i, i+1) = betanew(i)
  enddo
  if(idiag.eq.1) then
     call dsyev('V', 'U', lz, ulanc, nlmax, lambda(1,lz), work, lwork, info)
  elseif(idiag.eq.2)then
     call dsyevd('V', 'U', lz, ulanc, nlmax, lambda(1,lz), work, lwork, iwork, liwork,  info)
  endif
  write(*,*) ''
  write(*,*) ' ***** Reference results ***** ' 
  write(*,*) ''
  if(info == 0) then
     write(*,*) ''
     write(*,*) ' SubL = ', lz-1, ' E(i) = ', (lambda(i,lz), i = 1,min(6,lz))
     write(*,*) ''
  elseif(info < 0 ) then
     write(*,*) ' *** Error in dsyev !! **** '
     write(*,*) ' '
     write(*,*) ' *** The ',-info,'-th argument had an illegal value'
  else
     write(*,*) ' *** Error in dsyev !! **** '
     write(*,*) ' '
     write(*,*) ' *** The algorithm failed to converge;', info
     write(*,*) ' *** off-diagonal elements of an intermediate tridiagonal'
     write(*,*) ' *** form did not converge to zero.'
  endif

  eref = lambda(1,lz)   ! this is the best estimate of the ground-state energy we have

  ! compute Q0 using high-precision arithmetics
  
  call FMDP2M( ulanc(1,1), uMP) 
  qboundMP = (oneMP-uMP**2)/uMP**2                  ! this is the exact Q_0 in quadruple precision
  factor = 1.00d0     
  call FMDP2M( factor, factorMP ) 
  qboundMP = qboundMP * factorMP                    ! test check senssitivity!!!
  call FMM2DP(qboundMP,qex)                         ! this is the exact Q_0 in double precision
  if(newbound.and.test) usave(:) = ulanc(:,1)       ! when test=.true. save the exact ground-state vector

  ! this is the dp value of the exact Ebar1
  ebar1ex  = ( alphanew(1)-eref ) / ( 1.d0-ulanc(1,1)**2 ) + eref
  write(*,*) ' Exact Energy, Ebar1, Q0 = ', eref, ebar1ex, qex
  write(*,*) ''
  write(*,*) ' ******************** ' 
  write(*,*) ''
 
  ! compute spectrum of residual Hamiltonian for L=0
  vlanc = 0.0d0
  do i = 2, lz
     vlanc(i, i) = alphanew(i)
     if(i+1.le.lz) vlanc(i, i+1) = betanew(i)
  enddo
  if(idiag.eq.1) then
     call dsyev('N', 'U', lz, vlanc, nlmax, lambdarest, work, lwork, info)
  elseif(idiag.eq.2)then
     call dsyevd('N', 'U', lz, vlanc, nlmax, lambdarest, work, lwork, iwork, liwork,  info)
  endif
  if(info < 0 ) then
     write(*,*) ' *** Error in dsyev !! **** '
     write(*,*) ' '
     write(*,*) ' *** The ',-info,'-th argument had an illegal value'
  elseif( info > 0 ) then
     write(*,*) ' *** Error in dsyev !! **** '
     write(*,*) ' '
     write(*,*) ' *** The algorithm failed to converge;', info
     write(*,*) ' *** off-diagonal elements of an intermediate tridiagonal'
     write(*,*) ' *** form did not converge to zero.'
  endif

  if(.not.newbound)then
     ! ****************
     ! set flag for Q bound
     qflag=.false.
     write(*,*) ' We are almost ready to compute lower bound estimates..'
     write(*,*) '   Select Q-bound estimate' 
     write(*,*) '   Type <bare|improved> for bare or improved overlap  '
     read(*,*) answer
     if (trim(answer)=="bare" )  qflag = .true.
     
     ! if qfixed = .true.  uses energy-independent qbound  
     !           = .false. uses energy-dependent   qbound
     !     test  = .true. for using exact Q0 (enforces qfixed = .true. )
     !      scf  = .true. instruct the code to reach self-consistency
     qfixed=.false.
     write(*,*) '   Select energy dependent Q-bound ' 
     write(*,*) '   Type <T|F> for energy dependent or independent Q bound'
     read(*,*) qfixed
  endif
  write(*,*) '   Ask for self-consistency ' 
  write(*,*) '   Type <T|F> for self-consistent or bare estimate'
  read(*,*) scf
  write(*,*) '   Ask for tight self-consistency ' 
  write(*,*) '   Type <T|F> for tight self-consistent estimate'
  read(*,*) tscf
  if(tscf)then
     open(file='errors_tbound-phi0.dat',unit=26)
     open(file='terrors_scf-phi0.dat',unit=28)
     write(*,*) ' Opening file errors_tbound-phi0.dat      for writing'    ! these are for the (better) upper bound given by the approximate eigenstate 
     write(*,*) ' Opening file terrors_scf-phi0.dat        for writing'    ! save tbound iterative results for the lower bound
     if(.not.newbound)then
        open(file='qbounds-phi0.dat',unit=25)
        open(file='qbounds_tbound-phi0.dat',unit=27)
        write(*,*) ' Opening file qbounds-phi0.dat            for writing'    ! this contains the bound to the effective Q-ratio used in the calculation
        write(*,*) ' Opening file qbounds_tbound-phi0.dat     for writing'    ! this contains the bound to the effective Q-ratio used in the calculation
     endif
  endif
  write(*,*) ' Current status of flags newbound | qflag | qflag2 | QE | SC =', newbound, qflag, qflag2, qfixed, scf
  write(*,*) ' '

  if( test ) qfixed = .true.
  if( test ) scf = .false.
  if( test ) qflag = .false.
  if (test ) qflag2 = .false.
  if( test) then
     write(*,*) ' '
     write(*,*) ' !!!!! This is a test calculation using exact Q0 !!!!!!'
     write(*,*) ' '
  endif
  write(*,*) ' '
  write(*,*) '   Input Ebar1 estimate (Y/N)? '
  read(*,*) ans
  if ( ans == 'Y' ) then
     write(*,*) '   Default value from L=', n-1,' step =', ebar(n)
     write(*,*) '   Is that OK (Y/N)? '
     read(*,*) ans
     if ( ans == 'N' ) then
        write(*,*) '   Input E1bar..'
        read(*,*) ebar1
     else
        ebar1 = ebar(n)
     endif
     eread=.true.
  endif
  write(*,*)' '

  ! ****************  
  
  write(*,*) ''
  write(*,*) ' ******  Compute lower bounds for the desired target state, L = ', n-1
  write(*,*) ' '
  
  lambda(1,1) = alphanew(1)
  if(newbound.and.test)then
     v = usave(1)
     v=min(abs(v),1.0d0)
     call FMDP2M( v, uMP) 
     qboundMP = (oneMP-uMP**2)/uMP**2                  ! this is the exact Q_0 in quadruple precision
     call FMM2DP( qboundMP, qbound)
     qbound = max(qbound,1.d-16)
     ebar1ex  = lambda(1,n) + (lambda(1,n)-eref)/ qbound
  endif
  ! ************* add this for printing purposes only
  n=1
  ulanc(1,1) = 1.d0
  icyc = 1
  jscf = 0
  if(.not.eread) ebar1 = lambda(1,1) + eps
  write(*,*) ' SubL = ', n-1, ' E(i) = ', (lambda(i,n), i = 1,min(6,n))
  write(*,*) ''
  write(*,*) '                    StE(i) = ', ( abs( beta0(n)*ulanc(n,i) ), i = 1,min(6,n) )
  write(*,*) ' ResL = ', n-1, ' E(i) = ', (lambdarest(i), i = 1,min(6,lz-n))
  write(*,*) ' Lambda_0 (bar)        = ', ebar1
  if(newbound.and.test)write(*,*) ' Lambda_0 (bar) exact  = ', ebar1ex  
  write(*,*) ' E_lb (in)          = ', elower
  write(*,*) ' E_lb (out)         = ', elowernew
  if(scf)then
     write(*,*) ' # of self-cycles   = ', icyc-1
  endif
  if(tscf)then
     write(*,*) ' # of tself-cycles  = ', jscf-1
     write(*,*) ' E_lb (tight)       = ', epsl2(n)
  endif
  write(*,*) ' '
  write(16,*) n-1, ebar1
  write(19,*) n-1, (lambdarest(i),i=1,6)
  write(20,*) n-1, (lambda(i,n), i=1,6)
  write(21,*) n-1, elowernew
  write(22,*) n-1, (lambda(i,n)-abs(beta0(n)*ulanc(n,i)), i=1,6 )
  ! *****************

  do k = 2, min(lz,maxd)
     icyc = 1
     jscf = 0
     ulanc = 0.0d0
     vlanc = 0.0d0
     do i = 1, k
        ulanc(i, i) = alphanew(i)
        if(i+1.le.lz) ulanc(i, i+1) = betanew(i)
     enddo
     do i = k+1, lz
        vlanc(i, i) = alphanew(i)
        if(i+1.le.lz) vlanc(i, i+1) = betanew(i)
     enddo
     ! compute the spectrum of the projected Hamiltonian
     if(idiag.eq.1)then
        call dsyev('V', 'U', k, ulanc, nlmax, lambda(1,k), work, lwork, info)
     elseif(idiag.eq.2)then
        call dsyev('V', 'U', k, ulanc, nlmax, lambda(1,k), work, lwork, iwork, liwork, info)
     endif
     if(info == 0) then
        write(*,*) ' SubL = ', k-1, ' E(i) = ', (lambda(i,k), i = 1,min(6,k))
     elseif(info < 0 ) then
        write(*,*) ' *** Error in dsyev !! **** '
        write(*,*) ' '
        write(*,*) ' *** The ',-info,'-th argument had an illegal value'
     else
        write(*,*) ' *** Error in dsyev !! **** '
        write(*,*) ' '
        write(*,*) ' *** The algorithm failed to converge;', info
        write(*,*) ' *** off-diagonal elements of an intermediate tridiagonal'
        write(*,*) ' *** form did not converge to zero.'
     endif
     ! compute the spectrum of the residual Hamiltonian (notice we keep it of the same max dimension)
     if(idiag.eq.1)then
        call dsyev('V', 'U', lz, vlanc, nlmax, lambdarest, work, lwork, info)
     elseif(idiag.eq.2)then
        call dsyev('V', 'U', lz, vlanc, nlmax, lambdarest, work, lwork, iwork, liwork, info)
     endif
     if(info == 0) then
        write(*,*) ' ResL = ', k-1, ' E(i) = ', (lambdarest(i), i = 1,min(6,lz-k))
     elseif(info < 0 ) then
        write(*,*) ' *** Error in dsyev !! **** '
        write(*,*) ' '
        write(*,*) ' *** The ',-info,'-th argument had an illegal value'
     else
        write(*,*) ' *** Error in dsyev !! **** '
        write(*,*) ' '
        write(*,*) ' *** The algorithm failed to converge;', info
        write(*,*) ' *** off-diagonal elements of an intermediate tridiagonal'
        write(*,*) ' *** form did not converge to zero.'
     endif
     ! compute eigenstate variance, using high precision arithmetics
     do i = 1, ni
        call FMDP2M( ulanc(k,i), uMP )
        call FMDP2M( betanew(k), betaMP)
        sigma_L0MP(i) = abs(betaMP*uMP)
     enddo
     ! when test = .true. compute exact Q_0 for the new bound estimate
     if(newbound.and.test)then
        v = 0.0d0
        do i = 1, lz
           v = v + usave(i)*ulanc(i,1)
        enddo
        v=min(abs(v),1.0d0)
        call FMDP2M( v, uMP) 
        qboundMP = (oneMP-uMP**2)/uMP**2                  ! this is the exact Q_0 in quadruple precision
        call FMM2DP( qboundMP, qbound)
        qbound = max(qbound,1.d-16)
        ebar1ex  = lambda(1,k) + (lambda(1,k)-eref)/ qbound 
     endif
     if( k.lt. lz) then
        ! lower bound estimate for the first excited state: it can be larger than E1!
        ebar2 = lambda(2,k)-abs(betanew(k)*ulanc(k,2))
        if(ebar2.lt.lambda(1,k)) then   !!! To be carefully checked for crossing  !!
           if ( qflag ) then
              ebar1 = lambda(1,1) + eps
           else
              ebar1 = lambda(1,k) + eps
           endif
           ebarold = ebar2
        endif
        if(ebar2.lt.ebarold .and. ebar1.lt.ebarold ) then ! .and. ebar2.gt.lambda(1,k) ) then
           ! take ebar1 infinitesimally larger than E0
           if( qflag ) then
              if(.not.eread) ebar1 = lambda(1,1) + eps
           else
              if(.not.eread) ebar1 = lambda(1,k) + eps
           endif
           ebarold = ebar2
        else
        ! take best Weinstein lower bound estimate for the first excited state  
           if( ebar2.gt.lambda(1,k) ) ebar1 = max(ebar2, ebar1)
        endif

11      continue
        if(elowernew.gt.elower) elower = elowernew
        if(elower.gt.lambda(1,k)) elower = lambda(1,k) -abs(betanew(k)*ulanc(k,1))
        
        ! compute bound with high precision (see precision set in routine bounds)
        if(newbound)then
           qtmp = .false.
           if(k.gt.5) qtmp = qflag2
           call bounds_new( k, lambda, nlmax, sigma_L0MP, betanew, ebar1, elower, elowernew, & 
                qboundMP, test, qtmp, eref, fmprec )
           icyc = icyc + 1
           if( scf .and. (elowernew-elower.gt.sth) .and. icyc.lt.icycmax ) goto 11
        else
           if( qflag ) then
              if( .not. qfixed ) scf = .false.        ! nothig to make self-consistent if qfixed is false
              call qbounds( k, lambda, nlmax, sigma_L0MP, betanew, ebar1, elower, elowernew, euppernew, & 
                   epsm(k), qboundMP, rad, qfixed, test, eref, fmprec )
           else
              call bounds( k, lambda, nlmax, sigma_L0MP, betanew, ebar1, elower, elowernew, euppernew, & 
                   epsm(k), qboundMP, rad, qfixed, test, eref, fmprec )
           endif
           icyc = icyc + 1
           if(rad.lt.0.d0) then
              write(*,*) ' Warning!! Rho =', rad,' < 0.d0 !! &
                   Resetting lower bound estimate for self-consistency'
              ! this ensures that elower at the next step is Temple's lower bound
              elower = lambda(1,k) - 2.0d0* abs(betanew(k)*ulanc(k,1))**2/( ebar1-lambda(1,k) ) 
              elowernew = elower +  abs(betanew(k)*ulanc(k,1))**2/( ebar1-lambda(1,k) ) 
           endif
           if( scf .and. (elowernew-elower.gt.sth) .and. icyc.lt.icycmax .and. rad.gt.0.d0 ) goto 11
        endif

        elower2 = elowernew

        if(tscf)then
           if(newbound) then
               ! call tbounds_new
           else
              if ( .not.qflag ) then           
21               call tbounds( k, lambda, nlmax, sigma_L0MP, betanew, ebar1, elower2, elowernew2, euppernew2, & 
                      qboundMP2, rad2, qfixed, fmprec, etmp )
                 jscf = jscf + 1
                 if (rad2.lt.0.0d0) write(*,*) ' TBound warning!! Rho =', rad2,' < 0.0d0 !!'
                 if( jscf.le.jscfmax ) then
                    epstight(k, jscf) = elowernew2
                    elower2 = elowernew2
                    ! save first tbound results
                    if(jscf.eq.1) then
                       ! this is for the tbound
                       epsu2(k) = euppernew2
                       epsl2(k) = elowernew2
                       call FMM2DP(qboundMP2,qb2(k))    ! this is an upper bound to the effective Q-ratio used 
                       rho2(k) = rad2
                    endif
                    goto 21
                 endif
              endif
           endif
        endif
        
        write(*,*) '                    StE(i) = ', ( abs( betanew(k)*ulanc(k,i) ), i = 1,min(6,k) )   
        write(*,*) ' E1 (bar)           = ', ebar1
        if(newbound.and.test)write(*,*) ' E1 (bar) exact     = ', ebar1ex
        write(*,*) ' E_lb (in)          = ', elower
        write(*,*) ' E_lb (out)         = ', elowernew
        if(scf)then
           write(*,*) ' # of self-cycles   = ', icyc-1
        endif
        if(tscf)then
           write(*,*) ' # of tself-cycles  = ', jscf-1
           write(*,*) ' E_lb (tight)       = ', epsl2(k)
        endif
        write(*,*) ' '

        if( k.lt.lz) then
           write(16,*) k-1, ebar1
           write(19,*) k-1, (lambdarest(i),i=1,6)
           write(20,*) k-1, (lambda(i,k), i=1,6)
           write(21,*) k-1, elowernew
           write(22,*) k-1, (lambda(i,k)-abs(betanew(k)*ulanc(k,i)), i=1,6 )
           epsu(k) = euppernew
           epsl(k) = elowernew
           call FMM2DP(qboundMP,qb(k))    ! this is an upper bound to the effective Q-ratio used 
           rho(k) = rad                   ! output "rad" in qbounds.dat
        endif
     endif
  enddo

  do k = 2, lz-1
     write(24,*) k-1, eref-epsl(k), lambda(1,k)-eref, 0.5d0*( lambda(1,k)+epsl(k) )-eref    ! errors.dat
  enddo
  close(16)
  close(19)
  close(20)
  close(21)
  close(22)
  close(24)
  if(tscf)then
     do k = 2, lz-1
        write(26,*) k-1, eref-epsl2(k), lambda(1,k)-eref, 0.5d0*( lambda(1,k)+epsl2(k) )-eref  ! errors-TBOUND.dat
        write(28,*) k-1, (eref-epstight(k,jscf), jscf = 1,jscfmax)
     enddo
     close(26)
     close(28)
     if(.not.newbound)then
        do k = 2, lz -1
           write(27,*) k-1, qb2(k), qex, qb2(k)-qex, rho2(k)                                   ! qbounds-TBOUND.dat
        enddo
        close(27)
     endif
  endif
  if(.not.newbound)then
     do k = 2, lz -1
        write(23,*) k-1, eref-epsl(k), epsu(k)-eref, epsm(k)-eref                              ! errors_lanczos.dat
        write(25,*) k-1, qb(k), qex, qb(k)-qex, rho(k)                                         ! qbounds.dat
     enddo
     close(23)
     close(25)
  endif
  
end Program lowerb
