subroutine lanczos( H, nlmax, psi, lz, alpha, beta, iunit, n )

      implicit none

      !     this routine performs lz Lanczos iterations with the real symmetric hamiltonian matrix H.       
      !           psi(nlmax) : on input is the initial system vector array (UNCHANGED on exit)
      !          H(ld,nlmax) : Hamiltonian matrix, with leading dimension ld 
      !            alpha(lz) : diagonal of the Lanczos matrix
      !           beta(lz-1) : off-diagonal of the Lanczos matrix
      
      !     current revision date 22 07 2019

      integer(4) :: nlmax, ld, lz, iunit, n
      real(8) :: psi(nlmax)
      real(8) :: H(nlmax, nlmax)
      real(8) :: alpha(lz), beta(lz)

      ! temporary arrays  & variables
      integer(4) :: k, i, info, lwork
      real(8) :: workl
      real(8), allocatable :: v1(:), v2(:), v3(:)
      real(8), allocatable :: work(:), eigen(:), ulanc(:,:)
      
      allocate(v1(nlmax))
      allocate(v2(nlmax))
      allocate(v3(nlmax))
      allocate(eigen(lz))
      allocate(ulanc(lz,lz))
      eigen = 0.0d0
      
      !     normalize initial state

      workl = sqrt( dot_product(psi, psi))
      v1 = psi / workl
      
      !     Initialize Lanczos recursion

      k = 1
      v2 = matmul(H,v1)
      alpha(1) = dot_product(v1,v2)

      
      v2 = v2 - alpha(1)*v1
      beta(1) = sqrt(dot_product(v2,v2))
      v2 = v2 / beta(1)
      
      eigen(1) = alpha(1)           
      write(*,*) ' Stp =', k , ' E(i) = ', (eigen(i),i=1,6)
      write(iunit,*) k+n-1, (eigen(i),i=1,6)

      if(lz.lt.2) return
      
      ! set dimensions of work arrays at its optimal for lz
      
      call dsyev('N', 'U', lz, ulanc, lz, eigen, workl, -1, info)
      if(info == 0) then
         lwork = workl
      else
         write(*,*)' *** Check your code/data !! ***'
         write(*,*)' *** Something went wrong in dimensioning !'
      end if
      allocate(work(lwork))
      
      do k = 2, lz
         v3 = matmul(H,v2)
         alpha(k) = dot_product(v2,v3)
         
         ! check eigenvalues of H_k Lanczos matrix
         
         ulanc = 0.0d0
         do i = 1, k-1
            ulanc(i,i) = alpha(i)
            ulanc(i,i+1) = beta(i)
         enddo
         ulanc(k,k) = alpha(k)
         call dsyev('V', 'U', k, ulanc, lz, eigen, work, lwork, info)
         if(info == 0)then
            write(*,*) ' Stp =', k , ' E(i) = ', (eigen(i), i=1,6)
            write(iunit,*) k + n-1, (eigen(i),i=1,6)
            !write(iunit-1,*) k + n-1, ulanc(2,1)/ulanc(1,1)
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

         ! build next Lanczos vector
         
         v3 = v3 - alpha(k)*v2 - beta(k-1)*v1
         beta(k) = sqrt(dot_product(v3,v3))
         v1 = v2
         v2 = v3 / beta(k)
      enddo
      deallocate(v1)
      deallocate(v2)
      deallocate(v3)
      deallocate(eigen)
      deallocate(ulanc)
      return
    end subroutine lanczos
    
