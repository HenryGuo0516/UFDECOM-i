#include "solver_comm.h"

!==================================================
! SOR for Ax=b
!
      subroutine sor(A,IA,JA,b,x,itermax,iterr)
        use omp_lib
        use mod_solver

        implicit none
        INTEGER itermax,iterr
        cdouble :: A(nnz), b(N), x(N)
        int    :: IA(N+1), JA(nnz)
        cdouble :: D(N)

        int     :: i,j,k,idx_s,idx_e 
        ! cdouble :: err(N)  ! x^{(k+1)} - x^{(k)}
        cdouble :: err_nrm ! ||x^{(k+1)} - x^{(k)}||^2
        cdouble :: x_nrm   ! ||x^{(k+1)}||^2
        cdouble :: relerr  ! ||x^{(k+1)} - x^{(k)}|| / ||x^{(k+1)}||
        cdouble :: err_tmp 


!if (iter==0) then
!    write(*,'(A,F5.2)') 'omega: ', omega
!    write(*,'(A,I2)') 'kind(A(1)): ', kind(A(1))
!endif

#ifdef OMP
!$OMP PARALLEL PRIVATE(idx_s,idx_e,i,k) shared(IA,JA,A,D,b) 
!$OMP DO
#endif
! get the (inverse) diagonals of A
        do i = 1, N
            idx_s = IA(i)
            idx_e = IA(i+1)-1
            do k = idx_s, idx_e
                if (JA(k)==i) then
                    D(i) = 1.0/A(k)
				    if (A(k)==0) then
				        write(*,*) 'kkkkkk'
				        pause
				    endif
                    exit
                endif
            enddo
        enddo
#ifdef OMP
!$OMP END DO
#endif

    !if ( mod(iter,10000)==0) then
    !if ( omp_get_thread_num()==0) then
    !write(*,'(A,I2)') 'omp_get_num_threads()=', omp_get_num_threads()
    !endif
    !endif

#ifdef OMP
!$OMP DO
#endif
! scale A with D from left
        do i = 1, N
            idx_s = IA(i)
            idx_e = IA(i+1)-1
            do k = idx_s, idx_e
                A(k) = D(i)*A(k)  
            enddo
        enddo
#ifdef OMP
!$OMP END DO
#endif

#ifdef OMP
!$OMP DO
#endif
! scale b with D from left
        do i = 1, N  
            b(i) = D(i)*b(i)   
        enddo
#ifdef OMP
!$OMP END DO
!$OMP END PARALLEL 
#endif
        PRINT*,itermax
! start iteration
        do iterr = 1, itermax
            ! xold = x
            err_nrm = 0.0
            ! x_nrm   = 0.0
            do i = 1, N
                idx_s = IA(i)
                idx_e = IA(i+1)-1
                err_tmp = b(i) 
                do k = idx_s, idx_e
                    j = JA(k)
                    err_tmp = err_tmp - A(k)*x(j) 
                enddo 
                err_nrm = err_nrm + err_tmp*err_tmp
                x(i)    = x(i) + omega*err_tmp
                ! x_nrm   = x_nrm + x(i)*x(i) 
                ! write(*,'(A,F12.8,A,F12.8)') 'err:', err_tmp, ', x(i):', xnew(i)
            enddo
	      ! x = (1-omega)*xold + omega*x
	      ! xold = x - xold
	      ! err_nrm = dot_product(xold,xold)
	        x_nrm = dot_product(x,x)
            relerr = omega*sqrt(err_nrm)/sqrt(x_nrm) 
            ! write(*,'(A,I3,A,E12.5)') 'iter=', iterr,', relerr=', relerr 
            if (relerr<tol)  return 
        enddo

        return
      end 
!==================================================

