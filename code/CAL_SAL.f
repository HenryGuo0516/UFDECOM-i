#include "DEFS.h"
	SUBROUTINE CAL_SAL

	USE MOD_GLOBAL
	USE ADVT
	
	IMPLICIT NONE
	INTEGER NCON2
	INTEGER I,J,K,N,II
      NCON2 = 2
#ifdef HSIMT
      CALL ADVT_HSIMT_GPU(S,VF)  !WU HUI
#else
      PRINT*, 'Advection scheme undefined!'
	PAUSE
	STOP
#endif
      CALL PROFT_GPU(NCON2,VF,WSSURF)

      CALL BCOND(6) 
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 1,KB
c      DO I = 1,N_CTRD_AG
      DO CONCURRENT (I=1:N_CTRD_AG,K=1:KB)
      IF(FSM(I).GE.0.5) THEN
          S(I,K) = VF(I,K)
      ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO N = 1,NUMEBC
c      DO K = 1,KB
      DO CONCURRENT (N=1:NUMEBC,K=1:KB) LOCAL(II)
          II = NETA(N)
          S(II,K) = VF(II,K)
      ENDDO
c      ENDDO
c!$OMP END TARGET
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO N = 1,NUMQBC
c      DO K = 1,KB
      DO CONCURRENT (N=1:NUMQBC,K=1:KB) LOCAL(II)
          II = NQC(N)
          S(II,K) = VF(II,K)
      ENDDO
c      ENDDO
c!$OMP END TARGET

      RETURN  
	END SUBROUTINE CAL_SAL