#include "DEFS.h"
      SUBROUTINE CAL_TMP
      USE MOD_GLOBAL
      USE MOD_HEAT
      USE ADVT
      IMPLICIT NONE
      INTEGER NCON2
      INTEGER I,J,K,N,II
      NCON2 = 1
#if defined HEATFLUX_SPECI  ||  defined HEATFLUX_BULK
!$ACC UPDATE HOST (T(:,1)) !GPU Developer Code
      CALL SURFRAD_UNIFORM(T)
!$ACC UPDATE DEVICE (WTSURF,SWRAD) !GPU Developer Code
#endif
#ifdef HSIMT	  
      CALL ADVT_HSIMT_GPU(T,UF)  !WU HUI
#else   
      PRINT*, 'Advection scheme undefined!'
      PAUSE
      STOP
#endif
      CALL PROFT_GPU(NCON2,UF,WTSURF)
      CALL BCOND(6)
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 1,KB
c      DO I = 1,N_CTRD_AG
      DO CONCURRENT (I=1:N_CTRD_AG,K=1:KB)
      IF (FSM(I).GE.0.5) THEN
          T(I,K) =UF(I,K)
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
          T(II,K) =UF(II,K)
      ENDDO
c      ENDDO
c!$OMP END TARGET
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO N = 1,NUMQBC
c      DO K = 1,KB
      DO CONCURRENT (N=1:NUMQBC,K=1:KB) LOCAL(II)
          II = NQC(N)
          T(II,K) =UF(II,K)
      ENDDO
c      ENDDO        
c!$OMP END TARGET
      RETURN
      END SUBROUTINE CAL_TMP