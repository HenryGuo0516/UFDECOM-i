#include "DEFS.h"	
      SUBROUTINE REUV
C     VERSION(1/21/98)
C----------------------------------------------------------------------
C            VERSION PCGM-WD (9/15/97 MODIFIED BY CHEN)
C     MODIFIED BY WUHUI ,FOR THE PURPOSE OF THE MOVING TIDEFLAT NUMERICAL TEST
C----------------------------------------------------------------------
      USE MOD_GLOBAL
      INTEGER I,K
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 1, KB
c      DO I = 1, N_CTRDP1
      DO concurrent (I=1:N_CTRDP1,K=1:KB)
          UNN(I,K) = U(I,K)
          VNN(I,K) = V(I,K)
          U(I,K)=UF(I,K)
          V(I,K)=VF(I,K)
      ENDDO
c      ENDDO
c!$OMP END TARGET      
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO I = 1,N_CTRD_AG
c      DO K = 1,KB
      DO concurrent (I=1:N_CTRD_AG,K=1:KB) local(UIJK,VIJK)
      IF (FSM(I).GE.0.5) THEN
          UIJK = (U(I,K)+U(NIP1(I),K))*0.5
          VIJK = (V(I,K)+V(NJP1(I),K))*0.5
          IF (FSM(NIP1(I)).EQ.0. .AND. UIJK.GT.0) UIJK = 0.
          IF (FSM(NIM1(I)).EQ.0. .AND. UIJK.LT.0) UIJK = 0.
          IF (FSM(NJP1(I)).EQ.0. .AND. VIJK.GT.0) VIJK = 0.
          IF (FSM(NJM1(I)).EQ.0. .AND. VIJK.LT.0) VIJK = 0.
          UR(I,K) = YETA(I)*UIJK/H2(I)-YXI(I)*VIJK/H1(I)        
          VR(I,K) = -XETA(I)*UIJK/H2(I)+XXI(I)*VIJK/H1(I)
      ELSE 
          UR(I,K) = 0.
          VR(I,K) = 0.
      ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET      

      
      RETURN
      END SUBROUTINE REUV
