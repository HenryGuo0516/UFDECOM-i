#include "DEFS.h"	  
      SUBROUTINE VERTVL
      USE MOD_GLOBAL 
      
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 1, KB
c      DO I = 1, N_CTRDP1
      DO concurrent (I=1:N_CTRDP1,K=1:KB)
          WNN(I,K)=W(I,K)
          W(I,K)=0.0
      ENDDO
c      ENDDO
c!$OMP END TARGET      
      

c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I = 1,N_CTRD_AG
      DO concurrent (I=1:N_CTRD_AG)
      IF(FSM(I).GT.0.0) THEN
      DO K = 1,KBM1    
          W(I,K+1) = 
     * W(I,K)+DZ(K)*((XMFLUX(NIP1(I),K)-XMFLUX(I,K)
     * +YMFLUX(NJP1(I),K)-YMFLUX(I,K))/DJ(I)+(ELF(I)-EL(I))/DTI)  
      ENDDO 
      ENDIF
      ENDDO
c!$OMP END TARGET      

      RETURN
      END SUBROUTINE VERTVL