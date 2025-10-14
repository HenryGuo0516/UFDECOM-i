#INCLUDE "DEFS.h"
      
      SUBROUTINE FLUX_BALANCE
      USE MOD_GLOBAL
      IMPLICIT NONE
      REAL UC1,VC1
      INTEGER I,K,II
      
      XMFLUX=0.
      YMFLUX=0.

      DO K = 1,KBM1
#ifdef OMP
!$OMP PARALLEL DO PRIVATE(UC1,VC1)
#endif
      DO I = 1,N_CTRD
      IF (DUM(I).NE.0.0) THEN
          UC1 = 0.5*UF(I,K)*(H2(I)+H2(NIM1(I)))
     *    -0.25*(H3(I)+H3(NIM1(I)))/(H1(I)+H1(NIM1(I)))
     *    *(VF(I,K)+VF(NIM1(I),K)+VF(NJP1(I),K)+VF(NIM1JP1(I),K)) 
          XMFLUX(I,K) = UC1*DU(I)
      ENDIF
      IF (DVM(I).NE.0.0) THEN
          VC1 = 0.5*VF(I,K)*(H1(I)+H1(NJM1(I)))
     *    -0.25*(H3(I)+H3(NJM1(I)))/(H2(I)+H2(NJM1(I)))
     *    *(UF(I,K)+UF(NIP1(I),K)+UF(NJM1(I),K)+UF(NIP1JM1(I),K))
          YMFLUX(I,K) = VC1*DV(I)
      ENDIF
      ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif
      ENDDO
      DO 11 II=1,NNIAGCW
      I=NIAGCW(II)
      IF (DUM(NIP1(I)).GT.0.0) THEN 
      DO K=1,KBM1
      XMFLUX(NIP1(I),K)=XMFLUX(AIJ(I,6),K)+XMFLUX(AIJ(I,7),K)
     *+XMFLUX(AIJ(I,8),K)
      ENDDO
      ENDIF
11    CONTINUE
      DO 12 II=1,NNIAGCE
      I=NIAGCE(II)
      IF (DUM(I).GT.0.0) THEN 
      DO K=1,KBM1
      XMFLUX(I,K)=XMFLUX(AIJ(I,9),K)+XMFLUX(AIJ(I,10),K)
     *+XMFLUX(AIJ(I,11),K)
      ENDDO
      ENDIF
12    CONTINUE
      DO 13 II=1,NNIAGCN
      I=NIAGCN(II)
      IF (DVM(I).GT.0.0) THEN 
      DO K=1,KBM1
      YMFLUX(I,K)=YMFLUX(AIJ(I,9),K)+YMFLUX(AIJ(I,10),K)
     *+YMFLUX(AIJ(I,11),K)
      ENDDO
      ENDIF
13    CONTINUE 
      DO 14 II=1,NNIAGCS
      I=NIAGCS(II)
      IF (DVM(NJP1(I)).GT.0.0) THEN 
      DO K=1,KBM1
      YMFLUX(NJP1(I),K)=YMFLUX(AIJ(I,6),K)+YMFLUX(AIJ(I,7),K)
     *+YMFLUX(AIJ(I,8),K)
      ENDDO
      ENDIF
14    CONTINUE
      
      RETURN
      END SUBROUTINE FLUX_BALANCE