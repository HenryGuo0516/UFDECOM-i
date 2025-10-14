#include "DEFS.h"	  
      SUBROUTINE SMAG
      
      USE MOD_GLOBAL

C**********************************************************************
C**********************************************************************
C         0.01 < HORCON <0.20
C  HOR. VISC. = CONST.*H1*H2*SQRT((2.*DU/H1)**2+(2.*DV/H2)**2
C               +2.*(DU/H2+DV/H1)**2)
C**********************************************************************
C**********************************************************************
c      AAM=0.
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 1, KB
c      DO I = 1, N_CTRDP1
      DO concurrent (I=1:N_CTRDP1,K=1:KB)
            AAM(I,K)=0.0
      ENDDO
c      ENDDO
c!$OMP END TARGET
c#ifdef OMP
c!$OMP PARALLEL DO PRIVATE(HORCON1,AAA)
c#endif      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
cc      DO K = 1, KBM1
c      DO I = 1, N_CTRD
      DO concurrent (I=1:N_CTRD,K=1:KBM1) local_init(HORCON)
      IF (FSM(I).GT.0.0) THEN
c      IF (D(I).GE.3) THEN
c          HORCON1 = HORCON
c      ELSE
c          AAA = D(I)/3.
c          HORCON1 = AAA*HORCON+(1-AAA)*0.5*HORCON
c      ENDIF
c      HORCON1 = HORCON
      AAM(I,K) = HORCON*H1(I)*H2(I)*FSM(I)
     * *SQRT(4.0*(((U(NIP1(I),K)-U(I,K))/H1(I))**2
     * +((V(NJP1(I),K)-V(I,K))/H2(I))**2)
     * +0.125*(ABS(U(NJP1(I),K)-U(I,K))/DY1(I)
     * +ABS(U(I,K)-U(NJM1(I),K))/DY2(I)
     * +ABS(U(NIP1JP1(I),K)-U(NIP1(I),K))/DY3(I)
     * +ABS(U(NIP1(I),K)-U(NIP1JM1(I),K))/DY4(I)
     * +ABS(V(NJP1(I),K)-V(NIM1JP1(I),K))/DX1(I)
     * +ABS(V(NIP1JP1(I),K)-V(NJP1(I),K))/DX2(I)
     * +ABS(V(NIM1(I),K)-V(I,K))/DX3(I)
     * +ABS(V(NIP1(I),K)-V(I,K))/DX4(I))**2)
      ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET
c#ifdef OMP
c!$OMP END PARALLEL DO
c#endif

c      AAM(N_CTRDP1,:) = 0.
      IF (NUMQBC.NE.0) THEN
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO N = 1, NUMQBC
c      DO K = 1, KBM1
      DO concurrent (N=1:NUMQBC,K=1:KBM1) local(IC)
            IC = NQC(N)
            AAM(IC,K) = 0.0
      ENDDO
c      ENDDO
c!$OMP END TARGET
      ENDIF
      
      IF (NUMEBC.NE.0) THEN
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO N = 1, NUMEBC
      DO concurrent (N=1:NUMEBC) local(IC,IE,II,K) 
     * local_init(HORCON,KBM1)
          IE = NETA(N)
          IC = NCON(N)
          IF (IE.EQ.NIP1(IC)) THEN
              II = NIP1(IE)
          ELSEIF (IE.EQ.NIM1(IC)) THEN
              II = NIM1(IE)             
          ELSEIF (IE.EQ.NJP1(IC)) THEN
              II = NJP1(IE)
          ELSEIF (IE.EQ.NJM1(IC)) THEN
              II = NJM1(IE)
          ENDIF
          DO K = 1, KBM1
              AAM(IE,K) = AAM(IC,K)
              AAM(II,K) = HORCON
          ENDDO 
      ENDDO
c!$OMP END TARGET
      ENDIF

      RETURN
      END SUBROUTINE SMAG

