#include "DEFS.h"	
      SUBROUTINE VARBD_GPU(INDEX)
C---- VERSION: 23/10/2002,MADE BY WUHUI,FOR THE MOVING BOUNDARY STUDIES
C---------------------------------------------------------------------
      
      USE MOD_GLOBAL
      IMPLICIT NONE
      
      INTEGER INDEX,SLUICEPD !V2410
      REAL HR !V2410
      INTEGER I,J,K,II,JJ,KK,I1,J1,K1,N,IE,JE,IC,JC
	REAL,PARAMETER :: DD1 = 0.05
	REAL HUU,HVV,HUU1,HVV1
      
      SELECT CASE (INDEX)
      
      CASE (1)
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I = 1,N_CTRD
      DO concurrent (I=1:N_CTRD)
c          HU(I) = 0.5*(H(I)+H(NIM1(I))) !CBR: moved to CASE(11)
c          HV(I) = 0.5*(H(I)+H(NJM1(I))) !CBR: moved to CASE(11)
          D(I) = H(I)+ELF(I)
          FSM11(I) = FSM(I)
      ENDDO
c!$OMP END TARGET


c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I = 1,N_CTRDP1
      DO concurrent (I=1:N_CTRDP1)
          FSM(I)=0.0
          DUM(I)=0.0
          DVM(I)=0.0
      ENDDO
c!$OMP END TARGET
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I = 1,N_CTRD
      DO concurrent (I=1:N_CTRD) local(HUU,HVV,HUU1,HVV1)
     * local_init(DMIN)
          IF(D(I).LT.DMIN) THEN
          IF(H(I).GE.-20) THEN
              ELF(I) = AMAX1(ELF(I),-H(I)+0.8*DMIN)
              D(I) = H(I)+ELF(I)
          ENDIF
          ENDIF
          HUU = 0.5*(H(I)+H(NIM1(I)))
          HVV = 0.5*(H(I)+H(NJM1(I)))
          HUU1 = AMIN1(H(I),H(NIM1(I)))
          HVV1 = AMIN1(H(I),H(NJM1(I)))
          DU(I) = HUU+0.5*(ELF(I)+ELF(NIM1(I)))
          DV(I) = HVV+0.5*(ELF(I)+ELF(NJM1(I)))
          IF(HUU1.LT.-20.) DU(I) = 0.0
          IF(HVV1.LT.-20.) DV(I) = 0.0
      ENDDO
c!$OMP END TARGET

      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I = 1,N_CTRD
      DO concurrent (I=1:N_CTRD) local_init(DMIN)
      IF(DU(I).GE.DMIN) THEN
          IF(D(I).GE.DMIN.AND.D(NIM1(I)).GE.DMIN) THEN
              DUM(I) = 1.0
          ENDIF
          IF(D(I).GE.DMIN.AND.D(NIM1(I)).LT.DMIN.AND.
     *    (ELF(I)-ELF(NIM1(I))).GE.DD1) THEN
              DUM(I) = 1.0
          ENDIF
          IF(D(I).LT.DMIN.AND.D(NIM1(I)).GE.DMIN.AND.
     *    (ELF(NIM1(I))-ELF(I)).GE.DD1) THEN
              DUM(I) = 1.0
          ENDIF
      ENDIF

      IF(DV(I).GE.DMIN) THEN
          IF(D(I).GE.DMIN.AND.D(NJM1(I)).GE.DMIN) THEN
              DVM(I) = 1.0
          ENDIF
          IF(D(I).GE.DMIN.AND.D(NJM1(I)).LT.DMIN.AND.
     *    (ELF(I)-ELF(NJM1(I))).GE.DD1) THEN
              DVM(I) = 1.0
          ENDIF
          IF(D(I).LT.DMIN.AND.D(NJM1(I)).GE.DMIN.AND.
     *    (ELF(NJM1(I))-ELF(I)).GE.DD1) THEN
              DVM(I) = 1.0
          ENDIF
      ENDIF 
      IF(D(I).GE.DMIN) FSM(I) = 1.0
      ENDDO
c!$OMP END TARGET

      
      IF (NSTEP .GE. 0) THEN
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO N = 1,NUMEBC
      DO concurrent (N=1:NUMEBC) local(IE,IC,I1)
          IE = NETA(N)
          IC = NCON(N)
          IF(IE.EQ.NJP1(IC)) THEN !UP_OBC
              I1=IE
              DO WHILE(I1.NE.N_CTRDP1)
                  FSM(I1) = 0
c                  FSMADD(I1) = 0 !CBR: moved to CASE(11)
                  DUM(I1) = 0
                  DUM(NIP1(I1)) = 0
                  I1 = NJP1(I1)
              ENDDO
              I1 = NJP1(IE)
              DO WHILE(I1.NE.N_CTRDP1)
                  DVM(I1) = 0
                  I1 = NJP1(I1)
              ENDDO
          ENDIF
          IF(IE.EQ.NJM1(IC)) THEN
              I1=IE
              DO WHILE(I1.NE.N_CTRDP1)
                  FSM(I1) = 0
c                  FSMADD(I1) = 0 !CBR: moved to CASE(11)
                  DUM(I1) = 0
                  DUM(NIP1(I1)) = 0
                  DVM(I1) = 0.
                  I1 = NJM1(I1)
              ENDDO
          ENDIF
          IF(IE.EQ.NIP1(IC)) THEN
              I1=IE
              DO WHILE(I1.NE.N_CTRDP1)
                  FSM(I1) = 0
c                  FSMADD(I1) = 0 !CBR: moved to CASE(11)
                  DVM(I1) = 0
                  DVM(NJP1(I1)) = 0
                  I1 = NIP1(I1)
              ENDDO
              I1 = NIP1(IE)
              DO WHILE(I1.NE.N_CTRDP1)
                  DUM(I1) = 0.
                  I1 = NIP1(I1)
              ENDDO
          ENDIF
          IF(IE.EQ.NIM1(IC)) THEN
              I1=IE
              DO WHILE(I1.NE.N_CTRDP1)
                  FSM(I1) = 0.0
c                  FSMADD(I1) = 0.0 !CBR: moved to CASE(11)
                  DVM(I1) = 0.0
                  DVM(NJP1(I1)) = 0
                  DUM(I1) = 0.0
                  I1 = NIM1(I1)
              ENDDO
          ENDIF
#if defined TIDE_FLATHER
      FSM(IE) = 1
c      FSMADD(IE) = 1 !CBR: moved to CASE(11)
#endif
      ENDDO
c!$OMP END TARGET
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP PARALLEL DO
c      DO N = 1,NUMEBC
      DO concurrent (N=1:NUMEBC) local(II)
          II = NETA(N)
          FSM11(II) = 1.0
      ENDDO
c!$OMP END TARGET
      ENDIF
!======================================================================!V2410
! SLUICE  
#ifdef MODULE_SLUICE 
c --- to be developed based on CPU version VARBD.
#endif
! END: SLUICE   
!==============================================================================!V2410
	RETURN
!============================================================================== 
      CASE (2)
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 1,KBM1
c      DO I = 1,N_CTRD
      DO concurrent (I = 1:N_CTRD,K = 1:KBM1)
          IF(DUM(I).EQ.0.) THEN
              U(I,K) = 0.
          ENDIF
          IF(DVM(I).EQ.0.) THEN
              V(I,K) = 0.
          ENDIF
          IF(FSM(I).EQ.0.) THEN
              W(I,K) = 0.
          ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET
	RETURN

      CASE (11)
c ----------- CBR: EBC FSMADD moved from CASE(1):
      DO N = 1,NUMEBC
          IE = NETA(N)
          IC = NCON(N)
          IF(IE.EQ.NJP1(IC)) THEN !UP_OBC
              I1=IE
              DO WHILE(I1.NE.N_CTRDP1)
                  FSMADD(I1) = 0
                  I1 = NJP1(I1)
              ENDDO
          ENDIF
          IF(IE.EQ.NJM1(IC)) THEN
              I1=IE
              DO WHILE(I1.NE.N_CTRDP1)
                  FSMADD(I1) = 0
                  I1 = NJM1(I1)
              ENDDO
          ENDIF
          IF(IE.EQ.NIP1(IC)) THEN
              I1=IE
              DO WHILE(I1.NE.N_CTRDP1)
                  FSMADD(I1) = 0
                  I1 = NIP1(I1)
              ENDDO
          ENDIF
          IF(IE.EQ.NIM1(IC)) THEN
              I1=IE
              DO WHILE(I1.NE.N_CTRDP1)
                  FSMADD(I1) = 0.0
                  I1 = NIM1(I1)
              ENDDO
          ENDIF
#if defined TIDE_FLATHER
          FSMADD(IE) = 1
#endif
      ENDDO
c ----------- CBR: EBC FSMADD moved from CASE(1).
c ----------- CBR: HU, HV moved from CASE(1):
      DO I = 1,N_CTRD
          HU(I) = 0.5*(H(I)+H(NIM1(I)))
          HV(I) = 0.5*(H(I)+H(NJM1(I)))
      ENDDO
c ----------- CBR: HU, HV moved from CASE(1).
      
      
      END SELECT

	END SUBROUTINE VARBD_GPU
