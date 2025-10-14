#include "DEFS.h"	
      SUBROUTINE VARBD(INDEX)
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
#ifdef OMP
!$OMP PARALLEL DO
#endif      
      DO I = 1,N_CTRD
c          HU(I) = 0.5*(H(I)+H(NIM1(I))) !CBR: moved to CASE(11)
c          HV(I) = 0.5*(H(I)+H(NJM1(I))) !CBR: moved to CASE(11)
          D(I) = H(I)+ELF(I)
          FSM11(I) = FSM(I)
      ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif     
      FSM=0.0
      DUM=0.0
      DVM=0.0
#ifdef OMP
!$OMP PARALLEL DO PRIVATE(HUU,HVV,HUU1,HVV1)
#endif       
      DO I = 1,N_CTRD
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
#ifdef OMP
!$OMP END PARALLEL DO
#endif

#ifdef OMP
!$OMP PARALLEL DO
#endif  
      DO I = 1,N_CTRD
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
#ifdef OMP
!$OMP END PARALLEL DO
#endif      

      
      IF (NSTEP .GE. 0) THEN
      DO N = 1,NUMEBC
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
      FSMADD(IE) = 1
#endif
      ENDDO
      DO N = 1,NUMEBC
          II = NETA(N)
          FSM11(II) = 1.0
      ENDDO
      ENDIF
!======================================================================!V2410
! SLUICE   
#ifdef MODULE_SLUICE 
      HR = MOD(THOUR,24.)
      IF(THOUR.GE.0..AND.THOUR.LE.240) GOTO 113
      IF(THOUR.GE.363..AND.THOUR.LE.625) GOTO 113
      IF(THOUR.GE.747..AND.THOUR.LE.961) GOTO 113
      IF(THOUR.GE.1082..AND.THOUR.LE.1321) GOTO 113
      IF(THOUR.GE.1430) GOTO 113
      IF (HR.GE.6..AND.HR.LE.18.) GOTO 113
      
      DO I=1,N_SLUICE
      IF (TRIM(SLUICE_TYP(I)).EQ.'J_DIRECT') THEN
          IF (LOG_SLUICE_ON(I)) THEN
              SLUICEPD=0
          ELSE
              SLUICEPD=1
          ENDIF
          DO J=1,SLUICE_NUM(I)
              I1=SLUICE_GRD(I,J)
              IF (SLUICE_FSM(I,J).EQ.1) THEN
                  IF (FSM(I1).GT.0..AND.FSM(NIM1(I1)).GT.0..AND.
     *            (EL(I1)-EL(NIM1(I1))).gt.0.001.AND.
     *            UF(I1,1).LT.-0.01) THEN
                      LOG_SLUICE_ON(I)=.TRUE.
                      GO TO 111
                  ENDIF
                  IF (FSM(I1).GT.0..AND.FSM(NIM1(I1)).GT.0..AND.
     *            (EL(I1)-EL(NIM1(I1))).GT.-0.003) THEN
                      SLUICEPD=SLUICEPD+1
                  ENDIF
              ENDIF
          ENDDO
          IF (SLUICEPD.EQ.0) LOG_SLUICE_ON(I)=.FALSE.
111       CONTINUE
          DO J=1,SLUICE_NUM(I)
              I1=SLUICE_GRD(I,J)
              IF (SLUICE_FSM(I,J).EQ.0) THEN
	            DUM(I1)=0
              ELSE
                  IF ( LOG_SLUICE_ON(I)) THEN
                      DUM(I1)=0
                  ENDIF
              ENDIF
          ENDDO
      ELSE
          IF (LOG_SLUICE_ON(I)) THEN
              SLUICEPD=0
          ELSE
              SLUICEPD=1
          ENDIF
          DO J=1,SLUICE_NUM(I)
              I1=SLUICE_GRD(I,J)
              IF (SLUICE_FSM(I,J).EQ.1) THEN
                  IF (FSM(I1).GT.0..AND.FSM(NJM1(I1)).GT.0..AND.
     *            (EL(I1)-EL(NJM1(I1))).gt.0.0001.AND.
     *            V(I1,1).LT.-0.001) THEN
                      LOG_SLUICE_ON(I)=.TRUE.
                      GO TO 112
                  ENDIF
                  IF (FSM(I1).GT.0..AND.FSM(NJM1(I1)).GT.0..AND.
     *            (EL(I1)-EL(NJM1(I1))).GT.-0.003) THEN
                      SLUICEPD=SLUICEPD+1
                  ENDIF
              ENDIF
          ENDDO
          IF (SLUICEPD.EQ.0) LOG_SLUICE_ON(I)=.FALSE.
112       CONTINUE
          DO J=1,SLUICE_NUM(I)
              I1=SLUICE_GRD(I,J)
              IF (SLUICE_FSM(I,J).EQ.0) THEN
	            DVM(I1)=0
              ELSE
                  IF (LOG_SLUICE_ON(I)) THEN
	                DVM(I1)=0
                  ENDIF
              ENDIF
          ENDDO
      ENDIF	
      ENDDO
113   CONTINUE
#endif
! END: SLUICE   
!==============================================================================!V2410      
	RETURN
!============================================================================== 
      CASE (2)
#ifdef OMP
!$OMP PARALLEL DO PRIVATE(K)
#endif            
      DO I = 1,N_CTRD
      IF(DUM(I).EQ.0.) THEN
          DO K = 1,KBM1
              U(I,K) = 0.
          ENDDO
      ENDIF
      IF(DVM(I).EQ.0.) THEN
          DO K = 1,KBM1
              V(I,K) = 0.
          ENDDO
      ENDIF
      IF(FSM(I).EQ.0.) THEN
          DO K = 1,KBM1
              W(I,K) = 0.
          ENDDO
      ENDIF
      ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif        
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

	END SUBROUTINE VARBD
