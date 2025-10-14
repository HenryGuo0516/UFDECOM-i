#include "DEFS.h"
	SUBROUTINE SETDOM
C     INITIAL CONDITION PREPARATION PROGRAM FOR A GENERAL CIRCULATION
C     MODEL IN NON-ORTHOGONAL CURVILINEAR COORDINATES
      USE MOD_GLOBAL 
      REAL, PARAMETER :: PI = 3.1415926/180.
      REAL X0,Y0,Z0,TT,R1,R2,D1,D2,C1,C2,DD1,CC1 ! V2410
      REAL DHMIN ! V2410

      ISTART = 1
      IEND = 1
      NSTEP = 0
      WRITE (IUPRT,5000)
      READ (IUGRD,*)
      DO I=1,N_CTRDP1
          READ (IUGRD,*) II,XR(II),YR(II),LON(II),LAT(II)
     *    ,HTEMP(II),ANG(II),H1(II),H2(II),H3(II),DJ(II)
     *    ,XXI(II),XETA(II),YXI(II),YETA(II),XU(II),YU(II),XV(II),YV(II)
     *    ,LONU(II),LATU(II),LONV(II),LATV(II)
      ENDDO
      DO I=1,N_CTRD
          READ (IUGRD,*) II,NIM1(II),NIM2(II),NIM3(II),NIM4(II),XX,XX
     *    ,NIP1(II),NIP2(II),NIP3(II),NIP4(II),XX,XX
     *    ,NJM1(II),NJM2(II),NJM3(II),NJM4(II),XX,XX
     *    ,NJP1(II),NJP2(II),NJP3(II),NJP4(II),XX,XX
     *    ,NIM1JP1(II),NIM1JM1(II),NIP1JP1(II),NIP1JM1(II)
     *    ,NIM1JM2(II),NIM2JM1(II)
      ENDDO
      CLOSE (IUGRD)
      ALLOCATE (XNODE(4,N_CTRD))      ;XNODE = 0.0
      ALLOCATE (YNODE(4,N_CTRD))      ;YNODE = 0.0
      NB = 4
      IF (COORD.EQ.'XY') THEN
      FN = "./"//TRIM(IN_DIRE)//TRIM(XG)//"node"//TRIM(XG)
     *//"nodec_ag.txt"
      OPEN (IUNAG,FILE=FN,STATUS='OLD')
c      DO WHILE (.NOT.EOF(IUNAG))
      DO WHILE (1) !(.NOT.EOF(IUNAG))
       READ (IUNAG,*,end=1001) II,XNODE(1,II),YNODE(1,II),XNODE(2,II),
     *    YNODE(2,II),XNODE(3,II),YNODE(3,II),XNODE(4,II),YNODE(4,II)
      ENDDO
1001  continue      
      CLOSE (IUNAG)
      FN = "./"//TRIM(IN_DIRE)//TRIM(XG)//"node"//TRIM(XG)
     *//"nodec_evg.txt"
      OPEN (IUNEG,FILE=FN,STATUS='OLD')
c      DO WHILE (.NOT.EOF(IUNEG))
      DO WHILE (1) !(.NOT.EOF(IUNEG))
       READ (IUNEG,*,end=1002) II,XNODE(1,II),YNODE(1,II),XNODE(2,II),
     *    YNODE(2,II),XNODE(3,II),YNODE(3,II),XNODE(4,II),YNODE(4,II)
      ENDDO
1002  continue      
      CLOSE (IUNEG)
      FN = "./"//TRIM(IN_DIRE)//TRIM(XG)//"node"//TRIM(XG)
     *//"nodec_ivg.txt"
      OPEN (IUNIG,FILE=FN,STATUS='OLD')
c      DO WHILE (.NOT.EOF(IUNIG))
      DO WHILE (1) !(.NOT.EOF(IUNIG))
       READ (IUNIG,*,end=1003) II,XNODE(1,II),YNODE(1,II),XNODE(2,II),
     *    YNODE(2,II),XNODE(3,II),YNODE(3,II),XNODE(4,II),YNODE(4,II)
      ENDDO
1003  continue      
      CLOSE (IUNIG)
      ELSE
      FN = "./"//TRIM(IN_DIRE)//TRIM(XG)//"node"//TRIM(XG)
     *//"nodec_ag.dat"
      OPEN (IUNAG,FILE=FN,STATUS='OLD')
c      DO WHILE (.NOT.EOF(IUNAG))
      DO WHILE (1) !(.NOT.EOF(IUNIG))
       READ (IUNAG,*,end=1004) II,XNODE(1,II),YNODE(1,II),XNODE(2,II),
     *    YNODE(2,II),XNODE(3,II),YNODE(3,II),XNODE(4,II),YNODE(4,II)
      ENDDO
1004  continue      
      CLOSE (IUNAG)
      FN = "./"//TRIM(IN_DIRE)//TRIM(XG)//"node"//TRIM(XG)
     *//"nodec_evg.dat"
      OPEN (IUNEG,FILE=FN,STATUS='OLD')
c      DO WHILE (.NOT.EOF(IUNEG))
      DO WHILE (1) !(.NOT.EOF(IUNEG))
       READ (IUNEG,*,end=1005) II,XNODE(1,II),YNODE(1,II),XNODE(2,II),
     *    YNODE(2,II),XNODE(3,II),YNODE(3,II),XNODE(4,II),YNODE(4,II)
      ENDDO
1005  continue      
      CLOSE (IUNEG)
      FN = "./"//TRIM(IN_DIRE)//TRIM(XG)//"node"//TRIM(XG)
     *//"nodec_ivg.dat"
      OPEN (IUNIG,FILE=FN,STATUS='OLD')
c      DO WHILE (.NOT.EOF(IUNIG))
      DO WHILE (1) !(.NOT.EOF(IUNIG))
       READ (IUNIG,*,end=1006) II,XNODE(1,II),YNODE(1,II),XNODE(2,II),
     *    YNODE(2,II),XNODE(3,II),YNODE(3,II),XNODE(4,II),YNODE(4,II)
      ENDDO
1006  continue      
      CLOSE (IUNIG)
      ENDIF
! END: GRID BASIC INFO
!----------------------------------------------------------------------
      NUMCBCADJ=0
      IF (LOG_CBCADJ) THEN
c      DO WHILE(.NOT.EOF(IUCBCADJ)) 
      DO WHILE (1) !(.NOT.EOF(IUCBCADJ)) 
          READ(IUCBCADJ,*,end=1007)
          NUMCBCADJ=NUMCBCADJ+1
      ENDDO
1007  continue      
      ALLOCATE(CBCADJ(NUMCBCADJ))
      REWIND(IUCBCADJ)
      DO I=1,NUMCBCADJ
          READ(IUCBCADJ,*) CBCADJ(I)
      ENDDO
      ENDIF
      CLOSE(IUCBCADJ)  
!----------------------------------------------------------------------
! GRID BATHYMETRIC INFO
	READ (IUDEP,*)
	DO I = 1,N_CTRD
	    READ (IUDEP,*) II,H(II)
	    IF (HTEMP(II).LE.-10) H(II) = -99999.
	ENDDO
	CLOSE (IUDEP)     
      IF (LOG_DCHG) THEN
      DO I = 1,10000000
          READ (IUDCHG,*,ERR=1271,END=1272) II,H(II)
1271      CONTINUE
      ENDDO
1272  CONTINUE
	ENDIF
	CLOSE (IUDCHG)
      H(N_CTRDP1) = -99999.
!END: GRID BATHYMETRIC INFO
!----------------------------------------------------------------------
      DO I = 1,N_CTRDP1
          H1(I) = H1(I) + 1.E-10
          H2(I) = H2(I) + 1.E-10
          H3(I) = H3(I) + 1.E-10
          DJ(I) = DJ(I) + 1.E-10
          COR(I) = 2.*7.292E-5*SIN(LAT(I)*2.*3.14159/360.)  
      ENDDO
      DO I = 1,N_CTRD
          ARU(I) = 0.5*(DJ(I)+DJ(NIM1(I)))
          ARV(I) = 0.5*(DJ(I)+DJ(NJM1(I)))
      ENDDO 
#ifdef TIDE_EL
      IF (NUMEBC .NE. 0) THEN
      DO N = 1,NUMEBC
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
!II为沿IC至IE再向外延伸一格
          DJ(II) = DJ(IE)
          H3(II) = H3(IE)
          XXI(II) = XXI(IE)
          XETA(II) = XETA(IE)
          YXI(II) = YXI(IE)
          YETA(II) = YETA(IE)
          ARU(II) = ARU(IE)
          ARV(II) = ARV(IE)
      ENDDO
      ENDIF
#endif    
! END: GRID GEOMETRIC INFO
!----------------------------------------------------------------------
      
!----------------------------------------------------------------------
C     DEFINE MASK FOR FREE SURFACE HEIGHT = FSM
C     DEFINE MASK FOR (U,V) VELOCITY = (DUM,DVM)

      !MARK
      DO I = 1,N_CTRD
          HU(I) = AMIN1(H(I),H(NIM1(I)))
          HV(I) = AMIN1(H(I),H(NJM1(I)))
          IF(HU(I).GT.-10) HU(I) = 0.5*(H(I)+H(NIM1(I)))
          IF(HV(I).GT.-10) HV(I) = 0.5*(H(I)+H(NJM1(I)))
      ENDDO
      HU(N_CTRDP1) = -99999.
      HV(N_CTRDP1) = -99999.
      DO I = 1,N_CTRD
          DU(I) = HU(I)
          DV(I) = HV(I)
          D(I) = H(I)
		IF ((H(I).GE.-20) .AND. (H(I).LT.DMIN)) D(I) = 0.001   !WUHUI
      ENDDO
      DO I = 1,N_CTRD
          FSM(I) = 0.0
          DUM(I) = 0.0
          DVM(I) = 0.0
	    FSMADD(I) = 0.           !WUHUI
          IF (H(I).GE.-20) FSMADD(I) = 1.0         !WUHUI
          IF (H(I).GT.DMIN) FSM(I) = 1.0
          IF (H(I).GT.DMIN) FSM11(I) = 1.0  !WUHUI
          IF (HU(I).GT.DMIN) DUM(I) = 1.0
          IF (HV(I).GT.DMIN) DVM(I) = 1.0
      ENDDO
     !! DO I = 1,N_CTRDP1
     !! IF((FSM(I).LT.0.5) .AND. (FSMADD(I).GT.0.5)) THEN
	    !!IF (H(I) .GE. 0.) THEN
	    !!    EL(I) = 0.
	    !!    ELF(I) = 0.
	    !!ELSE
	    !!    EL(I) = -H(I)+0.001
	    !!    ELF(I) = EL(I)
	    !!ENDIF
     !! ENDIF
     !! ENDDO
      
      DO I = 1,N_CTRDP1
      IF((FSM(I).LT.0.5)) THEN
          IF  (FSMADD(I).GT.0.5) THEN
	    IF (H(I) .GE. 0.) THEN
	        EL(I) = 0.+EL_RISEUP !+0.4  c V2410
	        ELF(I) = 0.+EL_RISEUP !+0.4  c V2410
	    ELSE
	        EL(I) = -H(I)+0.001
	        ELF(I) = EL(I)
          ENDIF
          ENDIF
      ELSE
          EL(I) = 0.+EL_RISEUP !+0.4  c V2410
          ELF(I) = 0.+EL_RISEUP !+0.4  c V2410
      ENDIF
      ENDDO
      
!------------------------------------------------------------------------------

	
!------------------------------------------------------------------------------
! HORIZONTAL MIXING
#ifdef HORZMIX_CLOSURE
      DO 250 K = 1, KBM1
      DO 240 I = 1, N_CTRDP1
      IF (FSMADD(I).GT.0.0) THEN
          IF (D(I).GE.3) THEN
              HORCON1 = HORCON
          ELSE
              AAA = D(I)/3.
              HORCON1 = AAA*HORCON1+(1-AAA)*0.2*HORCON
          ENDIF
      ENDIF
      AAM(I,K) = HORCON* TPS(I) / FMIN   !WUHUI
240   CONTINUE
250   CONTINUE
#else
        AAM = HORCON   !WUHUI
#endif  
      IF (NUMEBC.NE.0) THEN
      ALLOCATE (AAMEBC(NUMEBC,KB))
      DO N = 1,NUMEBC
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
          DO  K = 1,KB
              AAMEBC(N,K) = AAM(II,K)
          ENDDO
        ENDDO
      ENDIF
!END: HORIZONTAL MIXING
!----------------------------------------------------------------------
      DO 310 K = 2, KBM1
      DO 300 I = 1, N_CTRD
      IF (FSM(I).GT.0.0) THEN
          Q2(I,K) = 1.E-12
          Q2L(I,K) = 1.E-12
          KM(I,K) = UMOL
          KQ(I,K) = UMOL
          KH(I,K) = UMOL / VPRNU
      ENDIF
300   CONTINUE
310   CONTINUE
       UMOL1 = UMOL
#ifdef VERTMIX_CONSTANT
      UMOL = 0.0
#endif
#ifdef MOD_ONLINE_VIEW
      IF (COORD.EQ.'XY') THEN
      XMIN=5000000
      XMAX=-1
      YMIN=10000000
      YMAX=-1
      DO I = 1,N_CTRD_AG
          IF (XR(I).LE.XMIN.AND.H(I).GE.-10) XMIN=XR(I)
          IF (XR(I).GE.XMAX.AND.H(I).GE.-10) XMAX=XR(I)
          IF (YR(I).LE.YMIN.AND.H(I).GE.-10) YMIN=YR(I)
          IF (YR(I).GE.YAMX.AND.H(I).GE.-10) YMAX=YR(I)
      ENDDO
      ELSE
      XMIN=190
      XMAX=-190
      YMIN=190
      YMAX=-190
      DO I = 1,N_CTRD_AG
          IF (LON(I).LE.XMIN.AND.H(I).GE.-10) XMIN=LON(I)
          IF (LON(I).GE.XMAX.AND.H(I).GE.-10) XMAX=LON(I)
          IF (LAT(I).LE.YMIN.AND.H(I).GE.-10) YMIN=LAT(I)
          IF (LAT(I).GE.YAMX.AND.H(I).GE.-10) YMAX=LAT(I)
      ENDDO 
      ENDIF  
      R1=XMAX-XMIN
      R2=YMAX-YMIN
      IF (R2.GE.R1) THEN
          XMIN=XMIN-(R2-R1)/2.
          XMAX=XMAX+(R2-R1)/2.
      ELSE
          YMIN=YMIN-(R1-R2)/2.
          YMAX=YMAX+(R1-R2)/2.
      ENDIF
      !XMIN=299000.
      !XMAX=462000.
      !YMIN=3401000.
      !YMAX=3548000.
#endif 

#ifdef MODULE_WAVE
      DO K=1,481
          WAVEFFP(K)=0.2+(K-1)*0.01
      ENDDO
#endif 

#ifdef OMP
!$OMP PARALLEL DO
#endif      
      DO I = 1, N_CTRD
          IF (FSM(NJP1(I)).EQ.1.) THEN
              DY1(I) = 0.5 * (H2(I)+H2(NJP1(I)))
          ELSE
              DY1(I) = H2(I)
          ENDIF
          IF (FSM(NJM1(I)).EQ.1.) THEN
              DY2(I) = 0.5 * (H2(I)+H2(NJM1(I)))
          ELSE
              DY2(I) = H2(I)
          ENDIF
          IF (FSM(NIP1(I)).EQ.1. .AND. FSM(NIP1JP1(I)).EQ.1.) THEN
              DY3(I) = 0.5 * (H2(NIP1(I))+H2(NIP1JP1(I)))
          ELSEIF (FSM(NIP1(I)).EQ.1. .AND. FSM(NIP1JP1(I)).EQ.0.0) THEN
              DY3(I) = H2(NIP1(I))
          ELSEIF (FSM(NIP1(I)).EQ.0.0 .AND. FSM(NIP1JP1(I)).EQ.1.) THEN
              DY3(I) = H2(NIP1JP1(I))
          ELSEIF (FSM(NIP1(I)).EQ.0.0.AND.FSM(NIP1JP1(I)).EQ.0.0) THEN
              DY3(I) = H2(I)
          ENDIF
          IF (FSM(NIP1(I)).EQ.1. .AND. FSM(NIP1JM1(I)).EQ.1.) THEN
              DY4(I) = 0.5 * (H2(NIP1(I))+H2(NIP1JM1(I)))
          ELSEIF (FSM(NIP1(I)).EQ.1. .AND. FSM(NIP1JM1(I)).EQ.0.0) THEN
              DY4(I) = H2(NIP1(I))
          ELSEIF (FSM(NIP1(I)).EQ.0.0 .AND. FSM(NIP1JM1(I)).EQ.1.) THEN
              DY4(I) = H2(NIP1JM1(I))
          ELSEIF (FSM(NIP1(I)).EQ.0.0.AND.FSM(NIP1JM1(I)).EQ.0.0) THEN
              DY4(I) = H2(I)
          ENDIF
          IF (FSM(NJP1(I)).EQ.1. .AND. FSM(NIM1JP1(I)).EQ.1.) THEN
              DX1(I) = 0.5 * (H1(NJP1(I))+H1(NIM1JP1(I)))
          ELSEIF (FSM(NJP1(I)).EQ.1. .AND. FSM(NIM1JP1(I)).EQ.0.0) THEN
              DX1(I) = H1(NJP1(I))
          ELSEIF (FSM(NJP1(I)).EQ.0. .AND. FSM(NIM1JP1(I)).EQ.1.) THEN
              DX1(I) = H1(NIM1JP1(I))
          ELSEIF (FSM(NJP1(I)).EQ.0. .AND. FSM(NIM1JP1(I)).EQ.0.) THEN
              DX1(I) = H1(I)
          ENDIF
          IF (FSM(NIP1JP1(I)).EQ.1. .AND. FSM(NJP1(I)).EQ.1.) THEN
              DX2(I) = 0.5 * (H1(NJP1(I))+H1(NIP1JP1(I)))
          ELSEIF (FSM(NIP1JP1(I)).EQ.1. .AND. FSM(NJP1(I)).EQ.0.0) THEN
              DX2(I) = H1(NIP1JP1(I))
          ELSEIF (FSM(NIP1JP1(I)).EQ.0. .AND. FSM(NJP1(I)).EQ.1.) THEN
              DX2(I) = H1(NJP1(I))
          ELSEIF (FSM(NIP1JP1(I)).EQ.0. .AND. FSM(NJP1(I)).EQ.0.) THEN
              DX2(I) = H1(I)
          ENDIF
          IF (FSM(NIM1(I)).EQ.1.) THEN
              DX3(I) = 0.5 * (H1(NIM1(I))+H1(I))
          ELSE
              DX3(I) = H1(I)
          ENDIF
          IF (FSM(NIP1(I)).EQ.1.) THEN
              DX4(I) = 0.5 * (H1(NIP1(I))+H1(I))
          ELSE
              DX4(I) = H1(I)
          ENDIF
      ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif

c --- Add NAJA:
      if (NUMEBC>0) then
      DO N = 1,NUMEBC !CBR: NETA - NCON - NAJA
          IE = NETA(N)
          IC = NCON(N)
          IF (IE.EQ.NIP1(IC)) THEN
              NAJA(N) = NIM1(IC)
          ELSEIF (IE.EQ.NIM1(IC)) THEN
              NAJA(N) = NIP1(IC)              
          ELSEIF (IE.EQ.NJP1(IC)) THEN
              NAJA(N) = NJM1(IC)
          ELSEIF (IE.EQ.NJM1(IC)) THEN
              NAJA(N) = NJP1(IC)
          ENDIF
      ENDDO
C --- Check EBC Extension:      
          open (1881,file="EBC_check.dat")
          DO N = 1,NUMEBC !CBR: NETA - NCON - NAJA
              write (1881,"(3I15)") NETA(N),NCON(N),NAJA(N)
          ENDDO
          close (1881)
      endif

!----------------------------------------------------------------------
      DO I = 1,N_CTRDP1
      IF (H(I).GT.0.0) THEN
          SUMH(I) = 2.0        ! 每平方米允许侵蚀多高的泥沙，单位M 
      ENDIF
      ENDDO
      
! --- V2410:
#ifdef MODULE_LAG  
      DO N=1,N_LAG
          DHMIN=1000000
          X0=PENDX_LAG(N)
	    Y0=PENDY_LAG(N)  
	    DO I=1,N_CTRD_AG
	        D1=SQRT((X0-XR(I))**2+(Y0-YR(I))**2)
              IF (D1.LE.DHMIN) THEN
                  DHMIN=D1
                  I0=I
              ENDIF
          ENDDO
          IF (FSM(I0).LT.0) THEN
              PRINT*,'PARTICLE ',N,' IS ON THE LAND AT START MOMENT'
              PAUSE
              STOP
          ENDIF
          PI_LAG(N)=I0         
          R1=X0-XR(I0)
          R2=Y0-YR(I0)
          D1=XR(NIP1(I0))-XR(I0)
          D2=YR(NIP1(I0))-YR(I0)
          DD1=(R1*D1+R2*D2)/SQRT(D1**2+D2**2)
          IF (DD1.GT.0) THEN
              PIX_LAG(N)=DD1/SQRT(D1**2+D2**2)
          ELSE
              D1=XR(NIM1(I0))-XR(I0)
              D2=YR(NIM1(I0))-YR(I0)
              DD1=(R1*D1+R2*D2)/SQRT(D1**2+D2**2)
              PIX_LAG(N)=-DD1/SQRT(D1**2+D2**2)
          ENDIF
          D1=XR(NJP1(I0))-XR(I0)
          D2=YR(NJP1(I0))-YR(I0)
          DD1=(R1*D1+R2*D2)/SQRT(D1**2+D2**2)
          IF (DD1.GT.0) THEN
              PIY_LAG(N)=DD1/SQRT(D1**2+D2**2)
          ELSE
              D1=XR(NJM1(I0))-XR(I0)
              D2=YR(NJM1(I0))-YR(I0)
              DD1=(R1*D1+R2*D2)/SQRT(D1**2+D2**2)
              PIY_LAG(N)=-DD1/SQRT(D1**2+D2**2)
          ENDIF
          IF(ABS(PIX_LAG(N)).GT.0.5.OR.ABS(PIY_LAG(N)).GT.0.5) THEN
              PRINT*,PIX_LAG(N),PIY_LAG(N),'WO ZHEN SHI MAN TOU DA HAN'
              PAUSE
              STOP
          ENDIF
      ENDDO

#endif   
!======================================================================
! SLUICE   
#ifdef MODULE_SLUICE 
      DO I=1,N_SLUICE
	    READ (IUSLUICE,*) SLUICE_NUM(I)  
	    DO J=1,SLUICE_NUM(I)
		    READ (IUSLUICE,*) SLUICE_GRD(I,J),SLUICE_FSM(I,J)
	    ENDDO
	    IF (SLUICE_GRD(I,2).EQ.NJP1(SLUICE_GRD(I,1))) THEN
		    SLUICE_TYP(I)='J_DIRECT'
              DO J=2,SLUICE_NUM(I)
                  IF (SLUICE_GRD(I,J).NE.NJP1(SLUICE_GRD(I,J-1))) THEN
                      PRINT*,'WRONG SLUICE SETTING UP!'
                      PAUSE
                  ENDIF
              ENDDO
          ELSEIF (SLUICE_GRD(I,2).EQ.NJM1(SLUICE_GRD(I,1))) THEN
              SLUICE_TYP(I)='J_DIRECT'
              DO J=2,SLUICE_NUM(I)
                  IF (SLUICE_GRD(I,J).NE.NJM1(SLUICE_GRD(I,J-1))) THEN
                      PRINT*,'WRONG SLUICE SETTING UP!'
                      PAUSE
                  ENDIF
              ENDDO
          ELSEIF (SLUICE_GRD(I,2).EQ.NIP1(SLUICE_GRD(I,1))) THEN
              SLUICE_TYP(I)='I_DIRECT'
              DO J=2,SLUICE_NUM(I)
                  IF (SLUICE_GRD(I,J).NE.NIP1(SLUICE_GRD(I,J-1))) THEN
                      PRINT*,'WRONG SLUICE SETTING UP!'
                      PAUSE
                  ENDIF
              ENDDO
          ELSEIF (SLUICE_GRD(I,2).EQ.NIM1(SLUICE_GRD(I,1))) THEN
              SLUICE_TYP(I)='I_DIRECT'
              DO J=2,SLUICE_NUM(I)
                  IF (SLUICE_GRD(I,J).NE.NIM1(SLUICE_GRD(I,J-1))) THEN
                      PRINT*,'WRONG SLUICE SETTING UP!'
                      PAUSE
                  ENDIF
              ENDDO 
	    ELSE
		    PRINT*, 'WRONG WEIR SETTING UP!'
		    PAUSE
		    STOP
          ENDIF
      ENDDO
#endif
! END: SLUICE   
!======================================================================  
! --- V2410.      
      
      RETURN
5000  FORMAT(/'...... MODEL STARTING UP FROM INITAL CONDITIONS .......')
	END SUBROUTINE SETDOM