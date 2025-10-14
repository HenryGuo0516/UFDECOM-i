#include "DEFS.h"
      
* VERSION (09/09/2009)
* This module integrate different advection schemes by Hui Wu
* (hwu@sklec.ecnu.edu.cn) from State Key Lab. of Estuarine and
* Coastal Research, East China Normal University, Shanghai, China. 
* One can choose following solutions as specified in "DEFS.h":
* a) Center scheme
* b) Upwind scheme
* b) TVD (Total Variation Diminishing) scheme, with superbee limiter
* c) TVD scheme, with van Leer monotonic limiter
* d) TVD scheme, with minimod limiter
* e) Kim scheme with 2nd Runge Kutta method
* f) MPDATA
* g) 3rd HSIMT-TVD scheme
*
* The Center scheme is the origional scheme in ECOM-si. The code for 
* MPDATA scheme, also named as smolar-scheme,is cut from POM with
* slightly modification. The codes for TVD family schemes, the Kim
* scheme,and the 3rd HSIMT-TVD scheme are developed by Hui Wu.
* Among these schemes, the 3rd HSIMT-TVD scheme is found the most
* accurate one, especially near the front region.      
*======================================================================

	
	MODULE ADVT

      SAVE
      CONTAINS

*======================================================================
*                     3RD HSIMT-TVD SCHEME (GPU version)
*======================================================================	
	SUBROUTINE ADVT_HSIMT_GPU(F,FF)
	
C     VERSION(03/02/90)
C
C     THIS SUBROUTINE INTEGRATES CONSERVATIVE CONSTITUENT EQUATIONS

      USE MOD_GLOBAL
      USE MOD_WEIR

      IMPLICIT NONE

      REAL,PARAMETER :: EPSON = 0.0001
      
      INTEGER I,J,K,II,JJ,KK,I0,J0,K0,I1,J1,K1,N,ID,JD,IC,JC,IE,JE
      REAL RL,RTR
      REAL SL,BETAL,COF3,SR,BETAR,SS,RD,RU,SD,BETAD,SU,BETAU
c      REAL EPSON
c      REAL*8 SW(KB)
c      REAL SW(KB)
      REAL SW(11)
      REAL HU1,HU2,HV1,HV2
      REAL ETAA
      REAL ALPHAL,ALPHAR,ALPHAD,ALPHAU,COFMPL
      REAL A1,B1
      REAL RKAL,RKAR,RKAD,RKAU
      REAL CFL,XFLWR,XFLWL,XUPR,XUPL,YFLWU,YFLWD,YUPU,YUPD
	REAL XFLW,XFUP,XFLMU,XFUPL,XFLMUS,WFLW,WFUP,WFLMU,WFLWR
	REAL R,RATIO,RAT,FFU,YFLW,YFUP,YFLMU,YFLWR,YFUPR,YFLMUS,FFV
	REAL AAMX1,AAMY1,WFUPR,WFLMUS,FFF,FFF1,FFF2,DUMT1,DUMT,DVMT,DVMT1
	REAL DTN,DTN1
      REAL ISIGN
c      REAL KAX(N_CTRDP1,KB),KAY(N_CTRDP1,KB),KAZ(N_CTRDP1,KB)
      REAL F(N_CTRDP1,KB),FF(N_CTRDP1,KB)
c      REAL AAMX(N_CTRD,KB),AAMY(N_CTRD,KB) !CBR: Never used.
      
      
c      EPSON = 0.0001
       
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c       DO I = 1,N_CTRD
       DO concurrent (I=1:N_CTRD)
         F(I,KB) = F(I,KBM1)
       ENDDO
c!$OMP END TARGET
        
  
c            XFLUX= 0.0
c            YFLUX= 0.0
c            FF= 0.0
c        KAX = 1.
c        KAY = 1.
c        KAZ = 1.
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO I = 1, N_CTRDP1
c      DO K = 1, KB
      DO concurrent (I=1:N_CTRDP1,K=1:KB)
          XFLUX(I,K)=0.
          YFLUX(I,K)=0.
          FF(I,K)=0.
          KAX(I,K)=1.
          KAY(I,K)=1.
          KAZ(I,K)=1.
      ENDDO
c      ENDDO
c!$OMP END TARGET
      

c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 1,KBM1
c        DO I = 1,N_CTRD
      DO concurrent (I=1:N_CTRD,K=1:KBM1) local(CFL) local_init(DTI)
          CFL = 0.5*(U(I,K)+UNN(I,K))*DTI
     *        /(0.5*(H1(I)+H1(NIM1(I)))+1.E-30)
          KAX(I,K) = 1.-ABS(CFL)

          CFL = 0.5*(V(I,K)+VNN(I,K))*DTI
     *        /(0.5*(H2(I)+H2(NJM1(I)))+1.E-30)
          KAY(I,K) = 1.-ABS(CFL)
        ENDDO
c       ENDDO
c!$OMP END TARGET
        
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c        DO K = 2,KBM1
c          DO I = 1,N_CTRD
      DO concurrent (I=1:N_CTRD,K=2:KBM1) local(CFL) local_init(DTI)
            CFL = 0.5*(W(I,K)+WNN(I,K))*DTI
     *          /(DZZ(K-1)*D(I)+1.E-30)
            KAZ(I,K) = 1.-ABS(CFL)
          ENDDO
c        ENDDO
c!$OMP END TARGET

C******* HORIZONTAL ADVECTION *****************************************
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 1,KBM1
c        DO I = 1,N_CTRD
      DO concurrent (I=1:N_CTRD,K=1:KBM1) local(RL,RKAL,RTR,RKAR,SL,
     *A1,B1,BETAL,SR,BETAR,RD,RKAD,RU,RKAU,SD,BETAD,SU,BETAU)
            
            !XFLUX
            IF(ABS(F(I,K)-F(NIM1(I),K)).LE.EPSON) THEN
              RL = 0.
            ELSE
              RL = (F(NIM1(I),K)-F(NIM2(I),K))/(F(I,K)-F(NIM1(I),K))
              RKAL = KAX(NIM1(I),K)/KAX(I,K)
            ENDIF
            IF(ABS(F(I,K)-F(NIM1(I),K)).LE.EPSON) THEN
              RTR = 0.
            ELSE
              RTR = (F(NIP1(I),K)-F(I,K))/(F(I,K)-F(NIM1(I),K))
              RKAR = KAX(NIP1(I),K)/KAX(I,K)
            ENDIF
            
            IF(U(I,K).GE.0.) THEN
              IF(FSM(NIM2(I)).EQ.0.) THEN
                SL = F(NIM1(I),K)
              ELSE
                IF(FSM(I).EQ.0.) THEN
                  SL = F(NIM1(I),K)
                ELSE
                  A1 = 1./4.*KAX(I,K)+1./2.-1./(12.*KAX(I,K))
                  B1 = -1./4.*KAX(I,K)+1./2.+1./(12.*KAX(I,K))
                  BETAL = A1+B1*RL
c  Note: COF3=AMAX1(0.,AMIN1(2.,2*R,BETA))   !MINMOD
                  SL = F(NIM1(I),K)
     * +0.5*AMAX1(0.,AMIN1(2.,2*RL*RKAL,BETAL))
     *     *(F(I,K)-F(NIM1(I),K))*KAX(I,K)
c     * +0.5*COF3(RL*RKAL,BETAL)*(F(I,K)-F(NIM1(I),K))*KAX(I,K)
                ENDIF
              ENDIF
              
              XFLUX(I,K) = XMFLUX(I,K)*SL
              
            ELSE
              
              IF(FSM(NIP1(I)).EQ.0.) THEN
                SR = F(I,K)
              ELSE
                IF(FSM(NIM1(I)).EQ.0.) THEN
                  SR = F(I,K)
                ELSE
                  A1 = 1./4.*KAX(I,K)+1./2.-1./(12.*KAX(I,K))
                  B1 = -1./4.*KAX(I,K)+1./2.+1./(12.*KAX(I,K))
                  BETAR = A1+B1*RTR
c  Note: COF3=AMAX1(0.,AMIN1(2.,2*R,BETA))   !MINMOD
                  SR = F(I,K)
     * -0.5*AMAX1(0.,AMIN1(2.,2*RTR*RKAR,BETAR))
     *     *(F(I,K)-F(NIM1(I),K))*KAX(I,K)
c     * -0.5*COF3(RTR*RKAR,BETAR)*(F(I,K)-F(NIM1(I),K))*KAX(I,K)
                ENDIF
              ENDIF
              
              XFLUX(I,K) = XMFLUX(I,K)*SR
              
            ENDIF
            
            !YFLUX
            IF(ABS(F(I,K)-F(NJM1(I),K)).LE.EPSON) THEN
              RD = 0.
            ELSE
              RD = (F(NJM1(I),K)-F(NJM2(I),K))/(F(I,K)-F(NJM1(I),K))
              RKAD = KAY(NJM1(I),K)/KAY(I,K)
            ENDIF
            IF(ABS(F(I,K)-F(NJM1(I),K)).LE.EPSON) THEN
              RU = 0.
            ELSE
              RU = (F(NJP1(I),K)-F(I,K))/(F(I,K)-F(NJM1(I),K))
              RKAU = KAY(NJP1(I),K)/KAY(I,K)
            ENDIF
            
            IF(V(I,K).GE.0.) THEN
              IF(FSM(NJM2(I)).EQ.0) THEN
                SD = F(NJM1(I),K)
              ELSE
                IF (FSM(I).EQ.0.) THEN
                  SD = F(NJM1(I),K)
                ELSE
                  A1 = 1./4.*KAY(I,K)+1./2.-1./(12.*KAY(I,K))
                  B1 = -1./4.*KAY(I,K)+1./2.+1./(12.*KAY(I,K))
                  BETAD = A1+B1*RD
c  Note: COF3=AMAX1(0.,AMIN1(2.,2*R,BETA))   !MINMOD
                  SD = F(NJM1(I),K)
     * +0.5*AMAX1(0.,AMIN1(2.,2*RD*RKAD,BETAD))
     *     *((F(I,K)-F(NJM1(I),K))*KAY(I,K))
c     * +0.5*COF3(RD*RKAD,BETAD)*((F(I,K)-F(NJM1(I),K))*KAY(I,K))
                ENDIF
              ENDIF
              
              YFLUX(I,K) = YMFLUX(I,K)*SD
              
            ELSE
              
              IF(FSM(NJP1(I)).EQ.0.) THEN
                SU = F(I,K)
              ELSE
                IF(FSM(NJM1(I)).EQ.0.) THEN
                  SU = F(I,K)
                ELSE
                  A1 = 1./4.*KAY(I,K)+1./2.-1./(12.*KAY(I,K))
                  B1 = -1./4.*KAY(I,K)+1./2.+1./(12.*KAY(I,K))
                  BETAU = A1+B1*RU
c  Note: COF3=AMAX1(0.,AMIN1(2.,2*R,BETA))   !MINMOD
                  SU = F(I,K)
     * -0.5*AMAX1(0.,AMIN1(2.,2*RU*RKAU,BETAU))
     *     *((F(I,K)-F(NJM1(I),K))*KAY(I,K))
c     * -0.5*COF3(RU*RKAU,BETAU)*((F(I,K)-F(NJM1(I),K))*KAY(I,K))
                ENDIF
              ENDIF
              
              YFLUX(I,K) = YMFLUX(I,K)*SU
              
            ENDIF
            
            IF(FSMADD(I)*FSMADD(NIM1(I)).EQ.0) XFLUX(I,K) = 0.
            IF(FSMADD(I)*FSMADD(NJM1(I)).EQ.0) YFLUX(I,K) = 0.
            
        ENDDO
c      ENDDO  
c!$OMP END TARGET
	
C******* HORIZONTAL ADVECTION *****************************************


     
C******  ADD DIFFUSIVE FLUXES *****************************************
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c       DO K = 1,KBM1
c         DO I = 1,N_CTRD
      DO concurrent (I=1:N_CTRD,K=1:KBM1) local(AAMX1,AAMY1)
     * local_init(HPRNU)
           AAMX1 = .5*(AAM(I,K)+AAM(NIM1(I),K))
           XFLUX(I,K) = XFLUX(I,K)
     * -AAMX1/HPRNU*(F(I,K)-F(NIM1(I),K))*DU(I)
     * *(0.5*H2(I)+0.5*H2(NIM1(I)))**2/(0.5*DJ(I)+0.5*DJ(NIM1(I)))

           AAMY1 = .5*(AAM(I,K)+AAM(NJM1(I),K))
           YFLUX(I,K) = YFLUX(I,K)
     * -AAMY1/HPRNU*(F(I,K)-F(NJM1(I),K))*DV(I)
     * *(0.5*H1(I)+0.5*H1(NJM1(I)))**2/(0.5*DJ(I)+0.5*DJ(NJM1(I)))
         ENDDO
c       ENDDO
c!$OMP END TARGET


C******* ADD RIVER INFLOWS   *******************************
C
C        WE NEED TO CHECK CAREFULLY IF THE FLUX FORM USED
C        HERE IS CORRECT
C

C
C       AT THE SOUTHERN AND NORTHERN BOUNDARIES
        
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c        DO N = 1, NUMQBC
c	  DO K = 1, KBM1
      DO concurrent (N=1:NUMQBC,K=1:KBM1) local(ID,IC)
        ID = NQD(N)
        IC = NQC(N)
!------------SOUTH-----------
	  IF(IC==NJM1(ID)) THEN
	    IF(V(ID,K) .GE. 0.0) THEN
	     YFLUX(ID,K) = F(IC,K)*YMFLUX(ID,K)
	    ELSE
	     YFLUX(ID,K) = F(ID,K)*YMFLUX(ID,K)
	    ENDIF
	  ENDIF
!------------NORTH-----------
	  IF(IC==NJP1(ID)) THEN
	    IF(V(IC,K) .GT. 0.0) THEN
	     YFLUX(IC,K) = F(ID,K)*YMFLUX(IC,K)
	    ELSE
	     YFLUX(IC,K) = F(IC,K)*YMFLUX(IC,K)
	    ENDIF
	  ENDIF
!------------EAST-----------
	  IF(IC==NIP1(ID)) THEN
	    IF(U(IC,K) .GT. 0.0) THEN
	     XFLUX(IC,K) = F(ID,K)*XMFLUX(IC,K)
	    ELSE
	     XFLUX(IC,K) = F(IC,K)*XMFLUX(IC,K)
	    ENDIF
	  ENDIF
!------------WEST-----------
	  IF(IC==NIM1(ID)) THEN
	    IF(U(ID,K) .GE. 0.0) THEN
	     XFLUX(ID,K) = F(IC,K)*XMFLUX(ID,K)
	    ELSE
	     XFLUX(ID,K) = F(ID,K)*XMFLUX(ID,K)
	    ENDIF
	  ENDIF
	ENDDO
c      ENDDO
c!$OMP END TARGET


!--------------------------------END ZHH-------------------------------



C-------------------------------WUHUI ADD------------------------------
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO N = 1, NUMEBC
c      DO K=1,KBM1
      DO concurrent (N=1:NUMEBC,K=1:KBM1) local(IE,IC)
      IE = NETA(N)
      IC = NCON(N)
#if defined TIDE_FLUX  || defined TIDE_FLATHER
      IF(IE.EQ.NJP1(IC)) THEN !NORTH
          IF (YMFLUX(IE,K).GT.0.0) THEN
              YFLUX(IE,K) = F(IC,K)*YMFLUX(IE,K)
          ELSE
              YFLUX(IE,K) = F(IE,K)*YMFLUX(IE,K)
          ENDIF
      ELSEIF (IE.EQ.NJM1(IC)) THEN !SOUTH
          IF (YMFLUX(IC,K).GT.0.0) THEN
              YFLUX(IC,K) = F(IE,K)*YMFLUX(IC,K)
          ELSE
              YFLUX(IC,K) = F(IC,K)*YMFLUX(IC,K)
          ENDIF
      ELSEIF (IE.EQ.NIP1(IC)) THEN !EAST
          IF (XMFLUX(IE,K).GT.0.0) THEN
              XFLUX(IE,K) = F(IC,K)*XMFLUX(IE,K)
          ELSE
              XFLUX(IE,K) = F(IE,K)*XMFLUX(IE,K)
          ENDIF
      ELSEIF (IE.EQ.NIM1(IC))THEN !WEST
          IF (XMFLUX(IC,K).GT.0.0) THEN
              XFLUX(IC,K) = F(IE,K)*XMFLUX(IC,K)
          ELSE
              XFLUX(IC,K) = F(IC,K)*XMFLUX(IC,K)
          ENDIF
      ENDIF
#elif defined TIDE_EL
      IF (IE==NIP1(IC)) THEN ! EAST SIDE !V2412
          XFLUX(IE,K) = F(IE,K)*XMFLUX(IE,K)
      ELSEIF (IE==NIM1(IC)) THEN ! WEST  SIDE !V2412
          XFLUX(IC,K) = F(IE,K)*XMFLUX(IC,K)
      ELSEIF (IE==NJP1(IC)) THEN ! NORTH  SIDE
          YFLUX(IE,K) = F(IE,K)*YMFLUX(IE,K)
      ELSEIF (IE==NJM1(IC)) THEN ! SOUTH  SIDE
          YFLUX(IC,K) = F(IE,K)*YMFLUX(IC,K)
      ENDIF
#endif
      ENDDO
c      ENDDO
c!$OMP END TARGET


!******************** THIS SUBROUTINE HAD BEING CHANGED ********************

#ifdef WEIR
      CALL ADDWEIR1(F)
#endif


      
C****** VERTICAL ADVECTION ********************************************
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO !PRIVATE(SW)
c      DO I = 1,N_CTRD_AG
      DO concurrent (I=1:N_CTRD_AG) local_init(KBM1) 
     * local(RD,RKAD,A1,B1,BETAD,K,RU,RKAU,BETAU)
        
        IF(FSM(I).EQ.1.) THEN
          
!          SW(1) = F(I,1)
          A(I,1) = F(I,1)
          
          IF (W(I,2).LT.0.) THEN
!            SW(2) = F(I,1)
             A(I,2) = F(I,1)
          ELSE
            IF (ABS(F(I,2)-F(I,1)).LE.EPSON) THEN
              RD = 0.
            ELSE
              RD = (F(I,3)-F(I,2))/(F(I,2)-F(I,1))
              RKAD = KAZ(I,3)/KAZ(I,2)
            ENDIF
            A1 = 1./4.*KAZ(I,2)+1./2.-1./(12.*KAZ(I,2))
            B1 = -1./4*KAZ(I,2)+1./2.+1./(12.*KAZ(I,2))   
            BETAD = A1+B1*RD
c  Note: COF3=AMAX1(0.,AMIN1(2.,2*R,BETA))   !MINMOD
!            SW(2) = F(I,2)
            A(I,2) = F(I,2)
     * +0.5*AMAX1(0.,AMIN1(2.,2*RD*RKAD,BETAD))
     *     *((F(I,1)-F(I,2))*KAZ(I,2))          
c     * +0.5*COF3(RD*RKAD,BETAD)*((F(I,1)-F(I,2))*KAZ(I,2))          
          ENDIF
          
          DO K = 3, KBM1
            IF (W(I,K).GE.0.) THEN
              IF (K.NE.KBM1) THEN
                IF (ABS(F(I,K)-F(I,K-1)).LE.EPSON) THEN
                  RD = 0.
                ELSE
                  RD = (F(I,K+1)-F(I,K))/(F(I,K)-F(I,K-1))
                  RKAD = KAZ(I,K+1)/KAZ(I,K)
                ENDIF
                A1 = 1./4.*KAZ(I,K)+1./2.-1./(12.*KAZ(I,K))
                B1 = -1./4*KAZ(I,K)+1./2.+1./(12.*KAZ(I,K))   
                BETAD = A1+B1*RD
c  Note: COF3=AMAX1(0.,AMIN1(2.,2*R,BETA))   !MINMOD
!                SW(K) = F(I,K)
                A(I,K) = F(I,K)
     * +0.5*AMAX1(0.,AMIN1(2.,2*RD*RKAD,BETAD))
     *     *((F(I,K-1)-F(I,K))*KAZ(I,K))
c     * +0.5*COF3(RD*RKAD,BETAD)*((F(I,K-1)-F(I,K))*KAZ(I,K))
              ELSE
!                SW(K) = F(I,K)
                A(I,K) = F(I,K)
              ENDIF 
          
            ELSE
              IF (ABS(F(I,K)-F(I,K-1)).LE.EPSON) THEN
                RU = 0.
              ELSE
                RU = (F(I,K-1)-F(I,K-2))/(F(I,K)-F(I,K-1))
                RKAU = KAZ(I,K-1)/KAZ(I,K)
              ENDIF
              A1 = 1./4.*KAZ(I,2)+1./2.-1./(12.*KAZ(I,2))
              B1 = -1./4*KAZ(I,2)+1./2.+1./(12.*KAZ(I,2))   
              BETAU = A1+B1*RU
c  Note: COF3=AMAX1(0.,AMIN1(2.,2*R,BETA))   !MINMOD
!              SW(K) = F(I,K-1)
              A(I,K) = F(I,K-1)
     * -0.5*AMAX1(0.,AMIN1(2.,2*RU*RKAU,BETAU))
     *     *((F(I,K-1)-F(I,K))*KAZ(I,K))
c     * -0.5*COF3(RU*RKAU,BETAU)*((F(I,K-1)-F(I,K))*KAZ(I,K))
            ENDIF 
          
          ENDDO
          
!          SW(KB) = F(I,KBM1)
          A(I,KB) = F(I,KBM1)
          
          DO K = 1,KBM1
!            FF(I,K) = DZR(K)*(SW(K)*W(I,K)-SW(K+1)*W(I,K+1))*DJ(I)
            FF(I,K) = DZR(K)*(A(I,K)*W(I,K)-A(I,K+1)*W(I,K+1))*DJ(I)
          ENDDO
          
        ENDIF
        
	ENDDO
c!$OMP END TARGET

       

C****** ADD NET HORIZONTAL FLUXES; THEN STEP FORWARD IN TIME ***********

c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I = 1,N_CTRD_AG
      DO concurrent (I=1:N_CTRD_AG) local_init(KBM1,DTI) 
     * local(DUMT,DUMT1,DVMT,DVMT1,DTN,DTN1,K)
      IF (FSM(I).GT.0.0) THEN
          DUMT  = FSMADD(I)*FSMADD(NIM1(I))
          DUMT1 = FSMADD(I)*FSMADD(NIP1(I))
          DVMT  = FSMADD(I)*FSMADD(NJM1(I))
          DVMT1 = FSMADD(I)*FSMADD(NJP1(I))
      
          DTN = ((H(I)+2*AMAX1(HU(I),-1.*EL(I))*DUMT)/(1.+2*DUMT)
     * +(H(I)+2*AMAX1(HU(NIP1(I)),-1.*EL(I))*DUMT1)/(1.+2*DUMT1)
     * +(H(I)+2*AMAX1(HV(I),-1.*EL(I))*DVMT)/(1.+2*DVMT)
     * +(H(I)+2*AMAX1(HV(NJP1(I)),-1.*EL(I))*DVMT1)/(1.+2*DVMT1))/4.
           
          DTN1 = ((H(I)+2*AMAX1(HU(I),-1.*ELF(I))*DUMT)/(1.+2*DUMT)
     * +(H(I)+2*AMAX1(HU(NIP1(I)),-1.*ELF(I))*DUMT1)/(1.+2*DUMT1)
     * +(H(I)+2*AMAX1(HV(I),-1.*ELF(I))*DVMT)/(1.+2*DVMT)
     * +(H(I)+2*AMAX1(HV(NJP1(I)),-1.*ELF(I))*DVMT1)/(1.+2*DVMT1))/4.

          DTN1=DTN1+ELF(I)
          DTN=(DTN+EL(I))/DTN1
          DTN1=DTI/(DTN1*DJ(I))

          DO K = 1,KBM1
             FF(I,K) = FF(I,K)
     * +XFLUX(NIP1(I),K)-XFLUX(I,K)+YFLUX(NJP1(I),K)-YFLUX(I,K)
             FF(I,K) = F(I,K)*DTN-FF(I,K)*DTN1
          ENDDO
          
       ENDIF
       ENDDO
c!$OMP END TARGET

c       DO K = 1,KBM1
c         DO I = 1,N_CTRD_AG
c           IF (FSM(I).GT.0.0) THEN

! Note 1:  DO K = 1,KBM1 here and move to Note 2
c             FF(I,K) = FF(I,K)
c     * +XFLUX(NIP1(I),K)-XFLUX(I,K)+YFLUX(NJP1(I),K)-YFLUX(I,K)
      
c             DUMT  = FSMADD(I)*FSMADD(NIM1(I))
c             DUMT1 = FSMADD(I)*FSMADD(NIP1(I))
c             DVMT  = FSMADD(I)*FSMADD(NJM1(I))
c             DVMT1 = FSMADD(I)*FSMADD(NJP1(I))
      
c      DTN = ((H(I)+2*AMAX1(HU(I),-1.*EL(I))*DUMT)/(1.+2*DUMT)
c     * +(H(I)+2*AMAX1(HU(NIP1(I)),-1.*EL(I))*DUMT1)/(1.+2*DUMT1)
c     * +(H(I)+2*AMAX1(HV(I),-1.*EL(I))*DVMT)/(1.+2*DVMT)
c     * +(H(I)+2*AMAX1(HV(NJP1(I)),-1.*EL(I))*DVMT1)/(1.+2*DVMT1))/4.
           
c      DTN1 = ((H(I)+2*AMAX1(HU(I),-1.*ELF(I))*DUMT)/(1.+2*DUMT)
c     * +(H(I)+2*AMAX1(HU(NIP1(I)),-1.*ELF(I))*DUMT1)/(1.+2*DUMT1)
c     * +(H(I)+2*AMAX1(HV(I),-1.*ELF(I))*DVMT)/(1.+2*DVMT)
c     * +(H(I)+2*AMAX1(HV(NJP1(I)),-1.*ELF(I))*DVMT1)/(1.+2*DVMT1))/4.
      
! Note 2:  DO K = 1,KBM1 should Better be here after simplified
      ! DTN1=DTN1+ELF(I)
      ! DTN=(DTN+EL(I))/DTN1
      ! DTN1=DTI/(DTN1*DJ(I))
      ! DO K = 1,KBM1
      !   FF(I,K) = F(I,K)*DTN-FF(I,K)*DTN1
      ! ENDDO
c             FF(I,K) = (F(I,K)*(DTN+EL(I))*DJ(I)-DTI*FF(I,K))
c     * /((DTN1+ELF(I))*DJ(I))      

c           ENDIF
c         ENDDO
c       ENDDO
       
        RETURN
        
      END SUBROUTINE ADVT_HSIMT_GPU
     
      
      
      
*======================================================================
*                     3RD HSIMT-TVD SCHEME
*======================================================================	
	SUBROUTINE ADVT_HSIMT(F,FF)
	
C     VERSION(03/02/90)
C
C     THIS SUBROUTINE INTEGRATES CONSERVATIVE CONSTITUENT EQUATIONS

      USE MOD_GLOBAL
      USE MOD_WEIR

      IMPLICIT NONE
      
      
      INTEGER I,J,K,II,JJ,KK,I0,J0,K0,I1,J1,K1,N,ID,JD,IC,JC,IE,JE
      REAL RL,RTR
      REAL SL,BETAL,COF3,SR,BETAR,SS,RD,RU,SD,BETAD,SU,BETAU
      REAL EPSON
c      REAL*8 SW(KB)
      REAL SW(KB)
      REAL HU1,HU2,HV1,HV2
      REAL ETAA
      REAL ALPHAL,ALPHAR,ALPHAD,ALPHAU,COFMPL
      REAL A1,B1
      REAL RKAL,RKAR,RKAD,RKAU
      REAL CFL,XFLWR,XFLWL,XUPR,XUPL,YFLWU,YFLWD,YUPU,YUPD
	REAL XFLW,XFUP,XFLMU,XFUPL,XFLMUS,WFLW,WFUP,WFLMU,WFLWR
	REAL R,RATIO,RAT,FFU,YFLW,YFUP,YFLMU,YFLWR,YFUPR,YFLMUS,FFV
	REAL AAMX1,AAMY1,WFUPR,WFLMUS,FFF,FFF1,FFF2,DUMT1,DUMT,DVMT,DVMT1
	REAL DTN,DTN1
      REAL ISIGN
c      REAL KAX(N_CTRDP1,KB),KAY(N_CTRDP1,KB),KAZ(N_CTRDP1,KB)
      REAL F(N_CTRDP1,KB),FF(N_CTRDP1,KB)
c      REAL AAMX(N_CTRD,KB),AAMY(N_CTRD,KB) !CBR: Never used.
      
      
      EPSON = 0.0001
       
       DO I = 1,N_CTRD
         F(I,KB) = F(I,KBM1)
       ENDDO
        
  
            XFLUX= 0.0
            YFLUX= 0.0
            FF= 0.0
	  
        
        KAX = 1.
        KAY = 1.
        KAZ = 1.

      DO K = 1,KBM1
#ifdef OMP
!$OMP PARALLEL DO PRIVATE(CFL)
#endif
        DO I = 1,N_CTRD
          CFL = 0.5*(U(I,K)+UNN(I,K))*DTI
     *        /(0.5*(H1(I)+H1(NIM1(I)))+1.E-30)
          KAX(I,K) = 1.-ABS(CFL)

          CFL = 0.5*(V(I,K)+VNN(I,K))*DTI
     *        /(0.5*(H2(I)+H2(NJM1(I)))+1.E-30)
          KAY(I,K) = 1.-ABS(CFL)
        ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif
       ENDDO
        
        DO K = 2,KBM1
#ifdef OMP
!$OMP PARALLEL DO PRIVATE(CFL)
#endif
          DO I = 1,N_CTRD
            CFL = 0.5*(W(I,K)+WNN(I,K))*DTI
     *          /(DZZ(K-1)*D(I)+1.E-30)
            KAZ(I,K) = 1.-ABS(CFL)
          ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif
        ENDDO
        

C******* HORIZONTAL ADVECTION *****************************************
      
      DO K = 1,KBM1
#ifdef OMP
!$OMP PARALLEL DO PRIVATE(RL,RKAL,RTR,RKAR,SL,A1,B1,BETAL,SR,BETAR
     * ,RD,RKAD,RU,RKAU,SD,BETAD,SU,BETAU)
#endif
        DO I = 1,N_CTRD
            
            !XFLUX
            IF(ABS(F(I,K)-F(NIM1(I),K)).LE.EPSON) THEN
              RL = 0.
            ELSE
              RL = (F(NIM1(I),K)-F(NIM2(I),K))/(F(I,K)-F(NIM1(I),K))
              RKAL = KAX(NIM1(I),K)/KAX(I,K)
            ENDIF
            IF(ABS(F(I,K)-F(NIM1(I),K)).LE.EPSON) THEN
              RTR = 0.
            ELSE
              RTR = (F(NIP1(I),K)-F(I,K))/(F(I,K)-F(NIM1(I),K))
              RKAR = KAX(NIP1(I),K)/KAX(I,K)
            ENDIF
            
            IF(U(I,K).GE.0.) THEN
              IF(FSM(NIM2(I)).EQ.0.) THEN
                SL = F(NIM1(I),K)
              ELSE
                IF(FSM(I).EQ.0.) THEN
                  SL = F(NIM1(I),K)
                ELSE
                  A1 = 1./4.*KAX(I,K)+1./2.-1./(12.*KAX(I,K))
                  B1 = -1./4.*KAX(I,K)+1./2.+1./(12.*KAX(I,K))
                  BETAL = A1+B1*RL
                  SL = F(NIM1(I),K)
     * +0.5*COF3(RL*RKAL,BETAL)*(F(I,K)-F(NIM1(I),K))*KAX(I,K)
                ENDIF
              ENDIF
              
              XFLUX(I,K) = XMFLUX(I,K)*SL
              
            ELSE
              
              IF(FSM(NIP1(I)).EQ.0.) THEN
                SR = F(I,K)
              ELSE
                IF(FSM(NIM1(I)).EQ.0.) THEN
                  SR = F(I,K)
                ELSE
                  A1 = 1./4.*KAX(I,K)+1./2.-1./(12.*KAX(I,K))
                  B1 = -1./4.*KAX(I,K)+1./2.+1./(12.*KAX(I,K))
                  BETAR = A1+B1*RTR
                  SR = F(I,K)
     * -0.5*COF3(RTR*RKAR,BETAR)*(F(I,K)-F(NIM1(I),K))*KAX(I,K)
                ENDIF
              ENDIF
              
              XFLUX(I,K) = XMFLUX(I,K)*SR
              
            ENDIF
            
            !YFLUX
            IF(ABS(F(I,K)-F(NJM1(I),K)).LE.EPSON) THEN
              RD = 0.
            ELSE
              RD = (F(NJM1(I),K)-F(NJM2(I),K))/(F(I,K)-F(NJM1(I),K))
              RKAD = KAY(NJM1(I),K)/KAY(I,K)
            ENDIF
            IF(ABS(F(I,K)-F(NJM1(I),K)).LE.EPSON) THEN
              RU = 0.
            ELSE
              RU = (F(NJP1(I),K)-F(I,K))/(F(I,K)-F(NJM1(I),K))
              RKAU = KAY(NJP1(I),K)/KAY(I,K)
            ENDIF
            
            IF(V(I,K).GE.0.) THEN
              IF(FSM(NJM2(I)).EQ.0) THEN
                SD = F(NJM1(I),K)
              ELSE
                IF (FSM(I).EQ.0.) THEN
                  SD = F(NJM1(I),K)
                ELSE
                  A1 = 1./4.*KAY(I,K)+1./2.-1./(12.*KAY(I,K))
                  B1 = -1./4.*KAY(I,K)+1./2.+1./(12.*KAY(I,K))
                  BETAD = A1+B1*RD
                  SD = F(NJM1(I),K)
     * +0.5*COF3(RD*RKAD,BETAD)*((F(I,K)-F(NJM1(I),K))*KAY(I,K))
                ENDIF
              ENDIF
              
              YFLUX(I,K) = YMFLUX(I,K)*SD
              
            ELSE
              
              IF(FSM(NJP1(I)).EQ.0.) THEN
                SU = F(I,K)
              ELSE
                IF(FSM(NJM1(I)).EQ.0.) THEN
                  SU = F(I,K)
                ELSE
                  A1 = 1./4.*KAY(I,K)+1./2.-1./(12.*KAY(I,K))
                  B1 = -1./4.*KAY(I,K)+1./2.+1./(12.*KAY(I,K))
                  BETAU = A1+B1*RU
                  SU = F(I,K)
     * -0.5*COF3(RU*RKAU,BETAU)*((F(I,K)-F(NJM1(I),K))*KAY(I,K))
                ENDIF
              ENDIF
              
              YFLUX(I,K) = YMFLUX(I,K)*SU
              
            ENDIF
            
            IF(FSMADD(I)*FSMADD(NIM1(I)).EQ.0) XFLUX(I,K) = 0.
            IF(FSMADD(I)*FSMADD(NJM1(I)).EQ.0) YFLUX(I,K) = 0.
            
        ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif
      ENDDO  
      
	
C******* HORIZONTAL ADVECTION *****************************************


     
C******  ADD DIFFUSIVE FLUXES *****************************************
      
       DO K = 1,KBM1
#ifdef OMP
!$OMP PARALLEL DO PRIVATE(AAMX1,AAMY1)
#endif
         DO I = 1,N_CTRD
           AAMX1 = .5*(AAM(I,K)+AAM(NIM1(I),K))
           XFLUX(I,K) = XFLUX(I,K)
     * -AAMX1/HPRNU*(F(I,K)-F(NIM1(I),K))*DU(I)
     * *(0.5*H2(I)+0.5*H2(NIM1(I)))**2/(0.5*DJ(I)+0.5*DJ(NIM1(I)))

           AAMY1 = .5*(AAM(I,K)+AAM(NJM1(I),K))
           YFLUX(I,K) = YFLUX(I,K)
     * -AAMY1/HPRNU*(F(I,K)-F(NJM1(I),K))*DV(I)
     * *(0.5*H1(I)+0.5*H1(NJM1(I)))**2/(0.5*DJ(I)+0.5*DJ(NJM1(I)))
         ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif
	 ENDDO
	
	 
C******* ADD RIVER INFLOWS   *******************************
C
C        WE NEED TO CHECK CAREFULLY IF THE FLUX FORM USED
C        HERE IS CORRECT
C

C
C       AT THE SOUTHERN AND NORTHERN BOUNDARIES
        
        DO  N = 1, NUMQBC
        ID = NQD(N)
        IC = NQC(N)
!------------SOUTH-----------
	  IF(IC==NJM1(ID)) THEN
	   DO K = 1, KBM1
	    IF(V(ID,K) .GE. 0.0) THEN
	     YFLUX(ID,K) = F(IC,K)*YMFLUX(ID,K)
	    ELSE
	     YFLUX(ID,K) = F(ID,K)*YMFLUX(ID,K)
	    ENDIF
	   ENDDO
	  ENDIF
!------------NORTH-----------
	  IF(IC==NJP1(ID)) THEN
	   DO K = 1, KBM1
	    IF(V(IC,K) .GT. 0.0) THEN
	     YFLUX(IC,K) = F(ID,K)*YMFLUX(IC,K)
	    ELSE
	     YFLUX(IC,K) = F(IC,K)*YMFLUX(IC,K)
	    ENDIF
	   ENDDO
	  ENDIF
!------------EAST-----------
	  IF(IC==NIP1(ID)) THEN
	   DO K = 1, KBM1
	    IF(U(IC,K) .GT. 0.0) THEN
	     XFLUX(IC,K) = F(ID,K)*XMFLUX(IC,K)
	    ELSE
	     XFLUX(IC,K) = F(IC,K)*XMFLUX(IC,K)
	    ENDIF
	   ENDDO
	  ENDIF
!------------WEST-----------
	  IF(IC==NIM1(ID)) THEN
	   DO K = 1, KBM1
	    IF(U(ID,K) .GE. 0.0) THEN
	     XFLUX(ID,K) = F(IC,K)*XMFLUX(ID,K)
	    ELSE
	     XFLUX(ID,K) = F(ID,K)*XMFLUX(ID,K)
	    ENDIF
	   ENDDO
	  ENDIF
      ENDDO


!--------------------------------END ZHH-------------------------------



C-------------------------------WUHUI ADD------------------------------
      
      DO 141 N = 1, NUMEBC
      IE = NETA(N)
      IC = NCON(N)
#if defined TIDE_FLUX  || defined TIDE_FLATHER
      DO 132 K=1,KBM1
      IF(IE.EQ.NJP1(IC)) THEN !NORTH
          IF (YMFLUX(IE,K).GT.0.0) THEN
              YFLUX(IE,K) = F(IC,K)*YMFLUX(IE,K)
          ELSE
              YFLUX(IE,K) = F(IE,K)*YMFLUX(IE,K)
          ENDIF
      ELSEIF (IE.EQ.NJM1(IC)) THEN !SOUTH
          IF (YMFLUX(IC,K).GT.0.0) THEN
              YFLUX(IC,K) = F(IE,K)*YMFLUX(IC,K)
          ELSE
              YFLUX(IC,K) = F(IC,K)*YMFLUX(IC,K)
          ENDIF
      ELSEIF (IE.EQ.NIP1(IC)) THEN !EAST
          IF (XMFLUX(IE,K).GT.0.0) THEN
              XFLUX(IE,K) = F(IC,K)*XMFLUX(IE,K)
          ELSE
              XFLUX(IE,K) = F(IE,K)*XMFLUX(IE,K)
          ENDIF
      ELSEIF (IE.EQ.NIM1(IC))THEN !WEST
          IF (XMFLUX(IC,K).GT.0.0) THEN
              XFLUX(IC,K) = F(IE,K)*XMFLUX(IC,K)
          ELSE
              XFLUX(IC,K) = F(IC,K)*XMFLUX(IC,K)
          ENDIF
      ENDIF
132   CONTINUE
#elif defined TIDE_EL
      DO 131 K = 1, KBM1
      IF (IE==NIP1(IC)) THEN ! EAST SIDE !V2412
          XFLUX(IE,K) = F(IE,K)*XMFLUX(IE,K)
      ELSEIF (IE==NIM1(IC)) THEN ! WEST  SIDE !V2412
          XFLUX(IC,K) = F(IE,K)*XMFLUX(IC,K)
      ELSEIF (IE==NJP1(IC)) THEN ! NORTH  SIDE
          YFLUX(IE,K) = F(IE,K)*YMFLUX(IE,K)
      ELSEIF (IE==NJM1(IC)) THEN ! SOUTH  SIDE
          YFLUX(IC,K) = F(IE,K)*YMFLUX(IC,K)
      ENDIF
131   CONTINUE
#endif
141   CONTINUE
!******************** THIS SUBROUTINE HAD BEING CHANGED ********************

#ifdef WEIR
      CALL ADDWEIR1(F)
#endif
      
      
C****** VERTICAL ADVECTION ********************************************
#ifdef OMP
!$OMP PARALLEL DO PRIVATE(SW,RD,RKAD,A1,B1,BETAD,RU,RKAU,BETAU,k)
#endif
      DO I = 1,N_CTRD_AG
        
        IF(FSM(I).EQ.1.) THEN
          
          SW(1) = F(I,1)
          
          IF (W(I,2).LT.0.) THEN
            SW(2) = F(I,1)
          ELSE
            IF (ABS(F(I,2)-F(I,1)).LE.EPSON) THEN
              RD = 0.
            ELSE
              RD = (F(I,3)-F(I,2))/(F(I,2)-F(I,1))
              RKAD = KAZ(I,3)/KAZ(I,2)
            ENDIF
            A1 = 1./4.*KAZ(I,2)+1./2.-1./(12.*KAZ(I,2))
            B1 = -1./4*KAZ(I,2)+1./2.+1./(12.*KAZ(I,2))   
            BETAD = A1+B1*RD
            SW(2) = F(I,2)
     * +0.5*COF3(RD*RKAD,BETAD)*((F(I,1)-F(I,2))*KAZ(I,2))          
          ENDIF
          
          DO K = 3, KBM1
            IF (W(I,K).GE.0.) THEN
              IF (K.NE.KBM1) THEN
                IF (ABS(F(I,K)-F(I,K-1)).LE.EPSON) THEN
                  RD = 0.
                ELSE
                  RD = (F(I,K+1)-F(I,K))/(F(I,K)-F(I,K-1))
                  RKAD = KAZ(I,K+1)/KAZ(I,K)
                ENDIF
                A1 = 1./4.*KAZ(I,K)+1./2.-1./(12.*KAZ(I,K))
                B1 = -1./4*KAZ(I,K)+1./2.+1./(12.*KAZ(I,K))   
                BETAD = A1+B1*RD
                SW(K) = F(I,K)
     * +0.5*COF3(RD*RKAD,BETAD)*((F(I,K-1)-F(I,K))*KAZ(I,K))
              ELSE
                SW(K) = F(I,K)
              ENDIF 
          
            ELSE
              IF (ABS(F(I,K)-F(I,K-1)).LE.EPSON) THEN
                RU = 0.
              ELSE
                RU = (F(I,K-1)-F(I,K-2))/(F(I,K)-F(I,K-1))
                RKAU = KAZ(I,K-1)/KAZ(I,K)
              ENDIF
              A1 = 1./4.*KAZ(I,2)+1./2.-1./(12.*KAZ(I,2))
              B1 = -1./4*KAZ(I,2)+1./2.+1./(12.*KAZ(I,2))   
              BETAU = A1+B1*RU
              SW(K) = F(I,K-1)
     * -0.5*COF3(RU*RKAU,BETAU)*((F(I,K-1)-F(I,K))*KAZ(I,K))
            ENDIF 
          
          ENDDO
          
          SW(KB) = F(I,KBM1)
          
          DO K = 1,KBM1
            FF(I,K) = DZR(K)*(SW(K)*W(I,K)-SW(K+1)*W(I,K+1))*DJ(I)
          ENDDO
          
        ENDIF
        
	ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif	

C****** ADD NET HORIZONTAL FLUXES; THEN STEP FORWARD IN TIME ***********
  
       DO K = 1,KBM1
#ifdef OMP
!$OMP PARALLEL DO PRIVATE(DUMT,DUMT1,DVMT,DVMT1,DTN,DTN1)
#endif  
         DO I = 1,N_CTRD_AG
           IF (FSM(I).GT.0.0) THEN

! Note 1:  DO K = 1,KBM1 here and move to Note 2
             FF(I,K) = FF(I,K)
     * +XFLUX(NIP1(I),K)-XFLUX(I,K)+YFLUX(NJP1(I),K)-YFLUX(I,K)
      
             DUMT  = FSMADD(I)*FSMADD(NIM1(I))
             DUMT1 = FSMADD(I)*FSMADD(NIP1(I))
             DVMT  = FSMADD(I)*FSMADD(NJM1(I))
             DVMT1 = FSMADD(I)*FSMADD(NJP1(I))
      
      DTN = ((H(I)+2*AMAX1(HU(I),-1.*EL(I))*DUMT)/(1.+2*DUMT)
     * +(H(I)+2*AMAX1(HU(NIP1(I)),-1.*EL(I))*DUMT1)/(1.+2*DUMT1)
     * +(H(I)+2*AMAX1(HV(I),-1.*EL(I))*DVMT)/(1.+2*DVMT)
     * +(H(I)+2*AMAX1(HV(NJP1(I)),-1.*EL(I))*DVMT1)/(1.+2*DVMT1))/4.
           
      DTN1 = ((H(I)+2*AMAX1(HU(I),-1.*ELF(I))*DUMT)/(1.+2*DUMT)
     * +(H(I)+2*AMAX1(HU(NIP1(I)),-1.*ELF(I))*DUMT1)/(1.+2*DUMT1)
     * +(H(I)+2*AMAX1(HV(I),-1.*ELF(I))*DVMT)/(1.+2*DVMT)
     * +(H(I)+2*AMAX1(HV(NJP1(I)),-1.*ELF(I))*DVMT1)/(1.+2*DVMT1))/4.
      
! Note 2:  DO K = 1,KBM1 should Better be here after simplified
      ! DTN1=DTN1+ELF(I)
      ! DTN=(DTN+EL(I))/DTN1
      ! DTN1=DTI/(DTN1*DJ(I))
      ! DO K = 1,KBM1
      !   FF(I,K) = F(I,K)*DTN-FF(I,K)*DTN1
      ! ENDDO
             FF(I,K) = (F(I,K)*(DTN+EL(I))*DJ(I)-DTI*FF(I,K))
     * /((DTN1+ELF(I))*DJ(I))      

           ENDIF
         ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif
       ENDDO
       
        RETURN
        
      END SUBROUTINE ADVT_HSIMT
     
	  
	  
	END MODULE ADVT
	  
	  