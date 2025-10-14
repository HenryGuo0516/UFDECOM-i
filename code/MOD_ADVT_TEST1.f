#include "DEFS.h"
      

	
	MODULE ADVT_TEST1

      SAVE
	CONTAINS

*======================================================================
*                     3RD HSIMT-TVD SCHEME
*======================================================================	
	SUBROUTINE ADVT_HSIMT_TEST1(F,FF)
	

      USE MOD_GLOBAL

      IMPLICIT NONE
      INTEGER I,J,K,II,JJ,KK,I0,J0,K0,I1,J1,K1,N,ID,JD,IC,JC,IE,JE
      INTEGER CNUM
      REAL LOC,X1,X2,X3,Y1,Y2,Y3,AEL,BEL,CEL,F1,F2,F3,F4,CNUM1,CNUM2
      REAL RL,RTR,CNUM3,CNUM4
      REAL SL,BETAL,COF3,SR,BETAR,SS,RD,RU,SD,BETAD,SU,BETAU
      REAL EPSON
      REAL*8 SW(KB)
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
      REAL AAMX(N_CTRD,KB),AAMY(N_CTRD,KB)
      
      
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
      DO K = 1,KBM1
#ifdef OMP
!$OMP PARALLEL DO PRIVATE(RL,RKAL,RTR,RKAR,SL,A1,B1,BETAL,SR,BETAR
     * ,RD,RKAD,RU,RKAU,SD,BETAD,SU,BETAU)
#endif
      DO I = 1,N_CTRD
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
                  SL = F(NIM1(I),K)+0.5*COF3(RL*RKAL,BETAL)
     *            *(F(I,K)-F(NIM1(I),K))*KAX(I,K)
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
                  SR = F(I,K)-0.5*COF3(RTR*RKAR,BETAR)
     *            *(F(I,K)-F(NIM1(I),K))*KAX(I,K)
              ENDIF
          ENDIF
          XFLUX(I,K) = XMFLUX(I,K)*SR
          ENDIF
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
                  SD = F(NJM1(I),K)+0.5*COF3(RD*RKAD,BETAD)
     *            *((F(I,K)-F(NJM1(I),K))*KAY(I,K))
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
                  SU = F(I,K)-0.5*COF3(RU*RKAU,BETAU)
     *            *((F(I,K)-F(NJM1(I),K))*KAY(I,K))
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
#ifdef OMP
!$OMP PARALLEL DO PRIVATE(SW,RD,RKAD,A1,B1,BETAD,RU,RKAU,BETAU,k)
#endif
      DO I = 1,N_CTRD
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
     *        +0.5*COF3(RD*RKAD,BETAD)*((F(I,1)-F(I,2))*KAZ(I,2))
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
     *            +0.5*COF3(RD*RKAD,BETAD)*((F(I,K-1)-F(I,K))*KAZ(I,K))
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
     *            -0.5*COF3(RU*RKAU,BETAU)*((F(I,K-1)-F(I,K))*KAZ(I,K))
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
      DO  N = 1, NUMQBC
      ID = NQD(N)
      IC = NQC(N)
      IF(IC==NJM1(ID)) THEN
	    DO K = 1, KBM1
              IF(V(ID,K) .GE. 0.0) THEN
	            YFLUX(ID,K) = F(IC,K)*YMFLUX(ID,K)
	        ELSE
	            YFLUX(ID,K) = F(ID,K)*YMFLUX(ID,K)
	        ENDIF
	    ENDDO
      ELSEIF(IC==NJP1(ID)) THEN
	    DO K = 1, KBM1
	        IF(V(IC,K) .GT. 0.0) THEN
	            YFLUX(IC,K) = F(ID,K)*YMFLUX(IC,K)
              ELSE
	            YFLUX(IC,K) = F(IC,K)*YMFLUX(IC,K)
	        ENDIF
	    ENDDO
      ELSEIF(IC==NIP1(ID)) THEN
	    DO K = 1, KBM1
	        IF(U(IC,K) .GT. 0.0) THEN
	            XFLUX(IC,K) = F(ID,K)*XMFLUX(IC,K)
	        ELSE
	            XFLUX(IC,K) = F(IC,K)*XMFLUX(IC,K)
	        ENDIF
	    ENDDO
      ELSEIF(IC==NIM1(ID)) THEN
	    DO K = 1, KBM1
	        IF(U(ID,K) .GE. 0.0) THEN
	            XFLUX(ID,K) = F(IC,K)*XMFLUX(ID,K)
	        ELSE
	            XFLUX(ID,K) = F(ID,K)*XMFLUX(ID,K)
	        ENDIF
	    ENDDO
      ENDIF
      ENDDO
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
      IF (IE==NIM1(IC)) THEN ! EAST SIDE
          XFLUX(IE,K) = F(IE,K)*XMFLUX(IE,K)
      ELSEIF (IE==NIP1(IC)) THEN ! WEST  SIDE
          XFLUX(IC,K) = F(IE,K)*XMFLUX(IC,K)
      ELSEIF (IE==NJP1(IC)) THEN ! NORTH  SIDE
          YFLUX(IE,K) = F(IE,K)*YMFLUX(IE,K)
      ELSEIF (IE==NJM1(IC)) THEN ! SOUTH  SIDE
          YFLUX(IC,K) = F(IE,K)*YMFLUX(IC,K)
      ENDIF
131   CONTINUE
#endif
141   CONTINUE
***********************************************************************
***                             MaRin                               ***
***                                                                 ***
      DO II=1,NNIVGFWALL  
      I=NIVGFWALL(II,1)
      IF (FSM(I).GT.0) THEN
      IC=NIVGFWALL(II,2)
      LOC=(YNODE(1,IC)+YNODE(2,IC))/2
      X1=((YNODE(1,NJM1(IC))+YNODE(2,NJM1(IC)))/2-LOC)/H1(I)
      X2=0.0
      X3=((YNODE(1,NJP1(IC))+YNODE(2,NJP1(IC)))/2-LOC)/H1(I)
      DO K=1,KBM1 
      IF (XFLUX(IC,K).NE.0) THEN
      Y1=XMFLUX(NJM1(IC),K)
      Y2=XMFLUX(IC,K)
      Y3=XMFLUX(NJP1(IC),K)
      BEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/(X2-X3)/
     *((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2))
      AEL=(Y1-Y2-BEL*(X1-X2))/(X1**2-X2**2)
      CEL=Y1-AEL*X1**2-BEL*X1
      IF (AEL.EQ.0..AND.BEL.EQ.0..AND.CEL.EQ.0.) CEL=1
      F1=(YNODE(1,IC)-LOC)/H1(I)
      F2=(YNODE(1,I)-LOC)/H1(I)
      F3=(YNODE(2,I)-LOC)/H1(I)
      F4=(YNODE(2,IC)-LOC)/H1(I)
      CNUM1=AEL/3*(F1**3-F4**3)+BEL/2*(F1**2-F4**2)+CEL*(F1-F4)
      CNUM2=AEL/3*(F2**3-F3**3)+BEL/2*(F2**2-F3**2)+CEL*(F2-F3)
      XFLUX(I,K)=CNUM2/CNUM1*XFLUX(IC,K)
      ELSE
      XFLUX(I,K)=0.
      ENDIF    
      ENDDO
      ENDIF
      ENDDO
      
      DO II=1,NNIVGFSALL  
      I=NIVGFSALL(II,1)
      IF (FSM(I).GT.0) THEN
      IC=NIVGFSALL(II,2)
      LOC=(XNODE(1,IC)+XNODE(4,IC))/2
      X1=((XNODE(1,NIM1(IC))+XNODE(4,NIM1(IC)))/2-LOC)/H2(I)
      X2=0.0
      X3=((XNODE(1,NIP1(IC))+XNODE(4,NIP1(IC)))/2-LOC)/H2(I)
      DO K=1,KBM1 
      IF (YFLUX(IC,K).NE.0) THEN
      Y1=YMFLUX(NIM1(IC),K)
      Y2=YMFLUX(IC,K)
      Y3=YMFLUX(NIP1(IC),K)
      BEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/(X2-X3)/
     *((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2))
      AEL=(Y1-Y2-BEL*(X1-X2))/(X1**2-X2**2)
      CEL=Y1-AEL*X1**2-BEL*X1
      IF (AEL.EQ.0..AND.BEL.EQ.0..AND.CEL.EQ.0.) CEL=1
      F1=(XNODE(1,IC)-LOC)/H2(I)
      F2=(XNODE(1,I)-LOC)/H2(I)
      F3=(XNODE(4,I)-LOC)/H2(I)
      F4=(XNODE(4,IC)-LOC)/H2(I)
      CNUM1=AEL/3*(F1**3-F4**3)+BEL/2*(F1**2-F4**2)+CEL*(F1-F4)
      CNUM2=AEL/3*(F2**3-F3**3)+BEL/2*(F2**2-F3**2)+CEL*(F2-F3)
      YFLUX(I,K)=CNUM2/CNUM1*YFLUX(IC,K)
      ELSE
      YFLUX(I,K)=0.
      ENDIF    
      ENDDO
      ENDIF
      ENDDO   
      DO II=1,N_CTRD_IVG
          I=NAG4VGC(II,1)
          DO K=1,KBM1
              XFLUX(I,K)=XFLUX(NAG4VGC(II,2),K)+XFLUX(NAG4VGC(II,3),K)+
     *        XFLUX(NAG4VGC(II,4),K)
              YFLUX(I,K)=YFLUX(NAG4VGC(II,5),K)+YFLUX(NAG4VGC(II,6),K)+
     *        YFLUX(NAG4VGC(II,7),K)
          ENDDO
      ENDDO
      DO II=1,NNIAGCE
      I=NIAGCE(II)
      IF (DUM(I).GT.0.0) THEN 
      DO K=1,KBM1
      XFLUX(I,K)=XFLUX(AIJ(I,9),K)+XFLUX(AIJ(I,10),K)
     *+XFLUX(AIJ(I,11),K)
      ENDDO
      ENDIF
      ENDDO
      DO II=1,NNIAGCN
      I=NIAGCN(II)
      IF (DVM(I).GT.0.0) THEN 
      DO K=1,KBM1
      YMFLUX(I,K)=YMFLUX(AIJ(I,9),K)+YMFLUX(AIJ(I,10),K)
     *+YMFLUX(AIJ(I,11),K)
      ENDDO
      ENDIF
      ENDDO
   
      DO K=1,KBM1
      DO I = 1,N_CTRD
      IF (FSM(I).GT.0.0) THEN    
      FF(I,K) = FF(I,K)
     * +XFLUX(NIP1(I),K)-XFLUX(I,K)+YFLUX(NJP1(I),K)-YFLUX(I,K)
      ENDIF
      ENDDO
      ENDDO
      
      DO II=1,NNIVGFEALL
      I=NIVGFEALL(II,1)
      IF (FSM(I).GT.0) THEN
      IC =NIP1(NIVGFEALL(II,2))
      LOC=(YNODE(1,IC)+YNODE(2,IC))/2
      X1=((YNODE(1,NJM1(IC))+YNODE(2,NJM1(IC)))/2-LOC)/H1(I)
      X2=0.0
      X3=((YNODE(1,NJP1(IC))+YNODE(2,NJP1(IC)))/2-LOC)/H1(I)
      DO K=1,KBM1
      IF (XFLUX(IC,K).NE.0.) THEN
      Y1=XFLUX(NJM1(IC),K)
      Y2=XFLUX(IC,K)
      Y3=XFLUX(NJP1(IC),K)
      BEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/(X2-X3)/
     *((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2))
      AEL=(Y1-Y2-BEL*(X1-X2))/(X1**2-X2**2)
      CEL=Y1-AEL*X1**2-BEL*X1
      IF (AEL.EQ.0..AND.BEL.EQ.0..AND.CEL.EQ.0.) CEL=1
      F1=(YNODE(1,IC)-LOC)/H1(I)
      F2=(YNODE(1,I)-LOC)/H1(I)
      F3=(YNODE(2,I)-LOC)/H1(I)
      F4=(YNODE(2,IC)-LOC)/H1(I)
      CNUM1=AEL/3*(F1**3-F4**3)+BEL/2*(F1**2-F4**2)+CEL*(F1-F4)
      CNUM2=AEL/3*(F2**3-F3**3)+BEL/2*(F2**2-F3**2)+CEL*(F2-F3)
      FF(I,K)=FF(I,K)+CNUM2/CNUM1*XFLUX(IC,K)
      ENDIF
      ENDDO
      ENDIF
      ENDDO
      
      DO II=1,NNIVGFNALL
      I=NIVGFNALL(II,1)
      IF (FSM(I).GT.0) THEN
      IC =NJP1(NIVGFNALL(II,2))
      LOC=(XNODE(1,IC)+XNODE(4,IC))/2
      X1=((XNODE(1,NIM1(IC))+XNODE(4,NIM1(IC)))/2-LOC)/H1(I)
      X2=0.0
      X3=((XNODE(1,NIP1(IC))+XNODE(4,NIP1(IC)))/2-LOC)/H1(I)
      DO K=1,KBM1
      IF (YFLUX(IC,K).NE.0.) THEN
      Y1=YFLUX(NIM1(IC),K)
      Y2=YFLUX(IC,K)
      Y3=YFLUX(NIP1(IC),K)
      BEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/(X2-X3)/
     *((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2))
      AEL=(Y1-Y2-BEL*(X1-X2))/(X1**2-X2**2)
      CEL=Y1-AEL*X1**2-BEL*X1
      IF (AEL.EQ.0..AND.BEL.EQ.0..AND.CEL.EQ.0.) CEL=1
      F1=(XNODE(1,IC)-LOC)/H1(I)
      F2=(XNODE(1,I)-LOC)/H1(I)
      F3=(XNODE(4,I)-LOC)/H1(I)
      F4=(XNODE(4,IC)-LOC)/H1(I)
      CNUM1=AEL/3*(F1**3-F4**3)+BEL/2*(F1**2-F4**2)+CEL*(F1-F4)
      CNUM2=AEL/3*(F2**3-F3**3)+BEL/2*(F2**2-F3**2)+CEL*(F2-F3)
      FF(I,K)=FF(I,K)+CNUM2/CNUM1*YFLUX(IC,K)
      ENDIF
      ENDDO
      ENDIF
      ENDDO
***                                                                 ***
***                                                                 ***
***********************************************************************
      DO K = 1,KBM1
#ifdef OMP
!$OMP PARALLEL DO PRIVATE(DUMT,DUMT1,DVMT,DVMT1,DTN,DTN1)
#endif  
      DO I = 1,N_CTRD
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
      
        
      FF(I,K) = (F(I,K)*(DTN+EL(I))*DJ(I)-DTI*FF(I,K))
     * /((DTN1+ELF(I))*DJ(I))      

      ENDIF
      ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif
      ENDDO
      DO II=1,N_CTRD_IVG
      I=NAG4VGC(II,1)
      IF (FSM(I).GT.0) THEN
      DO K=1,KBM1  
      FF(I,K)=0.
      DO J=1,9
          CNUM=NAG4VGC(II,J)
          FF(I,K)=FF(I,K)+DJ(CNUM)*D(CNUM)*DZ(K)*FSM(CNUM)*FF(CNUM,K)
      ENDDO
      FF(I,K)=FF(I,K)/DJ(I)/DZ(K)/D(I)
      ENDDO
      ENDIF
      ENDDO
      
      DO II=1,NNVGF4AG
      IC=NVGF4AG(II,1)
      IF(FSM(IC).GT.0.0) THEN
          DO J=1,9
              I=NVGF4AG(II,J+1)
              DO K=1,KBM1
                  FF(I,K)=FF(IC,K)
              ENDDO
          ENDDO
      ENDIF
      ENDDO
      
      RETURN
      END SUBROUTINE ADVT_HSIMT_TEST1
     
	END MODULE ADVT_TEST1
	  
	  