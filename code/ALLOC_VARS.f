#include "DEFS.h"
      
      SUBROUTINE ALLOC_VARS
      
      USE MOD_GLOBAL
      IMPLICIT NONE        
        
	
      !NEST
      ALLOCATE (NIM1(N_CTRD))          ;NIM1=0 
      ALLOCATE (NIM2(N_CTRD))          ;NIM2=0
      ALLOCATE (NIM3(N_CTRD))          ;NIM3=0
      ALLOCATE (NIM4(N_CTRD))          ;NIM4=0
      ALLOCATE (NIP1(N_CTRD))          ;NIP1=0 
      ALLOCATE (NIP2(N_CTRD))          ;NIP2=0
      ALLOCATE (NIP3(N_CTRD))          ;NIP3=0
      ALLOCATE (NIP4(N_CTRD))          ;NIP4=0
      ALLOCATE (NJM1(N_CTRD))          ;NJM1=0 
      ALLOCATE (NJM2(N_CTRD))          ;NJM2=0
      ALLOCATE (NJM3(N_CTRD))          ;NJM3=0
      ALLOCATE (NJM4(N_CTRD))          ;NJM4=0
      ALLOCATE (NJP1(N_CTRD))          ;NJP1=0 
      ALLOCATE (NJP2(N_CTRD))          ;NJP2=0
      ALLOCATE (NJP3(N_CTRD))          ;NJP3=0
      ALLOCATE (NJP4(N_CTRD))          ;NJP4=0
      ALLOCATE (NIM1JP1(N_CTRD))       ;NIM1JP1=0
      ALLOCATE (NIM1JM1(N_CTRD))       ;NIM1JM1=0
      ALLOCATE (NIP1JP1(N_CTRD))       ;NIP1JP1=0
      ALLOCATE (NIP1JM1(N_CTRD))       ;NIP1JM1=0
      ALLOCATE (NIM1JM2(N_CTRD))       ;NIM1JM2=0
      ALLOCATE (NIM2JM1(N_CTRD))       ;NIM2JM1=0
      
      ALLOCATE (SUMH(N_CTRDP1))       ;SUMH=0.0 !LXY
	ALLOCATE (TAUE0(N_CTRDP1))	   ;TAUE0=0.0
	ALLOCATE (TAUD0(N_CTRDP1))	   ;TAUD0=0.0
      ALLOCATE (H(N_CTRDP1))           ;H=0.0
      ALLOCATE (H1(N_CTRDP1))          ;H1=0.0
      ALLOCATE (H2(N_CTRDP1))          ;H2=0.0
      ALLOCATE (H3(N_CTRDP1))          ;H3=0.0
      ALLOCATE (DX1(N_CTRDP1))         ;DX1=0.
      ALLOCATE (DX2(N_CTRDP1))         ;DX2=0.
      ALLOCATE (DX3(N_CTRDP1))         ;DX3=0.
      ALLOCATE (DX4(N_CTRDP1))         ;DX4=0.
      ALLOCATE (DY1(N_CTRDP1))         ;DY1=0.
      ALLOCATE (DY2(N_CTRDP1))         ;DY2=0.
      ALLOCATE (DY3(N_CTRDP1))         ;DY3=0.
      ALLOCATE (DY4(N_CTRDP1))         ;DY4=0.
      ALLOCATE (D(N_CTRDP1))			;D=0.0
      ALLOCATE (ANG(N_CTRDP1))         ;ANG=0.0
      ALLOCATE (ARU(N_CTRDP1))        ;ARU=0.0
      ALLOCATE (ARV(N_CTRDP1))        ;ARV=0.0
      ALLOCATE (CBC_U(N_CTRDP1))      ;CBC_U=0.0
      ALLOCATE (CBC_V(N_CTRDP1))      ;CBC_V=0.0
      ALLOCATE (CBC_UV(N_CTRDP1))     ;CBC_UV=0.0
      ALLOCATE (CBC_COF(N_CTRDP1))    ;CBC_COF=0.0
      
      ALLOCATE (CBC_COF_A1(N_CTRDP1)) ;CBC_COF_A1=0.0
      ALLOCATE (CBC_COF_B1(N_CTRDP1)) ;CBC_COF_B1=0.0
      ALLOCATE (CBC_COF_A2(N_CTRDP1)) ;CBC_COF_A2=0.0
      ALLOCATE (CBC_COF_B2(N_CTRDP1)) ;CBC_COF_B2=0.0
      
      ALLOCATE (DUM(N_CTRDP1))         ;DUM=0.0
      ALLOCATE (DVM(N_CTRDP1))         ;DVM=0.0
      ALLOCATE (FSM(N_CTRDP1))         ;FSM=0.0
      ALLOCATE (FSM11(N_CTRDP1))       ;FSM11=0.0
      ALLOCATE (FSM12(N_CTRDP1))       ;FSM12=0.0
      ALLOCATE (FSMADD(N_CTRDP1))      ;FSMADD=0.0
      ALLOCATE (FSMTEMP(N_CTRDP1))     ;FSMTEMP = 1.
      ALLOCATE (COR(N_CTRDP1))         ;COR=0.0
      ALLOCATE (WTSURF(N_CTRD))        ;WTSURF=0.0
      ALLOCATE (WSSURF(N_CTRD))        ;WSSURF=0.0
      ALLOCATE (WUSURF(N_CTRDP1))      ;WUSURF=0.0
      ALLOCATE (WVSURF(N_CTRDP1))      ;WVSURF=0.0
      ALLOCATE (WUBOT(N_CTRDP1))       ;WUBOT=0.0
      ALLOCATE (WVBOT(N_CTRDP1))       ;WVBOT=0.0
      ALLOCATE (ELF(N_CTRDP1))         ;ELF=0.0
      ALLOCATE (EL(N_CTRDP1))          ;EL=0.0
      ALLOCATE (EL1(N_CTRDP1))         ;EL1=0.0
      ALLOCATE (RR(N_CTRDP1))          ;RR=0.0
      ALLOCATE (RRT(N_CTRDP1))         ;RRT=0.0
      ALLOCATE (AAUU(N_CTRDP1))        ;AAUU=0.0
      ALLOCATE (AAVV(N_CTRDP1))        ;AAVV=0.0
      ALLOCATE (AAZZ(N_CTRDP1))        ;AAZZ=0.0
      ALLOCATE (BBBB(N_CTRDP1))        ;BBBB=0.0
      ALLOCATE (HU(N_CTRDP1))             ;HU=0.0
      ALLOCATE (HV(N_CTRDP1))             ;HV=0.0
      ALLOCATE (DU(N_CTRDP1))             ;DU=0.0
      ALLOCATE (DV(N_CTRDP1))             ;DV=0.0
      ALLOCATE (DJ(N_CTRDP1))             ;DJ=0.0
      ALLOCATE (XXI(N_CTRDP1))            ;XXI=0.0
      ALLOCATE (YXI(N_CTRDP1))            ;YXI=0.0
      ALLOCATE (XETA(N_CTRDP1))           ;XETA=0.0
      ALLOCATE (YETA(N_CTRDP1))           ;YETA=0.0
      ALLOCATE (SWRAD(N_CTRDP1))          ;SWRAD=0.0
      ALLOCATE (HTEMP(N_CTRDP1))          ;HTEMP=0.0
      ALLOCATE (XR(N_CTRDP1))             ;XR=0.0
      ALLOCATE (YR(N_CTRDP1))             ;YR=0.0
      ALLOCATE (XU(N_CTRDP1))             ;XU=0.0
      ALLOCATE (YU(N_CTRDP1))             ;YU=0.0
      ALLOCATE (XV(N_CTRDP1))             ;XV=0.0
      ALLOCATE (YV(N_CTRDP1))             ;YV=0.0
      ALLOCATE (LON(N_CTRDP1))            ;LON=0.0
      ALLOCATE (LAT(N_CTRDP1))            ;LAT=0.0
      ALLOCATE (LONU(N_CTRDP1))           ;LONU=0.0
      ALLOCATE (LATU(N_CTRDP1))           ;LATU=0.0
      ALLOCATE (LONV(N_CTRDP1))           ;LONV=0.0
      ALLOCATE (LATV(N_CTRDP1))           ;LATV=0.0
      ALLOCATE (TPS(N_CTRDP1))            ;TPS=0.0
      ALLOCATE (DTXX(N_CTRDP1,2))         ;DTXX=0.0
      ALLOCATE (DTYY(N_CTRDP1,2))         ;DTYY=0.0
      ALLOCATE (A(N_CTRDP1,KB))           ;A=0.0
      ALLOCATE (C(N_CTRDP1,KB))           ;C=0.0
      ALLOCATE (VH(N_CTRDP1,KB))          ;VH=0.0
      ALLOCATE (VHP(N_CTRDP1,KB))         ;VHP=0.0
      ALLOCATE (PROD(N_CTRD,KB))          ;PROD=0.0
      ALLOCATE (DTEF(N_CTRDP1,KB))        ;DTEF=0.0
      ALLOCATE (T(N_CTRDP1,KB))           ;T=0.0
      ALLOCATE (S(N_CTRDP1,KB))           ;S=0.0
      ALLOCATE (KM(N_CTRDP1,KB))          ;KM=0.0
      ALLOCATE (KH(N_CTRDP1,KB))          ;KH=0.0
      ALLOCATE (KM_COF(N_CTRDP1))         ;KM_COF=0.0 !XK
      ALLOCATE (KH_COF(N_CTRDP1))         ;KH_COF=0.0 !XK
      ALLOCATE (KQ(N_CTRD,KB))            ;KQ=0.0
      ALLOCATE (L(N_CTRD_AG,KB))          ;L=0.0
      ALLOCATE (Q2(N_CTRDP1,KB))          ;Q2=0.0
      ALLOCATE (Q2L(N_CTRDP1,KB))         ;Q2L=0.0
      ALLOCATE (AAM(N_CTRDP1,KB))         ;AAM=0.0
      ALLOCATE (XMFLUX(N_CTRDP1,KB))      ;XMFLUX=0.0
      ALLOCATE (YMFLUX(N_CTRDP1,KB))      ;YMFLUX=0.0
      ALLOCATE (XMFLUX0(N_CTRDP1,KB))     ;XMFLUX0=0.0
      ALLOCATE (YMFLUX0(N_CTRDP1,KB))     ;YMFLUX0=0.0
      ALLOCATE (XFLUX(N_CTRDP1,KB))       ;XFLUX=0.0
      ALLOCATE (YFLUX(N_CTRDP1,KB))       ;YFLUX=0.0
      ALLOCATE (UF(N_CTRDP1,KB))          ;UF=0.0
      ALLOCATE (VF(N_CTRDP1,KB))          ;VF=0.0
      ALLOCATE (U(N_CTRDP1,KB))           ;U=0.0
      ALLOCATE (V(N_CTRDP1,KB))           ;V=0.0
      ALLOCATE (W(N_CTRDP1,KB))           ;W=0.0
      ALLOCATE (UR(N_CTRDP1,KB))          ;UR=0.0
      ALLOCATE (VR(N_CTRDP1,KB))          ;VR=0.0
      ALLOCATE (UC(N_CTRDP1,KB))          ;UC=0.0
      ALLOCATE (VC(N_CTRDP1,KB))          ;VC=0.0
      ALLOCATE (WR(N_CTRDP1,KB))          ;WR=0.0
      ALLOCATE (UNN(N_CTRDP1,KB))         ;UNN=0.0
      ALLOCATE (VNN(N_CTRDP1,KB))         ;VNN=0.0
      ALLOCATE (WNN(N_CTRDP1,KB))         ;WNN=0.0
      ALLOCATE (RHO(N_CTRDP1,KB))         ;RHO=0.0
      ALLOCATE (DRHOX(N_CTRD,KB))         ;DRHOX=0.0
      ALLOCATE (DRHOY(N_CTRD,KB))         ;DRHOY=0.0
      ALLOCATE (CURV4(N_CTRDP1,KB))       ;CURV4=0.0
      ALLOCATE (T0ZZ(N_CTRDP1,KSL))       ;T0ZZ=0.0
      ALLOCATE (S0ZZ(N_CTRDP1,KSL))       ;S0ZZ=0.0
      ALLOCATE (WINDU(N_CTRD))            ;WINDU=0.0
      ALLOCATE (WINDV(N_CTRD))            ;WINDV=0.0
      ALLOCATE (ATP(N_CTRD))              ;ATP=0.0
      ALLOCATE (RHM(N_CTRD))              ;RHM=0.0
      ALLOCATE (CLOUD(N_CTRD))            ;CLOUD=0.0
      ALLOCATE (APR(N_CTRD))              ;APR=0.0
      ALLOCATE (DATP(N_CTRD,2))           ;DATP=0.0
      ALLOCATE (DRHM(N_CTRD,2))           ;DRHM=0.0
      ALLOCATE (DCLOUD(N_CTRD,2))         ;DCLOUD=0.0
      ALLOCATE (DAPR(N_CTRD,2))           ;DAPR=0.0
      ALLOCATE (SED(N_CTRDP1,KB))         ;SED=0. 
      ALLOCATE (WFF(N_CTRDP1,KB))         ;WFF=0.
      ALLOCATE (DWS(N_CTRD,KB))           ;DWS=0.
      ALLOCATE (VNIAN(N_CTRD,KB))         ;VNIAN=0.
      ALLOCATE (USIN(N_CTRDP1))           ;USIN=0.
      ALLOCATE (TAU(N_CTRDP1))            ;TAU=0.
      ALLOCATE (TAU_WAVE(N_CTRDP1))       ;TAU_WAVE=0.
      ALLOCATE (TAU_TIDE(N_CTRDP1))       ;TAU_TIDE=0.
      ALLOCATE (TAU_U(N_CTRDP1))          ;TAU_U=0.
      ALLOCATE (TAU_V(N_CTRDP1))          ;TAU_V=0.
      ALLOCATE (TAUE(N_CTRDP1))           ;TAUE=0.
      ALLOCATE (TAUD(N_CTRDP1))           ;TAUD=0.
      ALLOCATE (QERO(N_CTRDP1))           ;QERO=0.
      ALLOCATE (QDEP(N_CTRDP1))           ;QDEP=0.
      ALLOCATE (QSI(N_CTRD_AG))           ;QSI=0.
      ALLOCATE (ZBED(N_CTRDP1))           ;ZBED=0.    ! DOUBLE PRECISION
      ALLOCATE (ZBEDD(N_CTRDP1))          ;ZBEDD=0. ! DOUBLE PRECISION  !!YR ADDED
      ALLOCATE (D50(N_CTRDP1))            ;D50=0.
      ALLOCATE (UBAR(N_CTRD,KB))          ;UBAR=0.
      ALLOCATE (VBAR(N_CTRD,KB))          ;VBAR=0.
      ALLOCATE (HSIG(N_CTRDP1))           ;HSIG=0.
      ALLOCATE (TSIG(N_CTRDP1))           ;TSIG=0.
      ALLOCATE (WAVEDIR(N_CTRDP1))        ;WAVEDIR=0.
      ALLOCATE (UBM(N_CTRDP1))            ;UBM=0.
      ALLOCATE (UBM_WB(N_CTRDP1))         ;UBM_WB=0.
      ALLOCATE (P(N_CTRDP1))              ;P=0.0
      ALLOCATE (AP(N_CTRDP1))             ;AP=0.0
      ALLOCATE (VARIAT1(N_CTRDP1))        ;VARIAT1=0.
      ALLOCATE (ACB_VARIAT(N_CTRDP1))     ;ACB_VARIAT=0.
      
	
      IF (JHM.NE.0) THEN
	  ALLOCATE (DTT22(N_CTRDP1,JHM))        ;DTT22=0.0
        ALLOCATE (DEI22(N_CTRDP1,JHM))        ;DEI22=0.0
	  ALLOCATE (ARCUX(N_CTRDP1,KB,JHM))     ;ARCUX=0.0
	  ALLOCATE (ARCVX(N_CTRDP1,KB,JHM))     ;ARCVX=0.0
	  ALLOCATE (ARCS(N_CTRDP1,KB,JHM))      ;ARCS=0.0
	  ALLOCATE (ARCET(N_CTRDP1,JHM))        ;ARCET=0.0
	  ALLOCATE (ARCUR(N_CTRDP1,KB,JHM))     ;ARCUR=0.0
	  ALLOCATE (ARCVR(N_CTRDP1,KB,JHM))     ;ARCVR=0.0
	  ALLOCATE (ARCSUX(N_CTRDP1,KB,JHM))    ;ARCSUX=0.0
	  ALLOCATE (ARCSVX(N_CTRDP1,KB,JHM))    ;ARCSVX=0.0
		
	  ALLOCATE (ARC_UA(N_CTRDP1,JHM))       ;ARC_UA=0.0
        ALLOCATE (ARC_CA(N_CTRDP1,JHM))       ;ARC_CA=0.0
        ALLOCATE (ARC_VA(N_CTRDP1,JHM))       ;ARC_VA=0.0
        ALLOCATE (ARC_UHH(N_CTRDP1,JHM))      ;ARC_UHH=0.0
        ALLOCATE (ARC_CHH(N_CTRDP1,JHM))      ;ARC_CHH=0.0
        ALLOCATE (ARC_UCHH(N_CTRDP1,JHM))     ;ARC_UCHH=0.0
        ALLOCATE (ARC_VHH(N_CTRDP1,JHM))      ;ARC_VHH=0.0
        ALLOCATE (ARC_VCHH(N_CTRDP1,JHM))     ;ARC_VCHH=0.0
        ALLOCATE (ARC_UCDH(N_CTRDP1,JHM))     ;ARC_UCDH=0.0
        ALLOCATE (ARC_VCDH(N_CTRDP1,JHM))     ;ARC_VCDH=0.0
        ALLOCATE (ARC_H(N_CTRDP1,JHM))        ;ARC_H=0.0
      
        ALLOCATE (FU1_WATER(N_CTRDP1,JHM))    ;FU1_WATER=0.0
        ALLOCATE (FU2_WATER(N_CTRDP1,JHM))    ;FU2_WATER=0.0
        ALLOCATE (FU1(N_CTRDP1,JHM))          ;FU1=0.0
        ALLOCATE (FU2(N_CTRDP1,JHM))          ;FU2=0.0
        ALLOCATE (FU3(N_CTRDP1,JHM))          ;FU3=0.0
        ALLOCATE (FU4(N_CTRDP1,JHM))          ;FU4=0.0
        ALLOCATE (FU5(N_CTRDP1,JHM))          ;FU5=0.0
        ALLOCATE (FV1_WATER(N_CTRDP1,JHM))    ;FV1_WATER=0.0
        ALLOCATE (FV2_WATER(N_CTRDP1,JHM))    ;FV2_WATER=0.0
        ALLOCATE (FV1(N_CTRDP1,JHM))          ;FV1=0.0
        ALLOCATE (FV2(N_CTRDP1,JHM))          ;FV2=0.0
        ALLOCATE (FV3(N_CTRDP1,JHM))          ;FV3=0.0
        ALLOCATE (FV4(N_CTRDP1,JHM))          ;FV4=0.0
        ALLOCATE (FV5(N_CTRDP1,JHM))          ;FV5=0.0
        ENDIF
 
        ALLOCATE (QS_RAD(N_CTRD_AG))       ;QS_RAD=0.0
        ALLOCATE (QB_RAD(N_CTRD_AG))       ;QB_RAD=0.0
        ALLOCATE (QE_RAD(N_CTRD_AG))       ;QE_RAD=0.0
        ALLOCATE (QH_RAD(N_CTRD_AG))       ;QH_RAD=0.0
        
      
 
      ALLOCATE (DWHT (N_CTRD,2))      ;DWHT =0.0
      ALLOCATE (DWPER(N_CTRD,2))      ;DWPER=0.0
      ALLOCATE (DWDIR(N_CTRD,2))      ;DWDIR=0.0
      ALLOCATE (DWLEN(N_CTRD,2))      ;DWLEN=0.0
      ALLOCATE (CONCENTRATION(N_CTRDP1,KB))  ;CONCENTRATION=0.0
     
      ALLOCATE (WAVEFFP(481))             ;WAVEFFP =0.0     

! CBR :      
!======================================================================
! INFO_EXCH 

#ifdef INT_UNI

! HSIMT PARABOLIC INTERPOLATION
#elif defined INT_HPI
      ALLOCATE (MASK4VAR(8,N_CTRD_EVG))  ;MASK4VAR=0.0
! HSIMT AEI (HSIMT ADVECTION-EQUIVALENT INTERPOLATION)
#elif defined INT_HAEI
      ALLOCATE (MASK4VAR(10,N_CTRD_EVG))  ;MASK4VAR=0.0
#endif

      
! END: INFO_EXCH
!======================================================================
! PROFQ
      ALLOCATE ( GH(N_CTRD,KB),SM(N_CTRD,KB),SH(N_CTRD,KB) )
      ALLOCATE ( KN(N_CTRD,KB),BOYGR(N_CTRD,KB) )
      GH=0.0
      SM=0.0
      SH=0.0
      KN=0.0
      BOYGR=0.0

! END: PROFQ
!======================================================================
! ADVT_HSIMT
      ALLOCATE ( KAX(N_CTRDP1,KB),KAY(N_CTRDP1,KB),KAZ(N_CTRDP1,KB) )
      KAX = 1.
      KAY = 1.
      KAZ = 1.
! END: ADVT_HSIMT
!======================================================================
! PROFT
      ALLOCATE ( RAD(N_CTRD_AG,KB) )
      RAD=0.0
! END: PROFT
!======================================================================
      
      
	RETURN
	END SUBROUTINE ALLOC_VARS