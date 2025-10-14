#include "DEFS.h"
	SUBROUTINE PRE_IE
	
!======================================================================
! This program is used to prepare for the info exchange process
!
! For update process, available choices include:
! UD_DR (directly-replacing), UD_IDW (inverse distance weighting interpolation), 
! UD_AVE (area-averaging), UD_SF (9-point Shapiro filtering),
! UD_FWO (full-weighting operator)
!
! For interpolation precess, available choices include:
! INT_UNI (0th order uniform), INT_IDW (inverse distance weighting), 
! INT_LINEAR (inverse bilinear interpolation), INT_QDT (quadratic interpolation),
! INT_HPI (HSIMT parabolic interpolation), INT_UAEI (upwind AEI), 
! INT_HAEI (HSIMT AEI)
! If choosing INT_QDT or INT_UAEI,it is recommended to set the refining multiple to odd numbers
!
! by Sylvers Ding, Dec. 30, 2019	
!======================================================================
	
	USE MOD_GLOBAL
	
	IMPLICIT NONE
	 
	INTEGER I,J,K,II
	INTEGER IND,NREF,N_AG
	INTEGER ITEMP(1)
      INTEGER FindCtrd
	REAL, PARAMETER :: PI = 3.1415926/180.
      REAL, PARAMETER :: EPSILON = 1E-10
	REAL DX, DY, DIS, DIS1, DIS2
	REAL XB(4),YB(4)
	REAL DISC(N_CTRD_AG), DISU(N_CTRD_AG), DISV(N_CTRD_AG)
	
	! SF !SMOOTHING ELEMENT shapiro:0.5 full-weighting:2./3.
#if defined UD_SF
	REAL, PARAMETER :: SE = 1./2. 
#elif defined UD_FWO
      REAL, PARAMETER :: SE = 2./3. 
#endif
	
#if defined INT_LINEAR
	REAL AT1,AT2,BT1,BT2,CT1,CT2
	REAL EI,EJ,FI,FJ,GI,GJ,HI,HJ
	REAL K0,K1,K2
	REAL UL,VL
#endif

c#if defined INT_QDT || defined INT_UAEI || \
c    defined INT_HPI || defined INT_HAEI
#if defined INT_QDT || defined INT_UAEI || defined INT_HPI || defined INT_HAEI
	INTEGER RFN
	REAL A1, A2, B1, B2, C1, C2
	REAL AT, BT, CT
	REAL XSEC, YSEC
	REAL EPS, ALPHA
#endif
	
#ifdef INT_HPI
      ALLOCATE (XIC(N_CTRD_EVG))     ;XIC = -9999.
      ALLOCATE (XIU(N_CTRD_EVG))     ;XIU = -9999.
      ALLOCATE (XIV(N_CTRD_EVG))     ;XIV = -9999.
#endif

#ifdef INT_HAEI
      ALLOCATE (EPSC(N_CTRD_EVG))   ;EPSC = -9999.
      ALLOCATE (EPSU(N_CTRD_EVG))   ;EPSU = -9999.
      ALLOCATE (EPSV(N_CTRD_EVG))   ;EPSV = -9999.
#endif

	PRINT*, 'PREPARE FOR THE INFO EXCHANGE...'
	
	VTP = 0
	
	ALLOCATE (NAGC_EVG(NMAX,N_CTRD_EVG))	  ;NAGC_EVG = 1 !DONT CHANGE
	ALLOCATE (NAGU_EVG(NMAX,N_CTRD_EVG))	  ;NAGU_EVG = 1
	ALLOCATE (NAGV_EVG(NMAX,N_CTRD_EVG))	  ;NAGV_EVG = 1
	ALLOCATE (WTC_EVG(NMAX,N_CTRD_EVG))     ;WTC_EVG = 0.
	ALLOCATE (WTU_EVG(NMAX,N_CTRD_EVG))     ;WTU_EVG = 0.
	ALLOCATE (WTV_EVG(NMAX,N_CTRD_EVG))     ;WTV_EVG = 0.
	ALLOCATE (NAGC_IVG(NMAX,N_CTRD_IVG))	  ;NAGC_IVG = 1
	ALLOCATE (NAGU_IVG(NMAX,N_CTRD_IVG))	  ;NAGU_IVG = 1
	ALLOCATE (NAGV_IVG(NMAX,N_CTRD_IVG))	  ;NAGV_IVG = 1
	ALLOCATE (WTC_IVG(NMAX,N_CTRD_IVG))     ;WTC_IVG = 0.
	ALLOCATE (WTU_IVG(NMAX,N_CTRD_IVG))     ;WTU_IVG = 0.
	ALLOCATE (WTV_IVG(NMAX,N_CTRD_IVG))     ;WTV_IVG = 0.

c#if defined UD_IDW || defined INT_IDW  || \
c    defined INT_QDT  || defined INT_LINEAR || \
c    defined INT_UAEI
#if defined UD_IDW || defined INT_IDW || defined INT_QDT  || defined INT_LINEAR || defined INT_UAEI
	  ALLOCATE (ELM(6)) ;ELM = 0.
	  NPOW = 1
	  SEARCHRADIUS = MAX(MAXVAL(H1),MAXVAL(H2))*1.5
#endif

!======================================================================
!                         INTERPOLATION
!======================================================================
	
!**********************************************************************
! 0TH ORDER UNIFORM
#if defined INT_UNI
	  DO I = 1,N_CTRD_EVG		
		DO J = 1,N_CTRD_AG
		  DO K = 1,4
			XB(K) = XNODE(K,J)
			YB(K) = YNODE(K,J)
		  ENDDO
		  
		  IF (COORD.EQ.'XY') THEN
			CALL INSIDE(XR(N_CTRD_AG+I),YR(N_CTRD_AG+I),XB,YB,NB,IND)
		  ELSE
			CALL INSIDE(LON(N_CTRD_AG+I),LAT(N_CTRD_AG+I),XB,YB,NB,IND)
		  ENDIF
		  
		  IF (IND.EQ.1) THEN
			NAGC_EVG(1,I) = J
			EXIT
		  ELSEIF (J.EQ.N_CTRD_AG) THEN
			PRINT*, 'CANNOT FIND MATCH GRID FOR EVG',N_CTRD_AG+I
			PAUSE
			STOP
		  ENDIF
		ENDDO
	  ENDDO
	  NAGU_EVG = NAGC_EVG
	  NAGV_EVG = NAGC_EVG
***********************************************************************
!INVERSE DISTANCE WEIGHTING INTERPOLATION
#elif defined INT_IDW
      IF (COORD .EQ. 'XY') THEN
      
      DO I = 1,N_CTRD_EVG
        CALL WEIGHT_IDW(XR(N_CTRD_AG+I),YR(N_CTRD_AG+I),XR,YR
     * ,WTC_EVG(:,I),NAGC_EVG(:,I))
        CALL WEIGHT_IDW(XU(N_CTRD_AG+I),YU(N_CTRD_AG+I),XU,YU
     * ,WTU_EVG(:,I),NAGU_EVG(:,I))
        CALL WEIGHT_IDW(XV(N_CTRD_AG+I),YV(N_CTRD_AG+I),XV,YV
     * ,WTV_EVG(:,I),NAGV_EVG(:,I))
      ENDDO
      
	ELSE
        
      DO I = 1,N_CTRD_EVG
        CALL WEIGHT_IDW(LON(N_CTRD_AG+I),LAT(N_CTRD_AG+I),LON,LAT
     * ,WTC_EVG(:,I),NAGC_EVG(:,I))
        CALL WEIGHT_IDW(LONU(N_CTRD_AG+I),LATU(N_CTRD_AG+I),LONU,LATU
     * ,WTU_EVG(:,I),NAGU_EVG(:,I))
        CALL WEIGHT_IDW(LONV(N_CTRD_AG+I),LATV(N_CTRD_AG+I),LONV,LATV
     * ,WTV_EVG(:,I),NAGV_EVG(:,I))
      ENDDO
        
	ENDIF
!**********************************************************************
! INVERSE BILINEAR INTERPOLATION
#elif defined INT_LINEAR
      DO I = 1,N_CTRD_EVG
	  !FIND IN WHICH AG IS THIS EVG
	  DO J = 1,N_CTRD_AG
		DO K = 1,4
		  XB(K) = XNODE(K,J)
		  YB(K) = YNODE(K,J)
		ENDDO
		  
		IF (COORD.EQ.'XY') THEN
		  CALL INSIDE(XR(N_CTRD_AG+I),YR(N_CTRD_AG+I),XB,YB,NB,IND)
		ELSE
		  CALL INSIDE(LON(N_CTRD_AG+I),LAT(N_CTRD_AG+I),XB,YB,NB,IND)
		ENDIF
		  
		IF (IND.EQ.1) THEN
		  NAGC_EVG(1,I) = J
		  NAGU_EVG(1,I) = J
		  NAGV_EVG(1,I) = J
		  EXIT
		ELSEIF (J.EQ.N_CTRD_AG) THEN
		  PRINT*, 'CANNOT FIND MATCH GRID FOR EVG',N_CTRD_AG+I
		  PAUSE
		  STOP
		ENDIF
	  ENDDO
	  
C----------------------------------------------------------------------
C NAGC & WTC
	  !DETERMINE AG NUMBERS FOR INTERPOLATION & WEIGHTS
	  IF ((NIP2(NAGC_EVG(1,I)).GT.N_CTRD_AG) .OR. 
     * (NIM2(NAGC_EVG(1,I)).GT.N_CTRD_AG)) THEN !WEST OR EAST
		
		IF (NIM1(NAGC_EVG(1,I)).EQ.N_CTRDP1) THEN
		  NAGC_EVG(:,I) = 1
		  NAGU_EVG(:,I) = 1
		  NAGV_EVG(:,I) = 1
		  CYCLE
		ENDIF
          
		IF (COORD.EQ.'XY') THEN
		  AT1 = YR(NIP1(NAGC_EVG(1,I)))-YR(NAGC_EVG(1,I))
		  BT1 = XR(NAGC_EVG(1,I))-XR(NIP1(NAGC_EVG(1,I)))
		  CT1 = 0.
		  DIS1 = (YR(N_CTRD_AG+I)-YR(NAGC_EVG(1,I)))
     * -(-AT1/BT1*(XR(N_CTRD_AG+I)-XR(NAGC_EVG(1,I)))-CT1/BT1)
		  
		  AT2 = YR(NJP1(NAGC_EVG(1,I)))-YR(NAGC_EVG(1,I))
		  BT2 = XR(NAGC_EVG(1,I))-XR(NJP1(NAGC_EVG(1,I)))
		  CT2 = 0.
		  DIS2 = (XR(N_CTRD_AG+I)-XR(NAGC_EVG(1,I)))
     * -(-BT2/AT2*(YR(N_CTRD_AG+I)-YR(NAGC_EVG(1,I)))-CT2/AT2)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT1 = LAT(NIP1(NAGC_EVG(1,I)))-LAT(NAGC_EVG(1,I))
		  BT1 = LON(NAGC_EVG(1,I))-LON(NIP1(NAGC_EVG(1,I)))
		  CT1 = 0.
		  DIS1 = (LAT(N_CTRD_AG+I)-LAT(NAGC_EVG(1,I)))
     * -(-AT1/BT1*(LON(N_CTRD_AG+I)-LON(NAGC_EVG(1,I)))-CT1/BT1)
		  
		  AT2 = LAT(NJP1(NAGC_EVG(1,I)))-LAT(NAGC_EVG(1,I))
		  BT2 = LON(NAGC_EVG(1,I))-LON(NJP1(NAGC_EVG(1,I)))
		  CT2 = 0.
		  DIS2 = (LON(N_CTRD_AG+I)-LON(NAGC_EVG(1,I)))
     * -(-BT2/AT2*(LAT(N_CTRD_AG+I)-LAT(NAGC_EVG(1,I)))-CT2/AT2)
		ENDIF
		
		IF (DIS1.GE.0 .AND. DIS2.GE.0) THEN
		  NAGC_EVG(2,I) = NJP1(NAGC_EVG(1,I))
		  NAGC_EVG(3,I) = NIP1(NAGC_EVG(1,I))
		  NAGC_EVG(4,I) = NJP1(NIP1(NAGC_EVG(1,I)))
		ELSEIF (DIS1.GE.0 .AND. DIS2.LT.0) THEN
		  NAGC_EVG(2,I) = NJP1(NAGC_EVG(1,I))
		  NAGC_EVG(3,I) = NIM1(NAGC_EVG(1,I))
		  NAGC_EVG(4,I) = NJP1(NIM1(NAGC_EVG(1,I)))
		ELSEIF (DIS1.LT.0 .AND. DIS2.GE.0) THEN
		  NAGC_EVG(2,I) = NJM1(NAGC_EVG(1,I))
		  NAGC_EVG(3,I) = NIP1(NAGC_EVG(1,I))
		  NAGC_EVG(4,I) = NJM1(NIP1(NAGC_EVG(1,I)))
		ELSEIF (DIS1.LT.0 .AND. DIS2.LT.0) THEN
		  NAGC_EVG(2,I) = NJM1(NAGC_EVG(1,I))
		  NAGC_EVG(3,I) = NIM1(NAGC_EVG(1,I))
		  NAGC_EVG(4,I) = NJM1(NIM1(NAGC_EVG(1,I)))
		ENDIF
		
        ELSE !NORTH, SOUTH OR OTHERS
     
          IF (COORD.EQ.'XY') THEN
		  AT1 = YR(NIP1(NAGC_EVG(1,I)))-YR(NAGC_EVG(1,I))
		  BT1 = XR(NAGC_EVG(1,I))-XR(NIP1(NAGC_EVG(1,I)))
		  CT1 = 0.
		  DIS1 = (YR(N_CTRD_AG+I)-YR(NAGC_EVG(1,I)))
     * -(-AT1/BT1*(XR(N_CTRD_AG+I)-XR(NAGC_EVG(1,I)))-CT1/BT1)
		  
		  AT2 = YR(NJP1(NAGC_EVG(1,I)))-YR(NAGC_EVG(1,I))
		  BT2 = XR(NAGC_EVG(1,I))-XR(NJP1(NAGC_EVG(1,I)))
		  CT2 = 0.
		  DIS2 = (XR(N_CTRD_AG+I)-XR(NAGC_EVG(1,I)))
     * -(-BT2/AT2*(YR(N_CTRD_AG+I)-YR(NAGC_EVG(1,I)))-CT2/AT2)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT1 = LAT(NIP1(NAGC_EVG(1,I)))-LAT(NAGC_EVG(1,I))
		  BT1 = LON(NAGC_EVG(1,I))-LON(NIP1(NAGC_EVG(1,I)))
		  CT1 = 0.
		  DIS1 = (LAT(N_CTRD_AG+I)-LAT(NAGC_EVG(1,I)))
     * -(-AT1/BT1*(LON(N_CTRD_AG+I)-LON(NAGC_EVG(1,I)))-CT1/BT1)
		  
		  AT2 = LAT(NJP1(NAGC_EVG(1,I)))-LAT(NAGC_EVG(1,I))
		  BT2 = LON(NAGC_EVG(1,I))-LON(NJP1(NAGC_EVG(1,I)))
		  CT2 = 0.
		  DIS2 = (LON(N_CTRD_AG+I)-LON(NAGC_EVG(1,I)))
     * -(-BT2/AT2*(LAT(N_CTRD_AG+I)-LAT(NAGC_EVG(1,I)))-CT2/AT2)
		ENDIF
		
		IF (DIS1.GE.0 .AND. DIS2.GE.0) THEN
		  NAGC_EVG(2,I) = NIP1(NAGC_EVG(1,I))
		  NAGC_EVG(3,I) = NJP1(NAGC_EVG(1,I))
		  NAGC_EVG(4,I) = NIP1(NJP1(NAGC_EVG(1,I)))
		ELSEIF (DIS1.GE.0 .AND. DIS2.LT.0) THEN
		  NAGC_EVG(2,I) = NIM1(NAGC_EVG(1,I))
		  NAGC_EVG(3,I) = NJP1(NAGC_EVG(1,I))
		  NAGC_EVG(4,I) = NIM1(NJP1(NAGC_EVG(1,I)))
		ELSEIF (DIS1.LT.0 .AND. DIS2.GE.0) THEN
		  NAGC_EVG(2,I) = NIP1(NAGC_EVG(1,I))
		  NAGC_EVG(3,I) = NJM1(NAGC_EVG(1,I))
		  NAGC_EVG(4,I) = NIP1(NJM1(NAGC_EVG(1,I)))
		ELSEIF (DIS1.LT.0 .AND. DIS2.LT.0) THEN
		  NAGC_EVG(2,I) = NIM1(NAGC_EVG(1,I))
		  NAGC_EVG(3,I) = NJM1(NAGC_EVG(1,I))
		  NAGC_EVG(4,I) = NIM1(NJM1(NAGC_EVG(1,I)))
		ENDIF
	  
	  ENDIF
	  
	    DO J = 2,4
		  IF (NAGC_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XR(NAGC_EVG(J,I)),YR(NAGC_EVG(J,I))
     * ,XR(1:N_CTRD_AG),YR(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LON(NAGC_EVG(J,I)),LAT(NAGC_EVG(J,I))
     * ,LON(1:N_CTRD_AG),LAT(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGC_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		IF (COORD.EQ.'XY') THEN
		  EI = XR(NAGC_EVG(2,I))-XR(NAGC_EVG(1,I))
		  EJ = YR(NAGC_EVG(2,I))-YR(NAGC_EVG(1,I))
		  FI = XR(NAGC_EVG(3,I))-XR(NAGC_EVG(1,I))
		  FJ = YR(NAGC_EVG(3,I))-YR(NAGC_EVG(1,I))
		  GI = XR(NAGC_EVG(1,I))-XR(NAGC_EVG(2,I))
     * +XR(NAGC_EVG(4,I))-XR(NAGC_EVG(3,I))
		  GJ = YR(NAGC_EVG(1,I))-YR(NAGC_EVG(2,I))
     * +YR(NAGC_EVG(4,I))-YR(NAGC_EVG(3,I))
		  HI = XR(N_CTRD_AG+I)-XR(NAGC_EVG(1,I))
		  HJ = YR(N_CTRD_AG+I)-YR(NAGC_EVG(1,I))
		  
		ELSE 
		  EI = LON(NAGC_EVG(2,I))-LON(NAGC_EVG(1,I))
		  EJ = LAT(NAGC_EVG(2,I))-LAT(NAGC_EVG(1,I))
		  FI = LON(NAGC_EVG(3,I))-LON(NAGC_EVG(1,I))
		  FJ = LAT(NAGC_EVG(3,I))-LAT(NAGC_EVG(1,I))
		  GI = LON(NAGC_EVG(1,I))-LON(NAGC_EVG(2,I))
     * +LON(NAGC_EVG(4,I))-LON(NAGC_EVG(3,I))
		  GJ = LAT(NAGC_EVG(1,I))-LAT(NAGC_EVG(2,I))
     * +LAT(NAGC_EVG(4,I))-LAT(NAGC_EVG(3,I))
		  HI = LON(N_CTRD_AG+I)-LON(NAGC_EVG(1,I))
		  HJ = LAT(N_CTRD_AG+I)-LAT(NAGC_EVG(1,I))
		  
		ENDIF
		
		K2 = GI*FJ-GJ*FI
		K1 = EI*FJ-EJ*FI+HI*GJ-HJ*GI
		K0 = HI*EJ-HJ*EI
        
        IF (K2.EQ.0.) THEN
          VL = -K0/K1
          IF (VL.GE.-0.1 .AND. VL.LE.1.1) THEN
            CONTINUE
          ELSE
				PRINT*, 'WRONG VL = -K0/K1 IN WTC!'
				PRINT*, VL
				PRINT*, N_CTRD_AG+I
				PAUSE
				STOP
          ENDIF
        ELSE
          VL = (-K1+SQRT(K1**2-4*K0*K2))/(2*K2)
          IF (VL.GE.-0.1 .AND. VL.LE.1.1) THEN
				CONTINUE
			ELSE
				VL = (-K1-SQRT(K1**2-4*K0*K2))/(2*K2)
            IF (VL.GE.-0.1 .AND. VL.LE.1.1) THEN
              CONTINUE
            ELSE
					PRINT*, 'WRONG VL IN WTC!'
					PRINT*, VL
					PRINT*, N_CTRD_AG+I
					PAUSE
					STOP
            ENDIF
			ENDIF
        ENDIF
     
		UL = (HI-FI*VL)/(EI+GI*VL)
        IF (UL.GE.-0.1 .AND. UL.LE.1.1) THEN
          CONTINUE
        ELSE
          PRINT*,'WRONG UL IN WTC!' 
          PRINT*, UL,VL
          PRINT*, HI,FI,EI,GI
			PRINT*, N_CTRD_AG+I
			PAUSE
			STOP
        ENDIF
		
        WTC_EVG(1,I) = 1-UL-VL+UL*VL
		WTC_EVG(2,I) = UL*(1-VL)
		WTC_EVG(3,I) = VL*(1-UL)
		WTC_EVG(4,I) = UL*VL
		
		
C----------------------------------------------------------------------
C NAGU & WTU
	  !DETERMINE AG NUMBERS FOR INTERPOLATION & WEIGHTS
	  IF ((NIP2(NAGU_EVG(1,I)).GT.N_CTRD_AG) .OR. 
     * (NIM2(NAGU_EVG(1,I)).GT.N_CTRD_AG)) THEN !WEST OR EAST
          
		IF (COORD.EQ.'XY') THEN
		  AT1 = YU(NIP1(NAGU_EVG(1,I)))-YU(NAGU_EVG(1,I))
		  BT1 = XU(NAGU_EVG(1,I))-XU(NIP1(NAGU_EVG(1,I)))
		  CT1 = 0.
		  DIS1 = (YU(N_CTRD_AG+I)-YU(NAGU_EVG(1,I)))
     * -(-AT1/BT1*(XU(N_CTRD_AG+I)-XU(NAGU_EVG(1,I)))-CT1/BT1)
		  
		  AT2 = YU(NJP1(NAGU_EVG(1,I)))-YU(NAGU_EVG(1,I))
		  BT2 = XU(NAGU_EVG(1,I))-XU(NJP1(NAGU_EVG(1,I)))
		  CT2 = 0.
		  DIS2 = (XU(N_CTRD_AG+I)-XU(NAGU_EVG(1,I)))
     * -(-BT2/AT2*(YU(N_CTRD_AG+I)-YU(NAGU_EVG(1,I)))-CT2/AT2)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT1 = LATU(NIP1(NAGU_EVG(1,I)))-LATU(NAGU_EVG(1,I))
		  BT1 = LONU(NAGU_EVG(1,I))-LONU(NIP1(NAGU_EVG(1,I)))
		  CT1 = 0.
		  DIS1 = (LATU(N_CTRD_AG+I)-LATU(NAGU_EVG(1,I)))
     * -(-AT1/BT1*(LONU(N_CTRD_AG+I)-LONU(NAGU_EVG(1,I)))-CT1/BT1)
		  
		  AT2 = LATU(NJP1(NAGU_EVG(1,I)))-LATU(NAGU_EVG(1,I))
		  BT2 = LONU(NAGU_EVG(1,I))-LONU(NJP1(NAGU_EVG(1,I)))
		  CT2 = 0.
		  DIS2 = (LONU(N_CTRD_AG+I)-LONU(NAGU_EVG(1,I)))
     * -(-BT2/AT2*(LATU(N_CTRD_AG+I)-LATU(NAGU_EVG(1,I)))-CT2/AT2)
		ENDIF
		
		IF (DIS1.GE.0 .AND. DIS2.GE.0) THEN
		  NAGU_EVG(2,I) = NJP1(NAGU_EVG(1,I))
		  NAGU_EVG(3,I) = NIP1(NAGU_EVG(1,I))
		  NAGU_EVG(4,I) = NJP1(NIP1(NAGU_EVG(1,I)))
		ELSEIF (DIS1.GE.0 .AND. DIS2.LT.0) THEN
		  NAGU_EVG(2,I) = NJP1(NAGU_EVG(1,I))
		  NAGU_EVG(3,I) = NIM1(NAGU_EVG(1,I))
		  NAGU_EVG(4,I) = NJP1(NIM1(NAGU_EVG(1,I)))
		ELSEIF (DIS1.LT.0 .AND. DIS2.GE.0) THEN
		  NAGU_EVG(2,I) = NJM1(NAGU_EVG(1,I))
		  NAGU_EVG(3,I) = NIP1(NAGU_EVG(1,I))
		  NAGU_EVG(4,I) = NJM1(NIP1(NAGU_EVG(1,I)))
		ELSEIF (DIS1.LT.0 .AND. DIS2.LT.0) THEN
		  NAGU_EVG(2,I) = NJM1(NAGU_EVG(1,I))
		  NAGU_EVG(3,I) = NIM1(NAGU_EVG(1,I))
		  NAGU_EVG(4,I) = NJM1(NIM1(NAGU_EVG(1,I)))
		ENDIF
		
        ELSE !NORTH, SOUTH OR OTHERS
     
        IF (COORD.EQ.'XY') THEN
		  AT1 = YU(NIP1(NAGU_EVG(1,I)))-YU(NAGU_EVG(1,I))
		  BT1 = XU(NAGU_EVG(1,I))-XU(NIP1(NAGU_EVG(1,I)))
		  CT1 = 0.
		  DIS1 = (YU(N_CTRD_AG+I)-YU(NAGU_EVG(1,I)))
     * -(-AT1/BT1*(XU(N_CTRD_AG+I)-XU(NAGU_EVG(1,I)))-CT1/BT1)
		  
		  AT2 = YU(NJP1(NAGU_EVG(1,I)))-YU(NAGU_EVG(1,I))
		  BT2 = XU(NAGU_EVG(1,I))-XU(NJP1(NAGU_EVG(1,I)))
		  CT2 = 0.
		  DIS2 = (XU(N_CTRD_AG+I)-XU(NAGU_EVG(1,I)))
     * -(-BT2/AT2*(YU(N_CTRD_AG+I)-YU(NAGU_EVG(1,I)))-CT2/AT2)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT1 = LATU(NIP1(NAGU_EVG(1,I)))-LATU(NAGU_EVG(1,I))
		  BT1 = LONU(NAGU_EVG(1,I))-LONU(NIP1(NAGU_EVG(1,I)))
		  CT1 = 0.
		  DIS1 = (LATU(N_CTRD_AG+I)-LATU(NAGU_EVG(1,I)))
     * -(-AT1/BT1*(LONU(N_CTRD_AG+I)-LONU(NAGU_EVG(1,I)))-CT1/BT1)
		  
		  AT2 = LATU(NJP1(NAGU_EVG(1,I)))-LATU(NAGU_EVG(1,I))
		  BT2 = LONU(NAGU_EVG(1,I))-LONU(NJP1(NAGU_EVG(1,I)))
		  CT2 = 0.
		  DIS2 = (LONU(N_CTRD_AG+I)-LONU(NAGU_EVG(1,I)))
     * -(-BT2/AT2*(LATU(N_CTRD_AG+I)-LATU(NAGU_EVG(1,I)))-CT2/AT2)
		ENDIF
		
		IF (DIS1.GE.0 .AND. DIS2.GE.0) THEN
		  NAGU_EVG(2,I) = NIP1(NAGU_EVG(1,I))
		  NAGU_EVG(3,I) = NJP1(NAGU_EVG(1,I))
		  NAGU_EVG(4,I) = NIP1(NJP1(NAGU_EVG(1,I)))
		ELSEIF (DIS1.GE.0 .AND. DIS2.LT.0) THEN
		  NAGU_EVG(2,I) = NIM1(NAGU_EVG(1,I))
		  NAGU_EVG(3,I) = NJP1(NAGU_EVG(1,I))
		  NAGU_EVG(4,I) = NIM1(NJP1(NAGU_EVG(1,I)))
		ELSEIF (DIS1.LT.0 .AND. DIS2.GE.0) THEN
		  NAGU_EVG(2,I) = NIP1(NAGU_EVG(1,I))
		  NAGU_EVG(3,I) = NJM1(NAGU_EVG(1,I))
		  NAGU_EVG(4,I) = NIP1(NJM1(NAGU_EVG(1,I)))
		ELSEIF (DIS1.LT.0 .AND. DIS2.LT.0) THEN
		  NAGU_EVG(2,I) = NIM1(NAGU_EVG(1,I))
		  NAGU_EVG(3,I) = NJM1(NAGU_EVG(1,I))
		  NAGU_EVG(4,I) = NIM1(NJM1(NAGU_EVG(1,I)))
		ENDIF
	  
	  ENDIF
	  
	    DO J = 2,4
		  IF (NAGU_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XU(NAGU_EVG(J,I)),YU(NAGU_EVG(J,I))
     * ,XU(1:N_CTRD_AG),YU(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LONU(NAGU_EVG(J,I)),LATU(NAGU_EVG(J,I))
     * ,LONU(1:N_CTRD_AG),LATU(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGU_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		IF (COORD.EQ.'XY') THEN
		  EI = XU(NAGU_EVG(2,I))-XU(NAGU_EVG(1,I))
		  EJ = YU(NAGU_EVG(2,I))-YU(NAGU_EVG(1,I))
		  FI = XU(NAGU_EVG(3,I))-XU(NAGU_EVG(1,I))
		  FJ = YU(NAGU_EVG(3,I))-YU(NAGU_EVG(1,I))
		  GI = XU(NAGU_EVG(1,I))-XU(NAGU_EVG(2,I))
     * +XU(NAGU_EVG(4,I))-XU(NAGU_EVG(3,I))
		  GJ = YU(NAGU_EVG(1,I))-YU(NAGU_EVG(2,I))
     * +YU(NAGU_EVG(4,I))-YU(NAGU_EVG(3,I))
		  HI = XU(N_CTRD_AG+I)-XU(NAGU_EVG(1,I))
		  HJ = YU(N_CTRD_AG+I)-YU(NAGU_EVG(1,I))
		  
		ELSE 
		  EI = LONU(NAGU_EVG(2,I))-LONU(NAGU_EVG(1,I))
		  EJ = LATU(NAGU_EVG(2,I))-LATU(NAGU_EVG(1,I))
		  FI = LONU(NAGU_EVG(3,I))-LONU(NAGU_EVG(1,I))
		  FJ = LATU(NAGU_EVG(3,I))-LATU(NAGU_EVG(1,I))
		  GI = LONU(NAGU_EVG(1,I))-LONU(NAGU_EVG(2,I))
     * +LONU(NAGU_EVG(4,I))-LONU(NAGU_EVG(3,I))
		  GJ = LATU(NAGU_EVG(1,I))-LATU(NAGU_EVG(2,I))
     * +LATU(NAGU_EVG(4,I))-LATU(NAGU_EVG(3,I))
		  HI = LONU(N_CTRD_AG+I)-LONU(NAGU_EVG(1,I))
		  HJ = LATU(N_CTRD_AG+I)-LATU(NAGU_EVG(1,I))
		  
		ENDIF
		
		K2 = GI*FJ-GJ*FI
		K1 = EI*FJ-EJ*FI+HI*GJ-HJ*GI
		K0 = HI*EJ-HJ*EI
		
		IF (K2.EQ.0.) THEN
          VL = -K0/K1
          IF (VL.GE.-0.1 .AND. VL.LE.1.1) THEN
            CONTINUE
          ELSE
				PRINT*, 'WRONG VL = -K0/K1 IN WTU!'
				PRINT*, VL
				PRINT*, N_CTRD_AG+I
				PAUSE
				STOP
          ENDIF
        ELSE
          VL = (-K1+SQRT(K1**2-4*K0*K2))/(2*K2)
          IF (VL.GE.-0.1 .AND. VL.LE.1.1) THEN
				CONTINUE
			ELSE
				VL = (-K1-SQRT(K1**2-4*K0*K2))/(2*K2)
            IF (VL.GE.-0.1 .AND. VL.LE.1.1) THEN
              CONTINUE
            ELSE
					PRINT*, 'WRONG VL IN WTU!'
					PRINT*, VL
					PRINT*, N_CTRD_AG+I
					PAUSE
					STOP
            ENDIF
			ENDIF
        ENDIF
     
		UL = (HI-FI*VL)/(EI+GI*VL)
        IF (UL.GE.-0.1 .AND. UL.LE.1.1) THEN
          CONTINUE
        ELSE
          PRINT*,'WRONG UL IN WTU!' 
          PRINT*, UL,VL
          PRINT*, HI,FI,EI,GI
			PRINT*, N_CTRD_AG+I
			PAUSE
			STOP
        ENDIF
		
        WTU_EVG(1,I) = 1-UL-VL+UL*VL
		WTU_EVG(2,I) = UL*(1-VL)
		WTU_EVG(3,I) = VL*(1-UL)
		WTU_EVG(4,I) = UL*VL
		
C----------------------------------------------------------------------
C NAGV & WTV
	  !DETERMINE AG NUMBERS FOR INTERPOLATION & WEIGHTS
	  IF ((NIP2(NAGV_EVG(1,I)).GT.N_CTRD_AG) .OR. 
     * (NIM2(NAGV_EVG(1,I)).GT.N_CTRD_AG)) THEN !WEST OR EAST
          
		IF (COORD.EQ.'XY') THEN
		  AT1 = YV(NIP1(NAGV_EVG(1,I)))-YV(NAGV_EVG(1,I))
		  BT1 = XV(NAGV_EVG(1,I))-XV(NIP1(NAGV_EVG(1,I)))
		  CT1 = 0.
		  DIS1 = (YV(N_CTRD_AG+I)-YV(NAGV_EVG(1,I)))
     * -(-AT1/BT1*(XV(N_CTRD_AG+I)-XV(NAGV_EVG(1,I)))-CT1/BT1)
		  
		  AT2 = YV(NJP1(NAGV_EVG(1,I)))-YV(NAGV_EVG(1,I))
		  BT2 = XV(NAGV_EVG(1,I))-XV(NJP1(NAGV_EVG(1,I)))
		  CT2 = 0.
		  DIS2 = (XV(N_CTRD_AG+I)-XV(NAGV_EVG(1,I)))
     * -(-BT2/AT2*(YV(N_CTRD_AG+I)-YV(NAGV_EVG(1,I)))-CT2/AT2)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT1 = LATV(NIP1(NAGV_EVG(1,I)))-LATV(NAGV_EVG(1,I))
		  BT1 = LONV(NAGV_EVG(1,I))-LONV(NIP1(NAGV_EVG(1,I)))
		  CT1 = 0.
		  DIS1 = (LATV(N_CTRD_AG+I)-LATV(NAGV_EVG(1,I)))
     * -(-AT1/BT1*(LONV(N_CTRD_AG+I)-LONV(NAGV_EVG(1,I)))-CT1/BT1)
		  
		  AT2 = LATV(NJP1(NAGV_EVG(1,I)))-LATV(NAGV_EVG(1,I))
		  BT2 = LONV(NAGV_EVG(1,I))-LONV(NJP1(NAGV_EVG(1,I)))
		  CT2 = 0.
		  DIS2 = (LONV(N_CTRD_AG+I)-LONV(NAGV_EVG(1,I)))
     * -(-BT2/AT2*(LATV(N_CTRD_AG+I)-LATV(NAGV_EVG(1,I)))-CT2/AT2)
		ENDIF
		
		IF (DIS1.GE.0 .AND. DIS2.GE.0) THEN
		  NAGV_EVG(2,I) = NJP1(NAGV_EVG(1,I))
		  NAGV_EVG(3,I) = NIP1(NAGV_EVG(1,I))
		  NAGV_EVG(4,I) = NJP1(NIP1(NAGV_EVG(1,I)))
		ELSEIF (DIS1.GE.0 .AND. DIS2.LT.0) THEN
		  NAGV_EVG(2,I) = NJP1(NAGV_EVG(1,I))
		  NAGV_EVG(3,I) = NIM1(NAGV_EVG(1,I))
		  NAGV_EVG(4,I) = NJP1(NIM1(NAGV_EVG(1,I)))
		ELSEIF (DIS1.LT.0 .AND. DIS2.GE.0) THEN
		  NAGV_EVG(2,I) = NJM1(NAGV_EVG(1,I))
		  NAGV_EVG(3,I) = NIP1(NAGV_EVG(1,I))
		  NAGV_EVG(4,I) = NJM1(NIP1(NAGV_EVG(1,I)))
		ELSEIF (DIS1.LT.0 .AND. DIS2.LT.0) THEN
		  NAGV_EVG(2,I) = NJM1(NAGV_EVG(1,I))
		  NAGV_EVG(3,I) = NIM1(NAGV_EVG(1,I))
		  NAGV_EVG(4,I) = NJM1(NIM1(NAGV_EVG(1,I)))
		ENDIF
		
        ELSE !NORTH, SOUTH OR OTHERS
     
          IF (COORD.EQ.'XY') THEN
		  AT1 = YV(NIP1(NAGV_EVG(1,I)))-YV(NAGV_EVG(1,I))
		  BT1 = XV(NAGV_EVG(1,I))-XV(NIP1(NAGV_EVG(1,I)))
		  CT1 = 0.
		  DIS1 = (YV(N_CTRD_AG+I)-YV(NAGV_EVG(1,I)))
     * -(-AT1/BT1*(XV(N_CTRD_AG+I)-XV(NAGV_EVG(1,I)))-CT1/BT1)
		  
		  AT2 = YV(NJP1(NAGV_EVG(1,I)))-YV(NAGV_EVG(1,I))
		  BT2 = XV(NAGV_EVG(1,I))-XV(NJP1(NAGV_EVG(1,I)))
		  CT2 = 0.
		  DIS2 = (XV(N_CTRD_AG+I)-XV(NAGV_EVG(1,I)))
     * -(-BT2/AT2*(YV(N_CTRD_AG+I)-YV(NAGV_EVG(1,I)))-CT2/AT2)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT1 = LATV(NIP1(NAGV_EVG(1,I)))-LATV(NAGV_EVG(1,I))
		  BT1 = LONV(NAGV_EVG(1,I))-LONV(NIP1(NAGV_EVG(1,I)))
		  CT1 = 0.
		  DIS1 = (LATV(N_CTRD_AG+I)-LATV(NAGV_EVG(1,I)))
     * -(-AT1/BT1*(LONV(N_CTRD_AG+I)-LONV(NAGV_EVG(1,I)))-CT1/BT1)
		  
		  AT2 = LATV(NJP1(NAGV_EVG(1,I)))-LATV(NAGV_EVG(1,I))
		  BT2 = LONV(NAGV_EVG(1,I))-LONV(NJP1(NAGV_EVG(1,I)))
		  CT2 = 0.
		  DIS2 = (LONV(N_CTRD_AG+I)-LONV(NAGV_EVG(1,I)))
     * -(-BT2/AT2*(LATV(N_CTRD_AG+I)-LATV(NAGV_EVG(1,I)))-CT2/AT2)
		ENDIF
		
		IF (DIS1.GE.0 .AND. DIS2.GE.0) THEN
		  NAGV_EVG(2,I) = NIP1(NAGV_EVG(1,I))
		  NAGV_EVG(3,I) = NJP1(NAGV_EVG(1,I))
		  NAGV_EVG(4,I) = NIP1(NJP1(NAGV_EVG(1,I)))
		ELSEIF (DIS1.GE.0 .AND. DIS2.LT.0) THEN
		  NAGV_EVG(2,I) = NIM1(NAGV_EVG(1,I))
		  NAGV_EVG(3,I) = NJP1(NAGV_EVG(1,I))
		  NAGV_EVG(4,I) = NIM1(NJP1(NAGV_EVG(1,I)))
		ELSEIF (DIS1.LT.0 .AND. DIS2.GE.0) THEN
		  NAGV_EVG(2,I) = NIP1(NAGV_EVG(1,I))
		  NAGV_EVG(3,I) = NJM1(NAGV_EVG(1,I))
		  NAGV_EVG(4,I) = NIP1(NJM1(NAGV_EVG(1,I)))
		ELSEIF (DIS1.LT.0 .AND. DIS2.LT.0) THEN
		  NAGV_EVG(2,I) = NIM1(NAGV_EVG(1,I))
		  NAGV_EVG(3,I) = NJM1(NAGV_EVG(1,I))
		  NAGV_EVG(4,I) = NIM1(NJM1(NAGV_EVG(1,I)))
		ENDIF
	  
	  ENDIF
	  
	    DO J = 2,4
		  IF (NAGV_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XV(NAGV_EVG(J,I)),YV(NAGV_EVG(J,I))
     * ,XV(1:N_CTRD_AG),YV(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LONV(NAGV_EVG(J,I)),LATV(NAGV_EVG(J,I))
     * ,LONV(1:N_CTRD_AG),LATV(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGV_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		IF (COORD.EQ.'XY') THEN
		  EI = XV(NAGV_EVG(2,I))-XV(NAGV_EVG(1,I))
		  EJ = YV(NAGV_EVG(2,I))-YV(NAGV_EVG(1,I))
		  FI = XV(NAGV_EVG(3,I))-XV(NAGV_EVG(1,I))
		  FJ = YV(NAGV_EVG(3,I))-YV(NAGV_EVG(1,I))
		  GI = XV(NAGV_EVG(1,I))-XV(NAGV_EVG(2,I))
     * +XV(NAGV_EVG(4,I))-XV(NAGV_EVG(3,I))
		  GJ = YV(NAGV_EVG(1,I))-YV(NAGV_EVG(2,I))
     * +YV(NAGV_EVG(4,I))-YV(NAGV_EVG(3,I))
		  HI = XV(N_CTRD_AG+I)-XV(NAGV_EVG(1,I))
		  HJ = YV(N_CTRD_AG+I)-YV(NAGV_EVG(1,I))
		  
		ELSE 
		  EI = LONV(NAGV_EVG(2,I))-LONV(NAGV_EVG(1,I))
		  EJ = LATV(NAGV_EVG(2,I))-LATV(NAGV_EVG(1,I))
		  FI = LONV(NAGV_EVG(3,I))-LONV(NAGV_EVG(1,I))
		  FJ = LATV(NAGV_EVG(3,I))-LATV(NAGV_EVG(1,I))
		  GI = LONV(NAGV_EVG(1,I))-LONV(NAGV_EVG(2,I))
     * +LONV(NAGV_EVG(4,I))-LONV(NAGV_EVG(3,I))
		  GJ = LATV(NAGV_EVG(1,I))-LATV(NAGV_EVG(2,I))
     * +LATV(NAGV_EVG(4,I))-LATV(NAGV_EVG(3,I))
		  HI = LONV(N_CTRD_AG+I)-LONV(NAGV_EVG(1,I))
		  HJ = LATV(N_CTRD_AG+I)-LATV(NAGV_EVG(1,I))
		  
		ENDIF
		
		K2 = GI*FJ-GJ*FI
		K1 = EI*FJ-EJ*FI+HI*GJ-HJ*GI
		K0 = HI*EJ-HJ*EI
		
		IF (K2.EQ.0.) THEN
          VL = -K0/K1
          IF (VL.GE.-0.1 .AND. VL.LE.1.1) THEN
            CONTINUE
          ELSE
				PRINT*, 'WRONG VL = -K0/K1 IN WTV!'
				PRINT*, VL
				PRINT*, N_CTRD_AG+I
				PAUSE
				STOP
          ENDIF
        ELSE
          VL = (-K1+SQRT(K1**2-4*K0*K2))/(2*K2)
          IF (VL.GE.-0.1 .AND. VL.LE.1.1) THEN
				CONTINUE
			ELSE
				VL = (-K1-SQRT(K1**2-4*K0*K2))/(2*K2)
            IF (VL.GE.-0.1 .AND. VL.LE.1.1) THEN
              CONTINUE
            ELSE
					PRINT*, 'WRONG VL IN WTV!'
					PRINT*, VL
					PRINT*, N_CTRD_AG+I
					PAUSE
					STOP
            ENDIF
			ENDIF
        ENDIF
     
		UL = (HI-FI*VL)/(EI+GI*VL)
        IF (UL.GE.-0.1 .AND. UL.LE.1.1) THEN
          CONTINUE
        ELSE
          PRINT*,'WRONG UL IN WTV!' 
          PRINT*, UL,VL
          PRINT*, HI,FI,EI,GI
			PRINT*, N_CTRD_AG+I
			PAUSE
			STOP
        ENDIF
		
        WTV_EVG(1,I) = 1-UL-VL+UL*VL
        WTV_EVG(2,I) = UL*(1-VL)
        WTV_EVG(3,I) = VL*(1-UL)
        WTV_EVG(4,I) = UL*VL
		
	ENDDO
	  
	  
!**********************************************************************
! QUADRATIC INTERPOLATION
#elif defined INT_QDT

      !FIND IN WHICH AG IS THIS EVG
      DO I = 1,N_CTRD_EVG
	  DO J = 1,N_CTRD_AG
		DO K = 1,4
		  XB(K) = XNODE(K,J)
		  YB(K) = YNODE(K,J)
		ENDDO
		  
		IF (COORD.EQ.'XY') THEN
		  CALL INSIDE(XR(N_CTRD_AG+I),YR(N_CTRD_AG+I),XB,YB,NB,IND)
		ELSE
		  CALL INSIDE(LON(N_CTRD_AG+I),LAT(N_CTRD_AG+I),XB,YB,NB,IND)
		ENDIF
		  
		IF (IND.EQ.1) THEN
		  NAGC_EVG(1,I) = J
		  NAGU_EVG(1,I) = J
		  NAGV_EVG(1,I) = J
		  EXIT
		ELSEIF (J.EQ.N_CTRD_AG) THEN
		  PRINT*, 'CANNOT FIND MATCH GRID FOR EVG',N_CTRD_AG+I
		  PAUSE
		  STOP
		ENDIF
	  ENDDO
	ENDDO
	
	
	DO I = 1,N_CTRD_EVG
C----------------------------------------------------------------------
C NAGC & WTC
	  !DETERMINE AG NUMBERS FOR INTERPOLATION & WEIGHTS
	  IF ((NIP2(NAGC_EVG(1,I)).GT.N_CTRD_AG) .OR. 
     * (NIM2(NAGC_EVG(1,I)).GT.N_CTRD_AG)) THEN !WEST OR EAST
		
		IF (NIM1(NAGC_EVG(1,I)).EQ.N_CTRDP1) THEN
		  NAGC_EVG(:,I) = 1
		  NAGU_EVG(:,I) = 1
		  NAGV_EVG(:,I) = 1
		  CYCLE
		ENDIF
		
		IF (COORD.EQ.'XY') THEN
		  AT = YR(NIP1(NAGC_EVG(1,I)))-YR(NAGC_EVG(1,I))
		  BT = XR(NAGC_EVG(1,I))-XR(NIP1(NAGC_EVG(1,I)))
		  !CT = -XR(NAGC_EVG(1,I))*YR(NIP1(NAGC_EVG(1,I)))
    ! * +XR(NIP1(NAGC_EVG(1,I)))*YR(NAGC_EVG(1,I))
		  CT = 0.
		  
		  DIS = (YR(N_CTRD_AG+I)-YR(NAGC_EVG(1,I)))
     * -(-AT/BT*(XR(N_CTRD_AG+I)-XR(NAGC_EVG(1,I)))-CT/BT)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT = LAT(NIP1(NAGC_EVG(1,I)))-LAT(NAGC_EVG(1,I))
		  BT = LON(NAGC_EVG(1,I))-LON(NIP1(NAGC_EVG(1,I)))
		  !CT = -LON(NAGC_EVG(1,I))*LAT(NIP1(NAGC_EVG(1,I)))
    ! * +LON(NIP1(NAGC_EVG(1,I)))*LAT(NAGC_EVG(1,I))
		  CT = 0.
		  
		  DIS = (LAT(N_CTRD_AG+I)-LAT(NAGC_EVG(1,I)))
     * -(-AT/BT*(LON(N_CTRD_AG+I)-LON(NAGC_EVG(1,I)))-CT/BT)
          ENDIF
		
		
          IF (DIS.GE.0.) THEN
		  NAGC_EVG(2,I) = NJP1(NAGC_EVG(1,I))
		  NAGC_EVG(3,I) = NIP1(NAGC_EVG(1,I))
		  NAGC_EVG(4,I) = NJP1(NIP1(NAGC_EVG(1,I)))
		  NAGC_EVG(5,I) = NIM1(NAGC_EVG(1,I))
		  NAGC_EVG(6,I) = NJP1(NIM1(NAGC_EVG(1,I)))
		ELSE
		  NAGC_EVG(2,I) = NJM1(NAGC_EVG(1,I))
		  NAGC_EVG(3,I) = NIP1(NAGC_EVG(1,I))
		  NAGC_EVG(4,I) = NJM1(NIP1(NAGC_EVG(1,I)))
		  NAGC_EVG(5,I) = NIM1(NAGC_EVG(1,I))
		  NAGC_EVG(6,I) = NJM1(NIM1(NAGC_EVG(1,I)))
		ENDIF
		
		!PERPENDICULAR LINE
		IF (NAGC_EVG(3,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YR(NAGC_EVG(3,I))-YR(NAGC_EVG(1,I))
		    B1 = XR(NAGC_EVG(1,I))-XR(NAGC_EVG(3,I))
		    C1 = XR(N_CTRD_AG+I)*(-A1)+YR(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LAT(NAGC_EVG(3,I))-LAT(NAGC_EVG(1,I))
		    B1 = LON(NAGC_EVG(1,I))-LON(NAGC_EVG(3,I))
		    C1 = LON(N_CTRD_AG+I)*(-A1)+LAT(N_CTRD_AG+I)*(-B1)
		  ENDIF
		  
		ELSEIF (NAGC_EVG(5,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YR(NAGC_EVG(5,I))-YR(NAGC_EVG(1,I))
		    B1 = XR(NAGC_EVG(1,I))-XR(NAGC_EVG(5,I))
		    C1 = XR(N_CTRD_AG+I)*(-A1)+YR(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LAT(NAGC_EVG(5,I))-LAT(NAGC_EVG(1,I))
		    B1 = LON(NAGC_EVG(1,I))-LON(NAGC_EVG(5,I))
		    C1 = LON(N_CTRD_AG+I)*(-A1)+LAT(N_CTRD_AG+I)*(-B1)
		  ENDIF
		ELSE
		  PRINT*, 'ERROR IN PRE_IE.FOR!'
		  PAUSE
		  STOP
		ENDIF
		
		DO J = 2,6
		  IF (NAGC_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XR(NAGC_EVG(J,I)),YR(NAGC_EVG(J,I))
     * ,XR(1:N_CTRD_AG),YR(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LON(NAGC_EVG(J,I)),LAT(NAGC_EVG(J,I))
     * ,LON(1:N_CTRD_AG),LAT(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGC_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		!WEIGHTS FOR LINEAR INTERPOLATION
		IF (COORD.EQ.'XY') THEN
		  DIS1 = ABS((A1*XR(NAGC_EVG(2,I))+B1*YR(NAGC_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
		  !DIS2 = SQRT((XR(NAGC_EVG(2,I))-XR(NAGC_EVG(1,I)))**2
    ! * +(YR(NAGC_EVG(2,I))-YR(NAGC_EVG(1,I)))**2)
            DIS2 = ABS((A1*XR(NAGC_EVG(1,I))+B1*YR(NAGC_EVG(1,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(1,I) = DIS1/(DIS1+DIS2+EPSILON)
		  WTC_EVG(2,I) = 1-WTC_EVG(1,I)
		  
		  DIS1 = ABS((A1*XR(NAGC_EVG(4,I))+B1*YR(NAGC_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = ABS((A1*XR(NAGC_EVG(3,I))+B1*YR(NAGC_EVG(3,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(3,I) = DIS1/(DIS1+DIS2+EPSILON)
		  WTC_EVG(4,I) = 1-WTC_EVG(3,I)
		  
		  DIS1 = ABS((A1*XR(NAGC_EVG(6,I))+B1*YR(NAGC_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = ABS((A1*XR(NAGC_EVG(5,I))+B1*YR(NAGC_EVG(5,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(5,I) = DIS1/(DIS1+DIS2+EPSILON)
		  WTC_EVG(6,I) = 1-WTC_EVG(5,I)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  DIS1 = ABS((A1*LON(NAGC_EVG(2,I))+B1*LAT(NAGC_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = ABS((A1*LON(NAGC_EVG(1,I))+B1*LAT(NAGC_EVG(1,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(1,I) = DIS1/(DIS1+DIS2+EPSILON)
		  WTC_EVG(2,I) = 1-WTC_EVG(1,I)
		  
		  DIS1 = ABS((A1*LON(NAGC_EVG(4,I))+B1*LAT(NAGC_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = ABS((A1*LON(NAGC_EVG(3,I))+B1*LAT(NAGC_EVG(3,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(3,I) = DIS1/(DIS1+DIS2)
		  WTC_EVG(4,I) = 1-WTC_EVG(3,I)
		  
		  DIS1 = ABS((A1*LON(NAGC_EVG(6,I))+B1*LAT(NAGC_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = ABS((A1*LON(NAGC_EVG(5,I))+B1*LAT(NAGC_EVG(5,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(5,I) = DIS1/(DIS1+DIS2+EPSILON)
		  WTC_EVG(6,I) = 1-WTC_EVG(5,I)
		  
		ENDIF
		
		!WEIGHTS FOR QUADRATIC INTERPOLATION
		!IF (COORD.EQ.'XY') THEN
		!  !PERPENDICULAR LINE
		!  A2 = YR(NAGC_EVG(2,I))-YR(NAGC_EVG(1,I))
		!  B2 = XR(NAGC_EVG(1,I))-XR(NAGC_EVG(2,I))
		!  C2 = -XR(NAGC_EVG(1,I))*YR(NAGC_EVG(2,I))
  !   * +XR(NAGC_EVG(2,I))*YR(NAGC_EVG(1,I))
		!  
		!  !INTERSECTION POINT BETWEEN TANGENTIAL LINE AND PERPENDICULAR LINE
		!  XSEC = (B1*C2-B2*C1)/(A1*B2-A2*B1)
		!  YSEC = (A2*C1-A1*C2)/(A1*B2-A2*B1)
		!  
		!  !DISTANCE BETWEEN THE INTERSECTION AND EVG
		!  DIS = SIGN(SQRT((XR(N_CTRD_AG+I)-XSEC)**2+(YR(N_CTRD_AG+I)-YSEC)**2)
  !   * ,(XR(NAGC_EVG(3,I))-XR(NAGC_EVG(1,I)))
  !   * /(XR(N_CTRD_AG+I)-XSEC+0.01))
		!  !EPSILON
		!  EPS = DIS/H1(NAGC_EVG(1,I))
		!  !ALPHA
		!  ALPHA = 1./24.*((H1(N_CTRD_AG+I)/H1(NAGC_EVG(1,I)))**2-1.)
		!  
		!  !WEIGHTS
		!  WTC_EVG(7,I) = (1-EPS**2)-2*ALPHA
		!  WTC_EVG(8,I) = EPS*(EPS+1)/2.+ALPHA
		!  WTC_EVG(9,I) = EPS*(EPS-1)/2.+ALPHA
		!  
		!ELSEIF (COORD.EQ.'BL') THEN
			
		  RFN = NINT(H1(NAGC_EVG(1,I))/H1(N_CTRD_AG+I)) !REFINING MULTIPLE
		  IF (NIP2(NAGC_EVG(1,I)).GT.N_CTRD_AG) THEN !WEST
		    K = 0
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NIP1(II)
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGC_EVG(1,II-N_CTRD_AG).NE.NAGC_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
			K = RFN-K
			
		  ELSE !EAST
			K = 1
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NIM1(II)
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGC_EVG(1,II-N_CTRD_AG).NE.NAGC_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
		  ENDIF
		  
		  !EPSILON
		  EPS = (2*K-1.)/2.*H1(N_CTRD_AG+I)/H1(NAGC_EVG(1,I))-0.5
		  !ALPHA
		  ALPHA = 1./24.*((H1(N_CTRD_AG+I)/H1(NAGC_EVG(1,I)))**2-1.)
		  
		  !WEIGHTS
		  WTC_EVG(7,I) = (1-EPS**2)-2*ALPHA
		  WTC_EVG(8,I) = EPS*(EPS+1)/2.+ALPHA
		  WTC_EVG(9,I) = EPS*(EPS-1)/2.+ALPHA
		  
		!ENDIF
		
		  
C        ELSEIF ((NJM2(NAGC_EVG(1,I)).GT.N_CTRD_AG) .OR. 
C     * (NJP2(NAGC_EVG(1,I)).GT.N_CTRD_AG)) THEN !NORTH OR SOUTH
	  ELSE !NORTH OR SOUTH OR OTHERS
		  
		IF (COORD.EQ.'XY') THEN
		  AT = YR(NJP1(NAGC_EVG(1,I)))-YR(NAGC_EVG(1,I))
		  BT = XR(NAGC_EVG(1,I))-XR(NJP1(NAGC_EVG(1,I)))
		  !CT = -XR(NAGC_EVG(1,I))*YR(NJP1(NAGC_EVG(1,I)))
    ! * +XR(NJP1(NAGC_EVG(1,I)))*YR(NAGC_EVG(1,I))
		  CT = 0
		  
		  DIS = (XR(N_CTRD_AG+I)-XR(NAGC_EVG(1,I)))
     * -(-BT/AT*(YR(N_CTRD_AG+I)-YR(NAGC_EVG(1,I)))-CT/AT)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT = LAT(NJP1(NAGC_EVG(1,I)))-LAT(NAGC_EVG(1,I))
		  BT = LON(NAGC_EVG(1,I))-LON(NJP1(NAGC_EVG(1,I)))
		  !CT = -LON(NAGC_EVG(1,I))*LAT(NJP1(NAGC_EVG(1,I)))
    ! * +LON(NJP1(NAGC_EVG(1,I)))*LAT(NAGC_EVG(1,I))
		  CT = 0.
		  
		  DIS = (LON(N_CTRD_AG+I)-LON(NAGC_EVG(1,I)))
     * -(-BT/AT*(LAT(N_CTRD_AG+I)-LAT(NAGC_EVG(1,I)))-CT/AT)
		  
          ENDIF
		  
		IF (DIS.GE.0.) THEN
		  NAGC_EVG(2,I) = NIP1(NAGC_EVG(1,I))
		  NAGC_EVG(3,I) = NJP1(NAGC_EVG(1,I))
		  NAGC_EVG(4,I) = NIP1(NJP1(NAGC_EVG(1,I)))
		  NAGC_EVG(5,I) = NJM1(NAGC_EVG(1,I))
		  NAGC_EVG(6,I) = NIP1(NJM1(NAGC_EVG(1,I)))
		ELSE
		  NAGC_EVG(2,I) = NIM1(NAGC_EVG(1,I))
		  NAGC_EVG(3,I) = NJP1(NAGC_EVG(1,I))
		  NAGC_EVG(4,I) = NIM1(NJP1(NAGC_EVG(1,I)))
		  NAGC_EVG(5,I) = NJM1(NAGC_EVG(1,I))
		  NAGC_EVG(6,I) = NIM1(NJM1(NAGC_EVG(1,I)))
		ENDIF
		
		IF (NAGC_EVG(3,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YR(NAGC_EVG(3,I))-YR(NAGC_EVG(1,I))
		    B1 = XR(NAGC_EVG(1,I))-XR(NAGC_EVG(3,I))
		    C1 = XR(N_CTRD_AG+I)*(-A1)+YR(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LAT(NAGC_EVG(3,I))-LAT(NAGC_EVG(1,I))
		    B1 = LON(NAGC_EVG(1,I))-LON(NAGC_EVG(3,I))
		    C1 = LON(N_CTRD_AG+I)*(-A1)+LAT(N_CTRD_AG+I)*(-B1)
		  ENDIF
		  
		ELSEIF (NAGC_EVG(5,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YR(NAGC_EVG(5,I))-YR(NAGC_EVG(1,I))
		    B1 = XR(NAGC_EVG(1,I))-XR(NAGC_EVG(5,I))
		    C1 = XR(N_CTRD_AG+I)*(-A1)+YR(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LAT(NAGC_EVG(5,I))-LAT(NAGC_EVG(1,I))
		    B1 = LON(NAGC_EVG(1,I))-LON(NAGC_EVG(5,I))
		    C1 = LON(N_CTRD_AG+I)*(-A1)+LAT(N_CTRD_AG+I)*(-B1)
		  ENDIF
		ELSE
		  PRINT*, 'ERROR IN PRE_IE.FOR!'
		  PAUSE
		  STOP
		ENDIF
		
		DO J = 2,6
		  IF (NAGC_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XR(NAGC_EVG(J,I)),YR(NAGC_EVG(J,I))
     * ,XR(1:N_CTRD_AG),YR(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LON(NAGC_EVG(J,I)),LAT(NAGC_EVG(J,I))
     * ,LON(1:N_CTRD_AG),LAT(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGC_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		!WEIGHTS FOR LINEAR INTERPOLATION
		IF (COORD.EQ.'XY') THEN
		  DIS1 = ABS((A1*XR(NAGC_EVG(2,I))+B1*YR(NAGC_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
		  !DIS2 = SQRT((XR(NAGC_EVG(2,I))-XR(NAGC_EVG(1,I)))**2
    ! * +(YR(NAGC_EVG(2,I))-YR(NAGC_EVG(1,I)))**2)
            DIS2 = ABS((A1*XR(NAGC_EVG(1,I))+B1*YR(NAGC_EVG(1,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(1,I) = DIS1/(DIS1+DIS2)
		  WTC_EVG(2,I) = 1-WTC_EVG(1,I)
		  
		  DIS1 = ABS((A1*XR(NAGC_EVG(4,I))+B1*YR(NAGC_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = ABS((A1*XR(NAGC_EVG(3,I))+B1*YR(NAGC_EVG(3,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(3,I) = DIS1/(DIS1+DIS2+EPSILON)
		  WTC_EVG(4,I) = 1-WTC_EVG(3,I)
		  
		  DIS1 = ABS((A1*XR(NAGC_EVG(6,I))+B1*YR(NAGC_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = ABS((A1*XR(NAGC_EVG(5,I))+B1*YR(NAGC_EVG(5,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(5,I) = DIS1/(DIS1+DIS2+EPSILON)
		  WTC_EVG(6,I) = 1-WTC_EVG(5,I)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  DIS1 = ABS((A1*LON(NAGC_EVG(2,I))+B1*LAT(NAGC_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = ABS((A1*LON(NAGC_EVG(1,I))+B1*LAT(NAGC_EVG(1,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(1,I) = DIS1/(DIS1+DIS2+EPSILON)
		  WTC_EVG(2,I) = 1-WTC_EVG(1,I)
		  
		  DIS1 = ABS((A1*LON(NAGC_EVG(4,I))+B1*LAT(NAGC_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = ABS((A1*LON(NAGC_EVG(3,I))+B1*LAT(NAGC_EVG(3,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(3,I) = DIS1/(DIS1+DIS2+EPSILON)
		  WTC_EVG(4,I) = 1-WTC_EVG(3,I)
		  
		  DIS1 = ABS((A1*LON(NAGC_EVG(6,I))+B1*LAT(NAGC_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = ABS((A1*LON(NAGC_EVG(5,I))+B1*LAT(NAGC_EVG(5,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(5,I) = DIS1/(DIS1+DIS2+EPSILON)
		  WTC_EVG(6,I) = 1-WTC_EVG(5,I)
		  
		ENDIF
		
		!WEIGHTS FOR QUADRATIC INTERPOLATION
		  RFN = NINT(H2(NAGC_EVG(1,I))/H2(N_CTRD_AG+I)) !REFINING MULTIPLE
		  IF (NJM4(N_CTRD_AG+I).EQ.N_CTRDP1) THEN !SOUTH
		    K = 0
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NJP1(II)
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGC_EVG(1,II-N_CTRD_AG).NE.NAGC_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
			K = RFN-K
			
		  ELSE !NORTH
			K = 1
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NJM1(II)
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGC_EVG(1,II-N_CTRD_AG).NE.NAGC_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
		  ENDIF
		  
		  !EPSILON
		  EPS = (2*K-1)/2.*H2(N_CTRD_AG+I)/H2(NAGC_EVG(1,I))-0.5
		  !ALPHA
		  ALPHA = 1./24.*((H2(N_CTRD_AG+I)/H2(NAGC_EVG(1,I)))**2-1.)
		  
		  !WEIGHTS
		  WTC_EVG(7,I) = (1-EPS**2)-2*ALPHA
		  WTC_EVG(8,I) = EPS*(EPS+1)/2.+ALPHA
		  WTC_EVG(9,I) = EPS*(EPS-1)/2.+ALPHA

	  ENDIF

		
!----------------------------------------------------------------------
! NAGU & WTU
	  !DETERMINE AG NUMBERS FOR INTERPOLATION & WEIGHTS
	  IF ((NIP2(NAGU_EVG(1,I)).GT.N_CTRD_AG) .OR. 
     * (NIM2(NAGU_EVG(1,I)).GT.N_CTRD_AG)) THEN !WEST OR EAST
		
		IF (COORD.EQ.'XY') THEN
		  AT = YU(NIP1(NAGU_EVG(1,I)))-YU(NAGU_EVG(1,I))
		  BT = XU(NAGU_EVG(1,I))-XU(NIP1(NAGU_EVG(1,I)))
		  !CT = -XU(NAGU_EVG(1,I))*YU(NIP1(NAGU_EVG(1,I)))
    ! * +XU(NIP1(NAGU_EVG(1,I)))*YU(NAGU_EVG(1,I))
		  CT = 0.
		  
		  DIS = (YU(N_CTRD_AG+I)-YU(NAGU_EVG(1,I)))
     * -(-AT/BT*(XU(N_CTRD_AG+I)-XU(NAGU_EVG(1,I)))-CT/BT)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT = LATU(NIP1(NAGU_EVG(1,I)))-LATU(NAGU_EVG(1,I))
		  BT = LONU(NAGU_EVG(1,I))-LONU(NIP1(NAGU_EVG(1,I)))
		  !CT = -LONU(NAGU_EVG(1,I))*LATU(NIP1(NAGU_EVG(1,I)))
    ! * +LONU(NIP1(NAGU_EVG(1,I)))*LATU(NAGU_EVG(1,I))
		  CT = 0.
		  
		  DIS = (LATU(N_CTRD_AG+I)-LATU(NAGU_EVG(1,I)))
     * -(-AT/BT*(LONU(N_CTRD_AG+I)-LONU(NAGU_EVG(1,I)))-CT/BT)
          ENDIF
		  
		IF (DIS.GE.0.) THEN
		  NAGU_EVG(2,I) = NJP1(NAGU_EVG(1,I))
		  NAGU_EVG(3,I) = NIP1(NAGU_EVG(1,I))
		  NAGU_EVG(4,I) = NJP1(NIP1(NAGU_EVG(1,I)))
		  NAGU_EVG(5,I) = NIM1(NAGU_EVG(1,I))
		  NAGU_EVG(6,I) = NJP1(NIM1(NAGU_EVG(1,I)))
		ELSE
		  NAGU_EVG(2,I) = NJM1(NAGU_EVG(1,I))
		  NAGU_EVG(3,I) = NIP1(NAGU_EVG(1,I))
		  NAGU_EVG(4,I) = NJM1(NIP1(NAGU_EVG(1,I)))
		  NAGU_EVG(5,I) = NIM1(NAGU_EVG(1,I))
		  NAGU_EVG(6,I) = NJM1(NIM1(NAGU_EVG(1,I)))
		ENDIF
		
		!PERPENDICULAR LINE
		IF (NAGU_EVG(3,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YU(NAGU_EVG(3,I))-YU(NAGU_EVG(1,I))
		    B1 = XU(NAGU_EVG(1,I))-XU(NAGU_EVG(3,I))
		    C1 = XU(N_CTRD_AG+I)*(-A1)+YU(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATU(NAGU_EVG(3,I))-LATU(NAGU_EVG(1,I))
		    B1 = LONU(NAGU_EVG(1,I))-LONU(NAGU_EVG(3,I))
		    C1 = LONU(N_CTRD_AG+I)*(-A1)+LATU(N_CTRD_AG+I)*(-B1)
		  ENDIF
		  
		ELSEIF (NAGU_EVG(5,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YU(NAGU_EVG(5,I))-YU(NAGU_EVG(1,I))
		    B1 = XU(NAGU_EVG(1,I))-XU(NAGU_EVG(5,I))
		    C1 = XU(N_CTRD_AG+I)*(-A1)+YU(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATU(NAGU_EVG(5,I))-LATU(NAGU_EVG(1,I))
		    B1 = LONU(NAGU_EVG(1,I))-LONU(NAGU_EVG(5,I))
		    C1 = LONU(N_CTRD_AG+I)*(-A1)+LATU(N_CTRD_AG+I)*(-B1)
		  ENDIF
		ELSE
		  PRINT*, 'ERROR IN PRE_IE.FOR!'
		  PAUSE
		  STOP
		ENDIF
		
		DO J = 2,6
		  IF (NAGU_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XU(NAGU_EVG(J,I)),YU(NAGU_EVG(J,I))
     * ,XU(1:N_CTRD_AG),YU(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LONU(NAGU_EVG(J,I)),LATU(NAGU_EVG(J,I))
     * ,LONU(1:N_CTRD_AG),LATU(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGU_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		!WEIGHTS FOR LINEAR INTERPOLATION
		IF (COORD.EQ.'XY') THEN
		  DIS1 = ABS((A1*XU(NAGU_EVG(2,I))+B1*YU(NAGU_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XU(NAGU_EVG(2,I))-XU(NAGU_EVG(1,I)))**2
     * +(YU(NAGU_EVG(2,I))-YU(NAGU_EVG(1,I)))**2)
		  WTU_EVG(1,I) = DIS1/(DIS2+EPSILON)
		  WTU_EVG(2,I) = 1-WTU_EVG(1,I)
		  
		  DIS1 = ABS((A1*XU(NAGU_EVG(4,I))+B1*YU(NAGU_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XU(NAGU_EVG(4,I))-XU(NAGU_EVG(3,I)))**2
     * +(YU(NAGU_EVG(4,I))-YU(NAGU_EVG(3,I)))**2)
		  WTU_EVG(3,I) = DIS1/(DIS2+EPSILON)
		  WTU_EVG(4,I) = 1-WTU_EVG(3,I)
		  
		  DIS1 = ABS((A1*XU(NAGU_EVG(6,I))+B1*YU(NAGU_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XU(NAGU_EVG(6,I))-XU(NAGU_EVG(5,I)))**2
     * +(YU(NAGU_EVG(6,I))-YU(NAGU_EVG(5,I)))**2)
		  WTU_EVG(5,I) = DIS1/(DIS2+EPSILON)
		  WTU_EVG(6,I) = 1-WTU_EVG(5,I)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  DIS1 = ABS((A1*LONU(NAGU_EVG(2,I))+B1*LATU(NAGU_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONU(NAGU_EVG(2,I))-LONU(NAGU_EVG(1,I)))**2
     * +(LATU(NAGU_EVG(2,I))-LATU(NAGU_EVG(1,I)))**2)
		  WTU_EVG(1,I) = DIS1/(DIS2+EPSILON)
		  WTU_EVG(2,I) = 1-WTU_EVG(1,I)
		  
		  DIS1 = ABS((A1*LONU(NAGU_EVG(4,I))+B1*LATU(NAGU_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONU(NAGU_EVG(4,I))-LONU(NAGU_EVG(3,I)))**2
     * +(LATU(NAGU_EVG(4,I))-LATU(NAGU_EVG(3,I)))**2)
		  WTU_EVG(3,I) = DIS1/(DIS2+EPSILON)
		  WTU_EVG(4,I) = 1-WTU_EVG(3,I)
		  
		  DIS1 = ABS((A1*LONU(NAGU_EVG(6,I))+B1*LATU(NAGU_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONU(NAGU_EVG(6,I))-LONU(NAGU_EVG(5,I)))**2
     * +(LATU(NAGU_EVG(6,I))-LATU(NAGU_EVG(5,I)))**2)
		  WTU_EVG(5,I) = DIS1/(DIS2+EPSILON)
		  WTU_EVG(6,I) = 1-WTU_EVG(5,I)
		  
		ENDIF
		
		!WEIGHTS FOR QUADRATIC INTERPOLATION	
		  RFN = NINT(H1(NAGU_EVG(1,I))/H1(N_CTRD_AG+I)) !REFINING MULTIPLE
		  IF (NIP2(NAGU_EVG(1,I)).GT.N_CTRD_AG) THEN !WEST
		    K = 0
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NIP1(II)
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGU_EVG(1,II-N_CTRD_AG).NE.NAGU_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
			K = RFN-K
			
		  ELSE !EAST
			K = 1
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NIM1(II)
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGU_EVG(1,II-N_CTRD_AG).NE.NAGU_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
		  ENDIF
		  
		  !EPSILON
		  EPS = (K-1)*H1(N_CTRD_AG+I)/H1(NAGU_EVG(1,I))
		  !ALPHA
		  ALPHA = 1./24.*((H1(N_CTRD_AG+I)/H1(NAGU_EVG(1,I)))**2-1.)
		  
		  !WEIGHTS
		  WTU_EVG(7,I) = (1-EPS**2)-2*ALPHA
		  WTU_EVG(8,I) = EPS*(EPS+1)/2.+ALPHA
		  WTU_EVG(9,I) = EPS*(EPS-1)/2.+ALPHA
		  
C        ELSEIF ((NJM2(NAGU_EVG(1,I)).GT.N_CTRD_AG) .OR. 
C     * (NJP2(NAGU_EVG(1,I)).GT.N_CTRD_AG)) THEN !NORTH OR SOUTH
	  ELSE !NORTH OR SOUTH OR OTHERS
		  
		IF (COORD.EQ.'XY') THEN
		  AT = YU(NJP1(NAGU_EVG(1,I)))-YU(NAGU_EVG(1,I))
		  BT = XU(NAGU_EVG(1,I))-XU(NJP1(NAGU_EVG(1,I)))
		  !CT = -XU(NAGU_EVG(1,I))*YU(NJP1(NAGU_EVG(1,I)))
    ! * +XU(NJP1(NAGU_EVG(1,I)))*YU(NAGU_EVG(1,I))
		  CT = 0.
		  
		  DIS = (XU(N_CTRD_AG+I)-XU(NAGU_EVG(1,I)))
     * -(-BT/AT*(YU(N_CTRD_AG+I)-YU(NAGU_EVG(1,I)))-CT/AT)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT = LATU(NJP1(NAGU_EVG(1,I)))-LATU(NAGU_EVG(1,I))
		  BT = LONU(NAGU_EVG(1,I))-LONU(NJP1(NAGU_EVG(1,I)))
		  !CT = -LONU(NAGU_EVG(1,I))*LATU(NJP1(NAGU_EVG(1,I)))
    ! * +LONU(NJP1(NAGU_EVG(1,I)))*LATU(NAGU_EVG(1,I))
		  CT = 0.
		  
		  DIS = (LONU(N_CTRD_AG+I)-LONU(NAGU_EVG(1,I)))
     * -(-BT/AT*(LATU(N_CTRD_AG+I)-LATU(NAGU_EVG(1,I)))-CT/AT)
          ENDIF
		  
		IF (DIS.GE.0.) THEN
		  NAGU_EVG(2,I) = NIP1(NAGU_EVG(1,I))
		  NAGU_EVG(3,I) = NJP1(NAGU_EVG(1,I))
		  NAGU_EVG(4,I) = NIP1(NJP1(NAGU_EVG(1,I)))
		  NAGU_EVG(5,I) = NJM1(NAGU_EVG(1,I))
		  NAGU_EVG(6,I) = NIP1(NJM1(NAGU_EVG(1,I)))
		ELSE
		  NAGU_EVG(2,I) = NIM1(NAGU_EVG(1,I))
		  NAGU_EVG(3,I) = NJP1(NAGU_EVG(1,I))
		  NAGU_EVG(4,I) = NIM1(NJP1(NAGU_EVG(1,I)))
		  NAGU_EVG(5,I) = NJM1(NAGU_EVG(1,I))
		  NAGU_EVG(6,I) = NIM1(NJM1(NAGU_EVG(1,I)))
		ENDIF
		
		IF (NAGU_EVG(3,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YU(NAGU_EVG(3,I))-YU(NAGU_EVG(1,I))
		    B1 = XU(NAGU_EVG(1,I))-XU(NAGU_EVG(3,I))
		    C1 = XU(N_CTRD_AG+I)*(-A1)+YU(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATU(NAGU_EVG(3,I))-LATU(NAGU_EVG(1,I))
		    B1 = LONU(NAGU_EVG(1,I))-LONU(NAGU_EVG(3,I))
		    C1 = LONU(N_CTRD_AG+I)*(-A1)+LATU(N_CTRD_AG+I)*(-B1)
		  ENDIF
		  
		ELSEIF (NAGU_EVG(5,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YU(NAGU_EVG(5,I))-YU(NAGU_EVG(1,I))
		    B1 = XU(NAGU_EVG(1,I))-XU(NAGU_EVG(5,I))
		    C1 = XU(N_CTRD_AG+I)*(-A1)+YU(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATU(NAGU_EVG(5,I))-LATU(NAGU_EVG(1,I))
		    B1 = LONU(NAGU_EVG(1,I))-LONU(NAGU_EVG(5,I))
		    C1 = LONU(N_CTRD_AG+I)*(-A1)+LATU(N_CTRD_AG+I)*(-B1)
		  ENDIF
		ELSE
		  PRINT*, 'ERROR IN PRE_IE.FOR!'
		  PAUSE
		  STOP
		ENDIF
		
		DO J = 2,6
		  IF (NAGU_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XU(NAGU_EVG(J,I)),YU(NAGU_EVG(J,I))
     * ,XU(1:N_CTRD_AG),YU(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LONU(NAGU_EVG(J,I)),LATU(NAGU_EVG(J,I))
     * ,LONU(1:N_CTRD_AG),LATU(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGU_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		!WEIGHTS FOR LINEAR INTERPOLATION
		IF (COORD.EQ.'XY') THEN
		  DIS1 = ABS((A1*XU(NAGU_EVG(2,I))+B1*YU(NAGU_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XU(NAGU_EVG(2,I))-XU(NAGU_EVG(1,I)))**2
     * +(YU(NAGU_EVG(2,I))-YU(NAGU_EVG(1,I)))**2)
		  WTU_EVG(1,I) = DIS1/(DIS2+EPSILON)
		  WTU_EVG(2,I) = 1-WTU_EVG(1,I)
		  
		  DIS1 = ABS((A1*XU(NAGU_EVG(4,I))+B1*YU(NAGU_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XU(NAGU_EVG(4,I))-XU(NAGU_EVG(3,I)))**2
     * +(YU(NAGU_EVG(4,I))-YU(NAGU_EVG(3,I)))**2)
		  WTU_EVG(3,I) = DIS1/(DIS2+EPSILON)
		  WTU_EVG(4,I) = 1-WTU_EVG(3,I)
		  
		  DIS1 = ABS((A1*XU(NAGU_EVG(6,I))+B1*YU(NAGU_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XU(NAGU_EVG(6,I))-XU(NAGU_EVG(5,I)))**2
     * +(YU(NAGU_EVG(6,I))-YU(NAGU_EVG(5,I)))**2)
		  WTU_EVG(5,I) = DIS1/(DIS2+EPSILON)
		  WTU_EVG(6,I) = 1-WTU_EVG(5,I)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  DIS1 = ABS((A1*LONU(NAGU_EVG(2,I))+B1*LATU(NAGU_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONU(NAGU_EVG(2,I))-LONU(NAGU_EVG(1,I)))**2
     * +(LATU(NAGU_EVG(2,I))-LATU(NAGU_EVG(1,I)))**2)
		  WTU_EVG(1,I) = DIS1/(DIS2+EPSILON)
		  WTU_EVG(2,I) = 1-WTU_EVG(1,I)
		  
		  DIS1 = ABS((A1*LONU(NAGU_EVG(4,I))+B1*LATU(NAGU_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONU(NAGU_EVG(4,I))-LONU(NAGU_EVG(3,I)))**2
     * +(LATU(NAGU_EVG(4,I))-LATU(NAGU_EVG(3,I)))**2)
		  WTU_EVG(3,I) = DIS1/(DIS2+EPSILON)
		  WTU_EVG(4,I) = 1-WTU_EVG(3,I)
		  
		  DIS1 = ABS((A1*LONU(NAGU_EVG(6,I))+B1*LATU(NAGU_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONU(NAGU_EVG(6,I))-LONU(NAGU_EVG(5,I)))**2
     * +(LATU(NAGU_EVG(6,I))-LATU(NAGU_EVG(5,I)))**2)
		  WTU_EVG(5,I) = DIS1/(DIS2+EPSILON)
		  WTU_EVG(6,I) = 1-WTU_EVG(5,I)
		  
		ENDIF
		
		!WEIGHTS FOR QUADRATIC INTERPOLATION
		  RFN = NINT(H2(NAGU_EVG(1,I))/H2(N_CTRD_AG+I)) !REFINING MULTIPLE
		  IF (NJM4(N_CTRD_AG+I).EQ.N_CTRDP1) THEN !SOUTH
		    K = 0
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NJP1(II)
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGU_EVG(1,II-N_CTRD_AG).NE.NAGU_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
			K = RFN-K
			
		  ELSE !NORTH
			K = 1
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NJM1(II)
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGU_EVG(1,II-N_CTRD_AG).NE.NAGU_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
		  ENDIF
		  
		  !EPSILON
		  EPS = (2*K-1)/2.*H2(N_CTRD_AG+I)/H2(NAGU_EVG(1,I))-0.5
		  !ALPHA
		  ALPHA = 1./24.*((H2(N_CTRD_AG+I)/H2(NAGU_EVG(1,I)))**2-1.)
		  
		  !WEIGHTS
		  WTU_EVG(7,I) = (1-EPS**2)-2*ALPHA
		  WTU_EVG(8,I) = EPS*(EPS+1)/2.+ALPHA
		  WTU_EVG(9,I) = EPS*(EPS-1)/2.+ALPHA
		
		ENDIF
		
!----------------------------------------------------------------------
! NAGV & WTV
	  !DETERMINE AG NUMBERS FOR INTERPOLATION & WEIGHTS
	  IF ((NIP2(NAGV_EVG(1,I)).GT.N_CTRD_AG) .OR. 
     * (NIM2(NAGV_EVG(1,I)).GT.N_CTRD_AG)) THEN !WEST OR EAST
		
		IF (COORD.EQ.'XY') THEN
		  AT = YV(NIP1(NAGV_EVG(1,I)))-YV(NAGV_EVG(1,I))
		  BT = XV(NAGV_EVG(1,I))-XV(NIP1(NAGV_EVG(1,I)))
		  !CT = -XV(NAGV_EVG(1,I))*YV(NIP1(NAGV_EVG(1,I)))
    ! * +XV(NIP1(NAGV_EVG(1,I)))*YV(NAGV_EVG(1,I))
		  CT = 0.
		  
		  DIS = (YV(N_CTRD_AG+I)-YV(NAGV_EVG(1,I)))
     * -(-AT/BT*(XV(N_CTRD_AG+I)-XV(NAGV_EVG(1,I)))-CT/BT)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT = LATV(NIP1(NAGV_EVG(1,I)))-LATV(NAGV_EVG(1,I))
		  BT = LONV(NAGV_EVG(1,I))-LONV(NIP1(NAGV_EVG(1,I)))
		  !CT = -LONV(NAGV_EVG(1,I))*LATV(NIP1(NAGV_EVG(1,I)))
    ! * +LONV(NIP1(NAGV_EVG(1,I)))*LATV(NAGV_EVG(1,I))
		  CT = 0.
		  
		  DIS = (LATV(N_CTRD_AG+I)-LATV(NAGV_EVG(1,I)))
     * -(-AT/BT*(LONV(N_CTRD_AG+I)-LONV(NAGV_EVG(1,I)))-CT/BT)
          ENDIF
		  
		IF (DIS.GE.0.) THEN
		  NAGV_EVG(2,I) = NJP1(NAGV_EVG(1,I))
		  NAGV_EVG(3,I) = NIP1(NAGV_EVG(1,I))
		  NAGV_EVG(4,I) = NJP1(NIP1(NAGV_EVG(1,I)))
		  NAGV_EVG(5,I) = NIM1(NAGV_EVG(1,I))
		  NAGV_EVG(6,I) = NJP1(NIM1(NAGV_EVG(1,I)))
		ELSE
		  NAGV_EVG(2,I) = NJM1(NAGV_EVG(1,I))
		  NAGV_EVG(3,I) = NIP1(NAGV_EVG(1,I))
		  NAGV_EVG(4,I) = NJM1(NIP1(NAGV_EVG(1,I)))
		  NAGV_EVG(5,I) = NIM1(NAGV_EVG(1,I))
		  NAGV_EVG(6,I) = NJM1(NIM1(NAGV_EVG(1,I)))
		ENDIF
		
		!PERPENDICULAR LINE
		IF (NAGV_EVG(3,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YV(NAGV_EVG(3,I))-YV(NAGV_EVG(1,I))
		    B1 = XV(NAGV_EVG(1,I))-XV(NAGV_EVG(3,I))
		    C1 = XV(N_CTRD_AG+I)*(-A1)+YV(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATV(NAGV_EVG(3,I))-LATV(NAGV_EVG(1,I))
		    B1 = LONV(NAGV_EVG(1,I))-LONV(NAGV_EVG(3,I))
		    C1 = LONV(N_CTRD_AG+I)*(-A1)+LATV(N_CTRD_AG+I)*(-B1)
		  ENDIF
		  
		ELSEIF (NAGV_EVG(5,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YV(NAGV_EVG(5,I))-YV(NAGV_EVG(1,I))
		    B1 = XV(NAGV_EVG(1,I))-XV(NAGV_EVG(5,I))
		    C1 = XV(N_CTRD_AG+I)*(-A1)+YV(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATV(NAGV_EVG(5,I))-LATV(NAGV_EVG(1,I))
		    B1 = LONV(NAGV_EVG(1,I))-LONV(NAGV_EVG(5,I))
		    C1 = LONV(N_CTRD_AG+I)*(-A1)+LATV(N_CTRD_AG+I)*(-B1)
		  ENDIF
		ELSE
		  PRINT*, 'ERROR IN PRE_IE.FOR!'
		  PAUSE
		  STOP
		ENDIF
		
		DO J = 2,6
		  IF (NAGV_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XV(NAGV_EVG(J,I)),YV(NAGV_EVG(J,I))
     * ,XV(1:N_CTRD_AG),YV(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LONV(NAGV_EVG(J,I)),LATV(NAGV_EVG(J,I))
     * ,LONV(1:N_CTRD_AG),LATV(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGV_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		!WEIGHTS FOR LINEAR INTERPOLATION
		IF (COORD.EQ.'XY') THEN
		  DIS1 = ABS((A1*XV(NAGV_EVG(2,I))+B1*YV(NAGV_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XV(NAGV_EVG(2,I))-XV(NAGV_EVG(1,I)))**2
     * +(YV(NAGV_EVG(2,I))-YV(NAGV_EVG(1,I)))**2)
		  WTV_EVG(1,I) = DIS1/(DIS2+EPSILON)
		  WTV_EVG(2,I) = 1-WTV_EVG(1,I)
		  
		  DIS1 = ABS((A1*XV(NAGV_EVG(4,I))+B1*YV(NAGV_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XV(NAGV_EVG(4,I))-XV(NAGV_EVG(3,I)))**2
     * +(YV(NAGV_EVG(4,I))-YV(NAGV_EVG(3,I)))**2)
		  WTV_EVG(3,I) = DIS1/(DIS2+EPSILON)
		  WTV_EVG(4,I) = 1-WTV_EVG(3,I)
		  
		  DIS1 = ABS((A1*XV(NAGV_EVG(6,I))+B1*YV(NAGV_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XV(NAGV_EVG(6,I))-XV(NAGV_EVG(5,I)))**2
     * +(YV(NAGV_EVG(6,I))-YV(NAGV_EVG(5,I)))**2)
		  WTV_EVG(5,I) = DIS1/(DIS2+EPSILON)
		  WTV_EVG(6,I) = 1-WTV_EVG(5,I)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  DIS1 = ABS((A1*LONV(NAGV_EVG(2,I))+B1*LATV(NAGV_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONV(NAGV_EVG(2,I))-LONV(NAGV_EVG(1,I)))**2
     * +(LATV(NAGV_EVG(2,I))-LATV(NAGV_EVG(1,I)))**2)
		  WTV_EVG(1,I) = DIS1/(DIS2+EPSILON)
		  WTV_EVG(2,I) = 1-WTV_EVG(1,I)
		  
		  DIS1 = ABS((A1*LONV(NAGV_EVG(4,I))+B1*LATV(NAGV_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONV(NAGV_EVG(4,I))-LONV(NAGV_EVG(3,I)))**2
     * +(LATV(NAGV_EVG(4,I))-LATV(NAGV_EVG(3,I)))**2)
		  WTV_EVG(3,I) = DIS1/(DIS2+EPSILON)
		  WTV_EVG(4,I) = 1-WTV_EVG(3,I)
		  
		  DIS1 = ABS((A1*LONV(NAGV_EVG(6,I))+B1*LATV(NAGV_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONV(NAGV_EVG(6,I))-LONV(NAGV_EVG(5,I)))**2
     * +(LATV(NAGV_EVG(6,I))-LATV(NAGV_EVG(5,I)))**2)
		  WTV_EVG(5,I) = DIS1/(DIS2+EPSILON)
		  WTV_EVG(6,I) = 1-WTV_EVG(5,I)
		  
		ENDIF
		
		!WEIGHTS FOR QUADRATIC INTERPOLATION
		  RFN = NINT(H1(NAGV_EVG(1,I))/H1(N_CTRD_AG+I)) !REFINING MULTIPLE
		  IF (NIP2(NAGV_EVG(1,I)).GT.N_CTRD_AG) THEN !WEST
		    K = 0
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NIP1(II)
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGV_EVG(1,II-N_CTRD_AG).NE.NAGV_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
			K = RFN-K
			
		  ELSE !EAST
			K = 1
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NIM1(II)
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGV_EVG(1,II-N_CTRD_AG).NE.NAGV_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
		  ENDIF
		  
		  !EPSILON
		  EPS = (2*K-1)/2.*H1(N_CTRD_AG+I)/H1(NAGV_EVG(1,I))-0.5
		  !ALPHA
		  ALPHA = 1./24.*((H1(N_CTRD_AG+I)/H1(NAGV_EVG(1,I)))**2-1.)
		  
		  !WEIGHTS
		  WTV_EVG(7,I) = (1-EPS**2)-2*ALPHA
		  WTV_EVG(8,I) = EPS*(EPS+1)/2.+ALPHA
		  WTV_EVG(9,I) = EPS*(EPS-1)/2.+ALPHA
		  
C        ELSEIF ((NJM2(NAGV_EVG(1,I)).GT.N_CTRD_AG) .OR. 
C     * (NJP2(NAGV_EVG(1,I)).GT.N_CTRD_AG)) THEN !NORTH OR SOUTH
	  ELSE !NORTH OR SOUTH OR OTHERS
		  
		IF (COORD.EQ.'XY') THEN
		  AT = YV(NJP1(NAGV_EVG(1,I)))-YV(NAGV_EVG(1,I))
		  BT = XV(NAGV_EVG(1,I))-XV(NJP1(NAGV_EVG(1,I)))
		  !CT = -XV(NAGV_EVG(1,I))*YV(NJP1(NAGV_EVG(1,I)))
    ! * +XV(NJP1(NAGV_EVG(1,I)))*YV(NAGV_EVG(1,I))
		  CT = 0.
		  
		  DIS = (XV(N_CTRD_AG+I)-XV(NAGV_EVG(1,I)))
     * -(-BT/AT*(YV(N_CTRD_AG+I)-YV(NAGV_EVG(1,I)))-CT/AT)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT = LATV(NJP1(NAGV_EVG(1,I)))-LATV(NAGV_EVG(1,I))
		  BT = LONV(NAGV_EVG(1,I))-LONV(NJP1(NAGV_EVG(1,I)))
		  !CT = -LONV(NAGV_EVG(1,I))*LATV(NJP1(NAGV_EVG(1,I)))
    ! * +LONV(NJP1(NAGV_EVG(1,I)))*LATV(NAGV_EVG(1,I))
		  CT = 0.
		  
		  DIS = (LONV(N_CTRD_AG+I)-LONV(NAGV_EVG(1,I)))
     * -(-BT/AT*(LATV(N_CTRD_AG+I)-LATV(NAGV_EVG(1,I)))-CT/AT)
          ENDIF
		  
		IF (DIS.GE.0.) THEN
		  NAGV_EVG(2,I) = NIP1(NAGV_EVG(1,I))
		  NAGV_EVG(3,I) = NJP1(NAGV_EVG(1,I))
		  NAGV_EVG(4,I) = NIP1(NJP1(NAGV_EVG(1,I)))
		  NAGV_EVG(5,I) = NJM1(NAGV_EVG(1,I))
		  NAGV_EVG(6,I) = NIP1(NJM1(NAGV_EVG(1,I)))
		ELSE
		  NAGV_EVG(2,I) = NIM1(NAGV_EVG(1,I))
		  NAGV_EVG(3,I) = NJP1(NAGV_EVG(1,I))
		  NAGV_EVG(4,I) = NIM1(NJP1(NAGV_EVG(1,I)))
		  NAGV_EVG(5,I) = NJM1(NAGV_EVG(1,I))
		  NAGV_EVG(6,I) = NIM1(NJM1(NAGV_EVG(1,I)))
		ENDIF
		
		IF (NAGV_EVG(3,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YV(NAGV_EVG(3,I))-YV(NAGV_EVG(1,I))
		    B1 = XV(NAGV_EVG(1,I))-XV(NAGV_EVG(3,I))
		    C1 = XV(N_CTRD_AG+I)*(-A1)+YV(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATV(NAGV_EVG(3,I))-LATV(NAGV_EVG(1,I))
		    B1 = LONV(NAGV_EVG(1,I))-LONV(NAGV_EVG(3,I))
		    C1 = LONV(N_CTRD_AG+I)*(-A1)+LATV(N_CTRD_AG+I)*(-B1)
		  ENDIF
		  
		ELSEIF (NAGV_EVG(5,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YV(NAGV_EVG(5,I))-YV(NAGV_EVG(1,I))
		    B1 = XV(NAGV_EVG(1,I))-XV(NAGV_EVG(5,I))
		    C1 = XV(N_CTRD_AG+I)*(-A1)+YV(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATV(NAGV_EVG(5,I))-LATV(NAGV_EVG(1,I))
		    B1 = LONV(NAGV_EVG(1,I))-LONV(NAGV_EVG(5,I))
		    C1 = LONV(N_CTRD_AG+I)*(-A1)+LATV(N_CTRD_AG+I)*(-B1)
		  ENDIF
		ELSE
		  PRINT*, 'ERROR IN PRE_IE.FOR!'
		  PAUSE
		  STOP
		ENDIF
		
		DO J = 2,6
		  IF (NAGV_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XV(NAGV_EVG(J,I)),YV(NAGV_EVG(J,I))
     * ,XV(1:N_CTRD_AG),YV(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LONV(NAGV_EVG(J,I)),LATV(NAGV_EVG(J,I))
     * ,LONV(1:N_CTRD_AG),LATV(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGV_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		!WEIGHTS FOR LINEAR INTERPOLATION
		IF (COORD.EQ.'XY') THEN
		  DIS1 = ABS((A1*XV(NAGV_EVG(2,I))+B1*YV(NAGV_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XV(NAGV_EVG(2,I))-XV(NAGV_EVG(1,I)))**2
     * +(YV(NAGV_EVG(2,I))-YV(NAGV_EVG(1,I)))**2)
		  WTV_EVG(1,I) = DIS1/(DIS2+EPSILON)
		  WTV_EVG(2,I) = 1-WTV_EVG(1,I)
		  
		  DIS1 = ABS((A1*XV(NAGV_EVG(4,I))+B1*YV(NAGV_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XV(NAGV_EVG(4,I))-XV(NAGV_EVG(3,I)))**2
     * +(YV(NAGV_EVG(4,I))-YV(NAGV_EVG(3,I)))**2)
		  WTV_EVG(3,I) = DIS1/(DIS2+EPSILON)
		  WTV_EVG(4,I) = 1-WTV_EVG(3,I)
		  
		  DIS1 = ABS((A1*XV(NAGV_EVG(6,I))+B1*YV(NAGV_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XV(NAGV_EVG(6,I))-XV(NAGV_EVG(5,I)))**2
     * +(YV(NAGV_EVG(6,I))-YV(NAGV_EVG(5,I)))**2)
		  WTV_EVG(5,I) = DIS1/(DIS2+EPSILON)
		  WTV_EVG(6,I) = 1-WTV_EVG(5,I)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  DIS1 = ABS((A1*LONV(NAGV_EVG(2,I))+B1*LATV(NAGV_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONV(NAGV_EVG(2,I))-LONV(NAGV_EVG(1,I)))**2
     * +(LATV(NAGV_EVG(2,I))-LATV(NAGV_EVG(1,I)))**2)
		  WTV_EVG(1,I) = DIS1/(DIS2+EPSILON)
		  WTV_EVG(2,I) = 1-WTV_EVG(1,I)
		  
		  DIS1 = ABS((A1*LONV(NAGV_EVG(4,I))+B1*LATV(NAGV_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONV(NAGV_EVG(4,I))-LONV(NAGV_EVG(3,I)))**2
     * +(LATV(NAGV_EVG(4,I))-LATV(NAGV_EVG(3,I)))**2)
		  WTV_EVG(3,I) = DIS1/(DIS2+EPSILON)
		  WTV_EVG(4,I) = 1-WTV_EVG(3,I)
		  
		  DIS1 = ABS((A1*LONV(NAGV_EVG(6,I))+B1*LATV(NAGV_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONV(NAGV_EVG(6,I))-LONV(NAGV_EVG(5,I)))**2
     * +(LATV(NAGV_EVG(6,I))-LATV(NAGV_EVG(5,I)))**2)
		  WTV_EVG(5,I) = DIS1/(DIS2+EPSILON)
		  WTV_EVG(6,I) = 1-WTV_EVG(5,I)
		  
		ENDIF
		
		!WEIGHTS FOR QUADRATIC INTERPOLATION
		  RFN = NINT(H2(NAGV_EVG(1,I))/H2(N_CTRD_AG+I)) !REFINING MULTIPLE
		  IF (NJM4(N_CTRD_AG+I).EQ.N_CTRDP1) THEN !SOUTH
		    K = 0
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NJP1(II)
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGV_EVG(1,II-N_CTRD_AG).NE.NAGV_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
			K = RFN-K
			
		  ELSE !NORTH
			K = 1
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NJM1(II)
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGV_EVG(1,II-N_CTRD_AG).NE.NAGV_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
		  ENDIF
		  
		  !EPSILON
		  EPS = (K-1)*H2(N_CTRD_AG+I)/H2(NAGV_EVG(1,I))
		  !ALPHA
		  ALPHA = 1./24.*((H2(N_CTRD_AG+I)/H2(NAGV_EVG(1,I)))**2-1.)
		  
		  !WEIGHTS
		  WTV_EVG(7,I) = (1-EPS**2)-2*ALPHA
		  WTV_EVG(8,I) = EPS*(EPS+1)/2.+ALPHA
		  WTV_EVG(9,I) = EPS*(EPS-1)/2.+ALPHA
		
		ENDIF
	  
        ENDDO

        
!**********************************************************************
! HPI (HSIMT PARABOLIC INTERPOLATION)
#elif defined INT_HPI
      !FIND IN WHICH AG IS THIS EVG
      DO I = 1,N_CTRD_EVG
      DO J = 1,N_CTRD_AG
		DO K = 1,4
		    XB(K) = XNODE(K,J)
		    YB(K) = YNODE(K,J)
          ENDDO 
		IF (COORD.EQ.'XY') THEN
		    CALL INSIDE(XR(N_CTRD_AG+I),YR(N_CTRD_AG+I),XB,YB,NB,IND)
		ELSE
		    CALL INSIDE(LON(N_CTRD_AG+I),LAT(N_CTRD_AG+I),XB,YB,NB,IND)
		ENDIF
		IF (IND.EQ.1) THEN
		    NAGC_EVG(1,I) = J
		    NAGU_EVG(1,I) = J
		    NAGV_EVG(1,I) = J
		    EXIT
		ELSEIF (J.EQ.N_CTRD_AG) THEN
		    PRINT*, 'CANNOT FIND MATCH GRID FOR EVG',N_CTRD_AG+I
		    STOP
		ENDIF
      ENDDO
      ENDDO
        
      DO I = 1,N_CTRD_EVG
C----------------------------------------------------------------------
C NAGC & WTC
      !DETERMINE AG NUMBERS FOR INTERPOLATION & WEIGHTS
      IF ((NIP2(NAGC_EVG(1,I)).GT.N_CTRD_AG) .OR. 
     * (NIM2(NAGC_EVG(1,I)).GT.N_CTRD_AG)) THEN !WEST OR EAST
          
      IF (NIM1(NAGC_EVG(1,I)).EQ.N_CTRDP1) THEN
          NAGC_EVG(:,I) = 1
          NAGU_EVG(:,I) = 1
          NAGV_EVG(:,I) = 1
          CYCLE
      ENDIF
*Judgment point EVG(I) is above or below the  line berween 
*NAGC_EVG(1,I) and NIP1(...)      ---- Notes by MaRin
      IF (COORD.EQ.'XY') THEN
          AT = YR(NIP1(NAGC_EVG(1,I)))-YR(NAGC_EVG(1,I))
          BT = XR(NAGC_EVG(1,I))-XR(NIP1(NAGC_EVG(1,I)))
          CT = 0.
          DIS = (YR(N_CTRD_AG+I)-YR(NAGC_EVG(1,I)))
     *    -(-AT/BT*(XR(N_CTRD_AG+I)-XR(NAGC_EVG(1,I)))-CT/BT)
      ELSEIF (COORD.EQ.'BL') THEN
          AT = LAT(NIP1(NAGC_EVG(1,I)))-LAT(NAGC_EVG(1,I))
          BT = LON(NAGC_EVG(1,I))-LON(NIP1(NAGC_EVG(1,I)))
          CT = 0. 
          DIS = (LAT(N_CTRD_AG+I)-LAT(NAGC_EVG(1,I)))
     *    -(-AT/BT*(LON(N_CTRD_AG+I)-LON(NAGC_EVG(1,I)))-CT/BT)
      ENDIF
      IF (DIS.GE.0.) THEN
          NAGC_EVG(2,I) = NJP1(NAGC_EVG(1,I))
          NAGC_EVG(3,I) = NIM2(NAGC_EVG(1,I))
          NAGC_EVG(4,I) = FindCtrd(NAGC_EVG(1,I),-2,1)
          NAGC_EVG(5,I) = NIM1(NAGC_EVG(1,I))
          NAGC_EVG(6,I) = FindCtrd(NAGC_EVG(1,I),-1,1)
          NAGC_EVG(7,I) = NIP1(NAGC_EVG(1,I))
          NAGC_EVG(8,I) = FindCtrd(NAGC_EVG(1,I),1,1)
      ELSE
          NAGC_EVG(2,I) = NJM1(NAGC_EVG(1,I))
          NAGC_EVG(3,I) = NIM2(NAGC_EVG(1,I))
          NAGC_EVG(4,I) = FindCtrd(NAGC_EVG(1,I),-2,-1)
          NAGC_EVG(5,I) = NIM1(NAGC_EVG(1,I))
          NAGC_EVG(6,I) = FindCtrd(NAGC_EVG(1,I),-1,-1)
          NAGC_EVG(7,I) = NIP1(NAGC_EVG(1,I))
          NAGC_EVG(8,I) = FindCtrd(NAGC_EVG(1,I),1,-1)
      ENDIF
        !PERPENDICULAR LINE
      IF (NAGC_EVG(5,I).LE.N_CTRD_AG) THEN
          IF (COORD.EQ.'XY') THEN
		    A1 = YR(NAGC_EVG(5,I))-YR(NAGC_EVG(1,I))
		    B1 = XR(NAGC_EVG(1,I))-XR(NAGC_EVG(5,I))
		    C1 = XR(N_CTRD_AG+I)*(-A1)+YR(N_CTRD_AG+I)*(-B1)
          ELSEIF (COORD.EQ.'BL') THEN
              A1 = LAT(NAGC_EVG(5,I))-LAT(NAGC_EVG(1,I))
		    B1 = LON(NAGC_EVG(1,I))-LON(NAGC_EVG(5,I))
		    C1 = LON(N_CTRD_AG+I)*(-A1)+LAT(N_CTRD_AG+I)*(-B1)
          ENDIF
      ELSEIF (NAGC_EVG(7,I).LE.N_CTRD_AG) THEN
          IF (COORD.EQ.'XY') THEN
		    A1 = YR(NAGC_EVG(7,I))-YR(NAGC_EVG(1,I))
		    B1 = XR(NAGC_EVG(1,I))-XR(NAGC_EVG(7,I))
		    C1 = XR(N_CTRD_AG+I)*(-A1)+YR(N_CTRD_AG+I)*(-B1)
          ELSEIF (COORD.EQ.'BL') THEN
              A1 = LAT(NAGC_EVG(7,I))-LAT(NAGC_EVG(1,I))
		    B1 = LON(NAGC_EVG(1,I))-LON(NAGC_EVG(7,I))
		    C1 = LON(N_CTRD_AG+I)*(-A1)+LAT(N_CTRD_AG+I)*(-B1)
          ENDIF
      ELSE
          PRINT*, 'ERROR IN PRE_IE.FOR!'
          PAUSE
          STOP
      ENDIF
      DO J = 2,8
      IF (NAGC_EVG(J,I).GT.N_CTRD_AG) THEN
      IF (COORD.EQ.'XY') THEN
          CALL FIND_NEAREST_AG(XR(NAGC_EVG(J,I)),YR(NAGC_EVG(J,I))
     *    ,XR(1:N_CTRD_AG),YR(1:N_CTRD_AG),N_AG)
      ELSEIF (COORD.EQ.'BL') THEN
          CALL FIND_NEAREST_AG(LON(NAGC_EVG(J,I)),LAT(NAGC_EVG(J,I))
     *    ,LON(1:N_CTRD_AG),LAT(1:N_CTRD_AG),N_AG)
      ENDIF
      NAGC_EVG(J,I) = N_AG
      ENDIF
      ENDDO
      !WEIGHTS FOR LINEAR INTERPOLATION
      IF (COORD.EQ.'XY') THEN
          DO J = 1,4
          DIS1 = 
     *    ABS((A1*XR(NAGC_EVG(2*J,I))+B1*YR(NAGC_EVG(2*J,I))+C1)
     *    /SQRT(A1**2+B1**2))
      
          DIS2 = 
     *    ABS((A1*XR(NAGC_EVG(2*J-1,I))+B1*YR(NAGC_EVG(2*J-1,I))+C1)
     *    /SQRT(A1**2+B1**2))
      
          WTC_EVG(2*J-1,I) = DIS1/(DIS1+DIS2+EPSILON)
          WTC_EVG(2*J,I) = 1-WTC_EVG(2*J-1,I)
          ENDDO
      ELSEIF (COORD.EQ.'BL') THEN
          DO J = 1,4
            DIS1 = 
     *    ABS((A1*LON(NAGC_EVG(2*J,I))+B1*LAT(NAGC_EVG(2*J,I))+C1)
     *    /SQRT(A1**2+B1**2))
          DIS2 = 
     *    ABS((A1*LON(NAGC_EVG(2*J-1,I))+B1*LAT(NAGC_EVG(2*J-1,I))+C1)
     *    /SQRT(A1**2+B1**2))
          WTC_EVG(2*J-1,I) = DIS1/(DIS1+DIS2+EPSILON)
          WTC_EVG(2*J,I) = 1-WTC_EVG(2*J-1,I)
          ENDDO
      ENDIF
      !DEFINE XI
      RFN = NINT(H1(NAGC_EVG(1,I))/H1(N_CTRD_AG+I)) !REFINING MULTIPLE
      IF (NIP2(NAGC_EVG(1,I)).GT.N_CTRD_AG) THEN !WEST
          K = 0
          II = N_CTRD_AG+I
          DO WHILE (.TRUE.)
              II = NIP1(II)
              IF (II.LE.N_CTRD_AG) THEN
                  EXIT
              ELSEIF (NAGC_EVG(1,II-N_CTRD_AG).NE.NAGC_EVG(1,I)) THEN
                  EXIT 
              ENDIF
              K = K+1
          ENDDO
          K = RFN-K
      ELSE !EAST
          K = 1
          II = N_CTRD_AG+I
          DO WHILE (.TRUE.)
              II = NIM1(II)
              IF (II.LE.N_CTRD_AG) THEN
                  EXIT
              ELSEIF (NAGC_EVG(1,II-N_CTRD_AG).NE.NAGC_EVG(1,I)) THEN
                  EXIT
              ENDIF
              K = K+1
          ENDDO
		ENDIF  
		!EPSILON
      EPS = (2*K-1.)/2.*H1(N_CTRD_AG+I)/H1(NAGC_EVG(1,I))-0.5
      !XI
      XIC(I) = EPS+0.5
      !WEIGHTS FOR UPWIND AEI
      WTC_EVG(9,I)   = 1.-EPS
      WTC_EVG(10,I)  = EPS
      WTC_EVG(11,I)  = 1.+EPS
      WTC_EVG(12,I)  = -EPS
      ELSE !NORTH OR SOUTH OR OTHERS

     
        IF (COORD.EQ.'XY') THEN
		  AT = YR(NJP1(NAGC_EVG(1,I)))-YR(NAGC_EVG(1,I))
		  BT = XR(NAGC_EVG(1,I))-XR(NJP1(NAGC_EVG(1,I)))
		  CT = 0
		  
		  DIS = (XR(N_CTRD_AG+I)-XR(NAGC_EVG(1,I)))
     * -(-BT/AT*(YR(N_CTRD_AG+I)-YR(NAGC_EVG(1,I)))-CT/AT)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT = LAT(NJP1(NAGC_EVG(1,I)))-LAT(NAGC_EVG(1,I))
		  BT = LON(NAGC_EVG(1,I))-LON(NJP1(NAGC_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (LON(N_CTRD_AG+I)-LON(NAGC_EVG(1,I)))
     * -(-BT/AT*(LAT(N_CTRD_AG+I)-LAT(NAGC_EVG(1,I)))-CT/AT)
		  
        ENDIF
        
        IF (DIS.GE.0.) THEN
		  NAGC_EVG(2,I) = NIP1(NAGC_EVG(1,I))
		  NAGC_EVG(3,I) = NJM2(NAGC_EVG(1,I))
          NAGC_EVG(4,I) = FindCtrd(NAGC_EVG(1,I),1,-2)
		  NAGC_EVG(5,I) = NJM1(NAGC_EVG(1,I))
          NAGC_EVG(6,I) = FindCtrd(NAGC_EVG(1,I),1,-1)
          NAGC_EVG(7,I) = NJP1(NAGC_EVG(1,I))
          NAGC_EVG(8,I) = FindCtrd(NAGC_EVG(1,I),1,1)
		ELSE
		  NAGC_EVG(2,I) = NIM1(NAGC_EVG(1,I))
		  NAGC_EVG(3,I) = NJM2(NAGC_EVG(1,I))
          NAGC_EVG(4,I) = FindCtrd(NAGC_EVG(1,I),-1,-2)
		  NAGC_EVG(5,I) = NJM1(NAGC_EVG(1,I))
          NAGC_EVG(6,I) = FindCtrd(NAGC_EVG(1,I),-1,-1)
          NAGC_EVG(7,I) = NJP1(NAGC_EVG(1,I))
          NAGC_EVG(8,I) = FindCtrd(NAGC_EVG(1,I),-1,1)
        ENDIF
        
        IF (NAGC_EVG(5,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YR(NAGC_EVG(5,I))-YR(NAGC_EVG(1,I))
		    B1 = XR(NAGC_EVG(1,I))-XR(NAGC_EVG(5,I))
		    C1 = XR(N_CTRD_AG+I)*(-A1)+YR(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LAT(NAGC_EVG(5,I))-LAT(NAGC_EVG(1,I))
		    B1 = LON(NAGC_EVG(1,I))-LON(NAGC_EVG(5,I))
		    C1 = LON(N_CTRD_AG+I)*(-A1)+LAT(N_CTRD_AG+I)*(-B1)
		  ENDIF
		  
		ELSEIF (NAGC_EVG(7,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YR(NAGC_EVG(7,I))-YR(NAGC_EVG(1,I))
		    B1 = XR(NAGC_EVG(1,I))-XR(NAGC_EVG(7,I))
		    C1 = XR(N_CTRD_AG+I)*(-A1)+YR(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LAT(NAGC_EVG(7,I))-LAT(NAGC_EVG(1,I))
		    B1 = LON(NAGC_EVG(1,I))-LON(NAGC_EVG(7,I))
		    C1 = LON(N_CTRD_AG+I)*(-A1)+LAT(N_CTRD_AG+I)*(-B1)
		  ENDIF
		ELSE
		  PRINT*, 'ERROR IN PRE_IE.FOR!'
		  PAUSE
		  STOP
		ENDIF
		
		DO J = 2,8
		  IF (NAGC_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XR(NAGC_EVG(J,I)),YR(NAGC_EVG(J,I))
     * ,XR(1:N_CTRD_AG),YR(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LON(NAGC_EVG(J,I)),LAT(NAGC_EVG(J,I))
     * ,LON(1:N_CTRD_AG),LAT(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGC_EVG(J,I) = N_AG
		  ENDIF
        ENDDO
        
        !WEIGHTS FOR LINEAR INTERPOLATION
		IF (COORD.EQ.'XY') THEN
          DO J = 1,4
            DIS1 = 
     * ABS((A1*XR(NAGC_EVG(2*J,I))+B1*YR(NAGC_EVG(2*J,I))+C1)
     * /SQRT(A1**2+B1**2))
            DIS2 = 
     * ABS((A1*XR(NAGC_EVG(2*J-1,I))+B1*YR(NAGC_EVG(2*J-1,I))+C1)
     * /SQRT(A1**2+B1**2))
				WTC_EVG(2*J-1,I) = DIS1/(DIS1+DIS2+EPSILON)
				WTC_EVG(2*J,I) = 1-WTC_EVG(2*J-1,I)
          ENDDO
		  
        ELSEIF (COORD.EQ.'BL') THEN
          DO J = 1,4
            DIS1 = 
     * ABS((A1*LON(NAGC_EVG(2*J,I))+B1*LAT(NAGC_EVG(2*J,I))+C1)
     * /SQRT(A1**2+B1**2))
				DIS2 = 
     * ABS((A1*LON(NAGC_EVG(2*J-1,I))+B1*LAT(NAGC_EVG(2*J-1,I))+C1)
     * /SQRT(A1**2+B1**2))
				WTC_EVG(2*J-1,I) = DIS1/(DIS1+DIS2+EPSILON)
				WTC_EVG(2*J,I) = 1-WTC_EVG(2*J-1,I)
          ENDDO
          
        ENDIF
        
        !DEFINE XI
        RFN = NINT(H2(NAGC_EVG(1,I))/H2(N_CTRD_AG+I)) !REFINING MULTIPLE
		IF (NJM4(N_CTRD_AG+I).EQ.N_CTRDP1) THEN !SOUTH
		  K = 0
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
				II = NJP1(II)
				IF (II.LE.N_CTRD_AG) THEN
					EXIT
				ELSEIF (NAGC_EVG(1,II-N_CTRD_AG).NE.NAGC_EVG(1,I)) THEN
					EXIT
				ENDIF
			  
				K = K+1
			ENDDO
			
			K = RFN-K
			
		ELSE !NORTH
			K = 1
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
				II = NJM1(II)
				IF (II.LE.N_CTRD_AG) THEN
				EXIT
				ELSEIF (NAGC_EVG(1,II-N_CTRD_AG).NE.NAGC_EVG(1,I)) THEN
				EXIT
				ENDIF
			  
				K = K+1
			ENDDO
			
		ENDIF
		  
		!EPSILON
		EPS = (2*K-1)/2.*H2(N_CTRD_AG+I)/H2(NAGC_EVG(1,I))-0.5
        !XIC
        XIC(I) = EPS+0.5
        !WEIGHTS FOR UPWIND UAEI
        WTC_EVG(9,I)   = 1.-EPS
		WTC_EVG(10,I)  = EPS
		WTC_EVG(11,I)  = 1.+EPS
		WTC_EVG(12,I)  = -EPS
        
      ENDIF  
      
C----------------------------------------------------------------------
C NAGU & WTU
      !DETERMINE AG NUMBERS FOR INTERPOLATION & WEIGHTS
	  IF ((NIP2(NAGU_EVG(1,I)).GT.N_CTRD_AG) .OR. 
     * (NIM2(NAGU_EVG(1,I)).GT.N_CTRD_AG)) THEN !WEST OR EAST
		
		IF (COORD.EQ.'XY') THEN
		  AT = YU(NIP1(NAGU_EVG(1,I)))-YU(NAGU_EVG(1,I))
		  BT = XU(NAGU_EVG(1,I))-XU(NIP1(NAGU_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (YU(N_CTRD_AG+I)-YU(NAGU_EVG(1,I)))
     * -(-AT/BT*(XU(N_CTRD_AG+I)-XU(NAGU_EVG(1,I)))-CT/BT)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT = LATU(NIP1(NAGU_EVG(1,I)))-LATU(NAGU_EVG(1,I))
		  BT = LONU(NAGU_EVG(1,I))-LONU(NIP1(NAGU_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (LATU(N_CTRD_AG+I)-LATU(NAGU_EVG(1,I)))
     * -(-AT/BT*(LONU(N_CTRD_AG+I)-LONU(NAGU_EVG(1,I)))-CT/BT)
          ENDIF
		  
		IF (DIS.GE.0.) THEN
          NAGU_EVG(2,I) = NJP1(NAGU_EVG(1,I))
          NAGU_EVG(3,I) = NIM2(NAGU_EVG(1,I))
          NAGU_EVG(4,I) = FindCtrd(NAGU_EVG(1,I),-2,1)
          NAGU_EVG(5,I) = NIM1(NAGU_EVG(1,I))
          NAGU_EVG(6,I) = FindCtrd(NAGU_EVG(1,I),-1,1)
          NAGU_EVG(7,I) = NIP1(NAGU_EVG(1,I))
          NAGU_EVG(8,I) = FindCtrd(NAGU_EVG(1,I),1,1)
        ELSE
          NAGU_EVG(2,I) = NJM1(NAGU_EVG(1,I))
          NAGU_EVG(3,I) = NIM2(NAGU_EVG(1,I))
          NAGU_EVG(4,I) = FindCtrd(NAGU_EVG(1,I),-2,-1)
          NAGU_EVG(5,I) = NIM1(NAGU_EVG(1,I))
          NAGU_EVG(6,I) = FindCtrd(NAGU_EVG(1,I),-1,-1)
          NAGU_EVG(7,I) = NIP1(NAGU_EVG(1,I))
          NAGU_EVG(8,I) = FindCtrd(NAGU_EVG(1,I),1,-1)
        ENDIF
		
		!PERPENDICULAR LINE
		IF (NAGU_EVG(5,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YU(NAGU_EVG(5,I))-YU(NAGU_EVG(1,I))
		    B1 = XU(NAGU_EVG(1,I))-XU(NAGU_EVG(5,I))
		    C1 = XU(N_CTRD_AG+I)*(-A1)+YU(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATU(NAGU_EVG(5,I))-LATU(NAGU_EVG(1,I))
		    B1 = LONU(NAGU_EVG(1,I))-LONU(NAGU_EVG(5,I))
		    C1 = LONU(N_CTRD_AG+I)*(-A1)+LATU(N_CTRD_AG+I)*(-B1)
		  ENDIF
		  
		ELSEIF (NAGU_EVG(7,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YU(NAGU_EVG(7,I))-YU(NAGU_EVG(1,I))
		    B1 = XU(NAGU_EVG(1,I))-XU(NAGU_EVG(7,I))
		    C1 = XU(N_CTRD_AG+I)*(-A1)+YU(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATU(NAGU_EVG(7,I))-LATU(NAGU_EVG(1,I))
		    B1 = LONU(NAGU_EVG(1,I))-LONU(NAGU_EVG(7,I))
		    C1 = LONU(N_CTRD_AG+I)*(-A1)+LATU(N_CTRD_AG+I)*(-B1)
		  ENDIF
		ELSE
		  PRINT*, 'ERROR IN PRE_IE.FOR!'
		  PAUSE
		  STOP
		ENDIF
		
		DO J = 2,8
		  IF (NAGU_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XU(NAGU_EVG(J,I)),YU(NAGU_EVG(J,I))
     * ,XU(1:N_CTRD_AG),YU(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LONU(NAGU_EVG(J,I)),LATU(NAGU_EVG(J,I))
     * ,LONU(1:N_CTRD_AG),LATU(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGU_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		!WEIGHTS FOR LINEAR INTERPOLATION
		IF (COORD.EQ.'XY') THEN
          DO J = 1,4
            DIS1 = 
     * ABS((A1*XU(NAGU_EVG(2*J,I))+B1*YU(NAGU_EVG(2*J,I))+C1)
     * /SQRT(A1**2+B1**2))
				DIS2 = SQRT((XU(NAGU_EVG(2*J,I))-XU(NAGU_EVG(2*J-1,I)))**2
     * +(YU(NAGU_EVG(2*J,I))-YU(NAGU_EVG(2*J-1,I)))**2)
				WTU_EVG(2*J-1,I) = DIS1/(DIS2+EPSILON)
				WTU_EVG(2*J,I) = 1-WTU_EVG(2*J-1,I)
          ENDDO
          
        ELSEIF (COORD.EQ.'BL') THEN
          DO J = 1,4
            DIS1 = 
     * ABS((A1*LONU(NAGU_EVG(2*J,I))+B1*LATU(NAGU_EVG(2*J,I))+C1)
     * /SQRT(A1**2+B1**2))
				DIS2 = SQRT((LONU(NAGU_EVG(2*J,I))-LONU(NAGU_EVG(2*J-1,I)))**2
     * +(LATU(NAGU_EVG(2*J,I))-LATU(NAGU_EVG(2*J-1,I)))**2)
				WTU_EVG(2*J-1,I) = DIS1/(DIS2+EPSILON)
				WTU_EVG(2*J,I) = 1-WTU_EVG(2*J-1,I)
          ENDDO
        ENDIF
        
        !DEFINE XI
        RFN = NINT(H1(NAGU_EVG(1,I))/H1(N_CTRD_AG+I)) !REFINING MULTIPLE
		IF (NIP2(NAGU_EVG(1,I)).GT.N_CTRD_AG) THEN !WEST
		  K = 0
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
				II = NIP1(II)
				IF (II.LE.N_CTRD_AG) THEN
					EXIT
				ELSEIF (NAGU_EVG(1,II-N_CTRD_AG).NE.NAGU_EVG(1,I)) THEN
					EXIT
				ENDIF
			  
				K = K+1
			ENDDO
			
			K = RFN-K
			
		ELSE !EAST
			K = 1
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
				II = NIM1(II)
				IF (II.LE.N_CTRD_AG) THEN
					EXIT
				ELSEIF (NAGU_EVG(1,II-N_CTRD_AG).NE.NAGU_EVG(1,I)) THEN
					EXIT
				ENDIF
			  
				K = K+1
			ENDDO
			
		ENDIF
		  
		!EPSILON
		EPS = (K-1)*H1(N_CTRD_AG+I)/H1(NAGU_EVG(1,I))
        !XIU
        XIU(I) = EPS+0.5
        !WEIGHTS FOR UPWIND AEI
        WTU_EVG(9,I)   = 1.-EPS
		WTU_EVG(10,I)  = EPS
		WTU_EVG(11,I)  = 1.+EPS
		WTU_EVG(12,I)  = -EPS
      
      ELSE !NORTH OR SOUTH OR OTHERS
		  
		IF (COORD.EQ.'XY') THEN
		  AT = YU(NJP1(NAGU_EVG(1,I)))-YU(NAGU_EVG(1,I))
		  BT = XU(NAGU_EVG(1,I))-XU(NJP1(NAGU_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (XU(N_CTRD_AG+I)-XU(NAGU_EVG(1,I)))
     * -(-BT/AT*(YU(N_CTRD_AG+I)-YU(NAGU_EVG(1,I)))-CT/AT)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT = LATU(NJP1(NAGU_EVG(1,I)))-LATU(NAGU_EVG(1,I))
		  BT = LONU(NAGU_EVG(1,I))-LONU(NJP1(NAGU_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (LONU(N_CTRD_AG+I)-LONU(NAGU_EVG(1,I)))
     * -(-BT/AT*(LATU(N_CTRD_AG+I)-LATU(NAGU_EVG(1,I)))-CT/AT)
        ENDIF
		  
		IF (DIS.GE.0.) THEN
		  NAGU_EVG(2,I) = NIP1(NAGU_EVG(1,I))
		  NAGU_EVG(3,I) = NJM2(NAGU_EVG(1,I))
          NAGU_EVG(4,I) = FindCtrd(NAGU_EVG(1,I),1,-2)
		  NAGU_EVG(5,I) = NJM1(NAGU_EVG(1,I))
          NAGU_EVG(6,I) = FindCtrd(NAGU_EVG(1,I),1,-1)
          NAGU_EVG(7,I) = NJP1(NAGU_EVG(1,I))
          NAGU_EVG(8,I) = FindCtrd(NAGU_EVG(1,I),1,1)
		ELSE
		  NAGU_EVG(2,I) = NIM1(NAGU_EVG(1,I))
		  NAGU_EVG(3,I) = NJM2(NAGU_EVG(1,I))
          NAGU_EVG(4,I) = FindCtrd(NAGU_EVG(1,I),-1,-2)
		  NAGU_EVG(5,I) = NJM1(NAGU_EVG(1,I))
          NAGU_EVG(6,I) = FindCtrd(NAGU_EVG(1,I),-1,-1)
          NAGU_EVG(7,I) = NJP1(NAGU_EVG(1,I))
          NAGU_EVG(8,I) = FindCtrd(NAGU_EVG(1,I),-1,1)
        ENDIF
		
		IF (NAGU_EVG(5,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YU(NAGU_EVG(5,I))-YU(NAGU_EVG(1,I))
		    B1 = XU(NAGU_EVG(1,I))-XU(NAGU_EVG(5,I))
		    C1 = XU(N_CTRD_AG+I)*(-A1)+YU(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATU(NAGU_EVG(5,I))-LATU(NAGU_EVG(1,I))
		    B1 = LONU(NAGU_EVG(1,I))-LONU(NAGU_EVG(5,I))
		    C1 = LONU(N_CTRD_AG+I)*(-A1)+LATU(N_CTRD_AG+I)*(-B1)
		  ENDIF
		  
		ELSEIF (NAGU_EVG(7,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YU(NAGU_EVG(7,I))-YU(NAGU_EVG(1,I))
		    B1 = XU(NAGU_EVG(1,I))-XU(NAGU_EVG(7,I))
		    C1 = XU(N_CTRD_AG+I)*(-A1)+YU(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATU(NAGU_EVG(7,I))-LATU(NAGU_EVG(1,I))
		    B1 = LONU(NAGU_EVG(1,I))-LONU(NAGU_EVG(7,I))
		    C1 = LONU(N_CTRD_AG+I)*(-A1)+LATU(N_CTRD_AG+I)*(-B1)
		  ENDIF
		ELSE
		  PRINT*, 'ERROR IN PRE_IE.FOR!'
		  PAUSE
		  STOP
		ENDIF
		
		DO J = 2,8
		  IF (NAGU_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XU(NAGU_EVG(J,I)),YU(NAGU_EVG(J,I))
     * ,XU(1:N_CTRD_AG),YU(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LONU(NAGU_EVG(J,I)),LATU(NAGU_EVG(J,I))
     * ,LONU(1:N_CTRD_AG),LATU(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGU_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		!WEIGHTS FOR LINEAR INTERPOLATION
		IF (COORD.EQ.'XY') THEN
          DO J = 1,4
            DIS1 = 
     * ABS((A1*XU(NAGU_EVG(2*J,I))+B1*YU(NAGU_EVG(2*J,I))+C1)
     * /SQRT(A1**2+B1**2))
				DIS2 = SQRT((XU(NAGU_EVG(2*J,I))-XU(NAGU_EVG(2*J-1,I)))**2
     * +(YU(NAGU_EVG(2*J,I))-YU(NAGU_EVG(2*J-1,I)))**2)
				WTU_EVG(2*J-1,I) = DIS1/(DIS2+EPSILON)
				WTU_EVG(2*J,I) = 1-WTU_EVG(2*J-1,I)
          ENDDO
		  
        ELSEIF (COORD.EQ.'BL') THEN
          DO J = 1,4
            DIS1 = 
     * ABS((A1*LONU(NAGU_EVG(2*J,I))+B1*LATU(NAGU_EVG(2*J,I))+C1)
     * /SQRT(A1**2+B1**2))
				DIS2 = 
     * SQRT((LONU(NAGU_EVG(2*J,I))-LONU(NAGU_EVG(2*J-1,I)))**2
     * +(LATU(NAGU_EVG(2*J,I))-LATU(NAGU_EVG(2*J-1,I)))**2)
				WTU_EVG(2*J-1,I) = DIS1/(DIS2+EPSILON)
				WTU_EVG(2*J,I) = 1-WTU_EVG(2*J-1,I)
          ENDDO
        ENDIF
        
        !DEFINE XI
        RFN = NINT(H2(NAGU_EVG(1,I))/H2(N_CTRD_AG+I)) !REFINING MULTIPLE
		IF (NJM4(N_CTRD_AG+I).EQ.N_CTRDP1) THEN !SOUTH
		  K = 0
		II = N_CTRD_AG+I
		DO WHILE (.TRUE.)
			II = NJP1(II)
			IF (II.LE.N_CTRD_AG) THEN
			EXIT
			ELSEIF (NAGU_EVG(1,II-N_CTRD_AG).NE.NAGU_EVG(1,I)) THEN
			EXIT
			ENDIF
			  
			K = K+1
		ENDDO
			
		K = RFN-K
			
		ELSE !NORTH
			K = 1
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
				II = NJM1(II)
				IF (II.LE.N_CTRD_AG) THEN
					EXIT
				ELSEIF (NAGU_EVG(1,II-N_CTRD_AG).NE.NAGU_EVG(1,I)) THEN
					EXIT
				ENDIF
			  
				K = K+1
			ENDDO
			
		ENDIF
		  
		!EPSILON
		EPS = (2*K-1)/2.*H2(N_CTRD_AG+I)/H2(NAGU_EVG(1,I))-0.5
        !XIU
        XIU(I) = EPS+0.5
        !WEIGHTS FOR UPWIND AEI
        WTU_EVG(9,I)   = 1.-EPS
		WTU_EVG(10,I)  = EPS
		WTU_EVG(11,I)  = 1.+EPS
		WTU_EVG(12,I)  = -EPS
      
      ENDIF
      
      
C----------------------------------------------------------------------
C NAGV & WTV
      
      !DETERMINE AG NUMBERS FOR INTERPOLATION & WEIGHTS
	  IF ((NIP2(NAGV_EVG(1,I)).GT.N_CTRD_AG) .OR. 
     * (NIM2(NAGV_EVG(1,I)).GT.N_CTRD_AG)) THEN !WEST OR EAST
		
		IF (COORD.EQ.'XY') THEN
		  AT = YV(NIP1(NAGV_EVG(1,I)))-YV(NAGV_EVG(1,I))
		  BT = XV(NAGV_EVG(1,I))-XV(NIP1(NAGV_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (YV(N_CTRD_AG+I)-YV(NAGV_EVG(1,I)))
     * -(-AT/BT*(XV(N_CTRD_AG+I)-XV(NAGV_EVG(1,I)))-CT/BT)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT = LATV(NIP1(NAGV_EVG(1,I)))-LATV(NAGV_EVG(1,I))
		  BT = LONV(NAGV_EVG(1,I))-LONV(NIP1(NAGV_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (LATV(N_CTRD_AG+I)-LATV(NAGV_EVG(1,I)))
     * -(-AT/BT*(LONV(N_CTRD_AG+I)-LONV(NAGV_EVG(1,I)))-CT/BT)
        ENDIF
		  
		IF (DIS.GE.0.) THEN
          NAGV_EVG(2,I) = NJP1(NAGV_EVG(1,I))
          NAGV_EVG(3,I) = NIM2(NAGV_EVG(1,I))
          NAGV_EVG(4,I) = FindCtrd(NAGV_EVG(1,I),-2,1)
          NAGV_EVG(5,I) = NIM1(NAGV_EVG(1,I))
          NAGV_EVG(6,I) = FindCtrd(NAGV_EVG(1,I),-1,1)
          NAGV_EVG(7,I) = NIP1(NAGV_EVG(1,I))
          NAGV_EVG(8,I) = FindCtrd(NAGV_EVG(1,I),1,1)
        ELSE
          NAGV_EVG(2,I) = NJM1(NAGV_EVG(1,I))
          NAGV_EVG(3,I) = NIM2(NAGV_EVG(1,I))
          NAGV_EVG(4,I) = FindCtrd(NAGV_EVG(1,I),-2,-1)
          NAGV_EVG(5,I) = NIM1(NAGV_EVG(1,I))
          NAGV_EVG(6,I) = FindCtrd(NAGV_EVG(1,I),-1,-1)
          NAGV_EVG(7,I) = NIP1(NAGV_EVG(1,I))
          NAGV_EVG(8,I) = FindCtrd(NAGV_EVG(1,I),1,-1)
        ENDIF
		
		!PERPENDICULAR LINE
		IF (NAGV_EVG(5,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YV(NAGV_EVG(5,I))-YV(NAGV_EVG(1,I))
		    B1 = XV(NAGV_EVG(1,I))-XV(NAGV_EVG(5,I))
		    C1 = XV(N_CTRD_AG+I)*(-A1)+YV(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATV(NAGV_EVG(5,I))-LATV(NAGV_EVG(1,I))
		    B1 = LONV(NAGV_EVG(1,I))-LONV(NAGV_EVG(5,I))
		    C1 = LONV(N_CTRD_AG+I)*(-A1)+LATV(N_CTRD_AG+I)*(-B1)
		  ENDIF
		  
		ELSEIF (NAGV_EVG(7,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YV(NAGV_EVG(7,I))-YV(NAGV_EVG(1,I))
		    B1 = XV(NAGV_EVG(1,I))-XV(NAGV_EVG(7,I))
		    C1 = XV(N_CTRD_AG+I)*(-A1)+YV(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATV(NAGV_EVG(7,I))-LATV(NAGV_EVG(1,I))
		    B1 = LONV(NAGV_EVG(1,I))-LONV(NAGV_EVG(7,I))
		    C1 = LONV(N_CTRD_AG+I)*(-A1)+LATV(N_CTRD_AG+I)*(-B1)
		  ENDIF
		ELSE
		  PRINT*, 'ERROR IN PRE_IE.FOR!'
		  PAUSE
		  STOP
		ENDIF
		
		DO J = 2,8
		  IF (NAGV_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XV(NAGV_EVG(J,I)),YV(NAGV_EVG(J,I))
     * ,XV(1:N_CTRD_AG),YV(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LONV(NAGV_EVG(J,I)),LATV(NAGV_EVG(J,I))
     * ,LONV(1:N_CTRD_AG),LATV(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGV_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		!WEIGHTS FOR LINEAR INTERPOLATION
		IF (COORD.EQ.'XY') THEN
          DO J = 1,4
            DIS1 = 
     * ABS((A1*XV(NAGV_EVG(2*J,I))+B1*YV(NAGV_EVG(2*J,I))+C1)
     * /SQRT(A1**2+B1**2))
				DIS2 = SQRT((XV(NAGV_EVG(2*J,I))-XV(NAGV_EVG(2*J-1,I)))**2
     * +(YV(NAGV_EVG(2*J,I))-YV(NAGV_EVG(2*J-1,I)))**2)
				WTV_EVG(2*J-1,I) = DIS1/(DIS2+EPSILON)
				WTV_EVG(2*J,I) = 1-WTV_EVG(2*J-1,I)
          ENDDO
          
        ELSEIF (COORD.EQ.'BL') THEN
          DO J = 1,4
            DIS1 = 
     * ABS((A1*LONV(NAGV_EVG(2*J,I))+B1*LATV(NAGV_EVG(2*J,I))+C1)
     * /SQRT(A1**2+B1**2))
				DIS2 = 
     * SQRT((LONV(NAGV_EVG(2*J,I))-LONV(NAGV_EVG(2*J-1,I)))**2
     * +(LATV(NAGV_EVG(2*J,I))-LATV(NAGV_EVG(2*J-1,I)))**2)
				WTV_EVG(2*J-1,I) = DIS1/(DIS2+EPSILON)
				WTV_EVG(2*J,I) = 1-WTV_EVG(2*J-1,I)
          ENDDO
        ENDIF
        
        !DEFINE XI
        RFN = NINT(H1(NAGV_EVG(1,I))/H1(N_CTRD_AG+I)) !REFINING MULTIPLE
		IF (NIP2(NAGV_EVG(1,I)).GT.N_CTRD_AG) THEN !WEST
		  K = 0
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
				II = NIP1(II)
				IF (II.LE.N_CTRD_AG) THEN
					EXIT
				ELSEIF (NAGV_EVG(1,II-N_CTRD_AG).NE.NAGV_EVG(1,I)) THEN
					EXIT
				ENDIF
			  
				K = K+1
			ENDDO
			
			K = RFN-K
			
		ELSE !EAST
			K = 1
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
				II = NIM1(II)
				IF (II.LE.N_CTRD_AG) THEN
					EXIT
				ELSEIF (NAGV_EVG(1,II-N_CTRD_AG).NE.NAGV_EVG(1,I)) THEN
					EXIT
				ENDIF
			  
				K = K+1
			ENDDO
			
		ENDIF
		  
		!EPSILON
		EPS = (2*K-1)/2.*H1(N_CTRD_AG+I)/H1(NAGV_EVG(1,I))-0.5
        !XIV
        XIV(I) = EPS+0.5
        !WEIGHTS FOR UPWIND AEI
        WTV_EVG(9,I)   = 1.-EPS
		WTV_EVG(10,I)  = EPS
		WTV_EVG(11,I)  = 1.+EPS
		WTV_EVG(12,I)  = -EPS
      
      ELSE !NORTH OR SOUTH OR OTHERS
		  
		IF (COORD.EQ.'XY') THEN
		  AT = YV(NJP1(NAGV_EVG(1,I)))-YV(NAGV_EVG(1,I))
		  BT = XV(NAGV_EVG(1,I))-XV(NJP1(NAGV_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (XV(N_CTRD_AG+I)-XV(NAGV_EVG(1,I)))
     * -(-BT/AT*(YV(N_CTRD_AG+I)-YV(NAGV_EVG(1,I)))-CT/AT)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT = LATV(NJP1(NAGV_EVG(1,I)))-LATV(NAGV_EVG(1,I))
		  BT = LONV(NAGV_EVG(1,I))-LONV(NJP1(NAGV_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (LONV(N_CTRD_AG+I)-LONV(NAGV_EVG(1,I)))
     * -(-BT/AT*(LATV(N_CTRD_AG+I)-LATV(NAGV_EVG(1,I)))-CT/AT)
          ENDIF
		  
		IF (DIS.GE.0.) THEN
		  NAGV_EVG(2,I) = NIP1(NAGV_EVG(1,I))
		  NAGV_EVG(3,I) = NJM2(NAGV_EVG(1,I))
          NAGV_EVG(4,I) = FindCtrd(NAGV_EVG(1,I),1,-2)
		  NAGV_EVG(5,I) = NJM1(NAGV_EVG(1,I))
          NAGV_EVG(6,I) = FindCtrd(NAGV_EVG(1,I),1,-1)
          NAGV_EVG(7,I) = NJP1(NAGV_EVG(1,I))
          NAGV_EVG(8,I) = FindCtrd(NAGV_EVG(1,I),1,1)
		ELSE
		  NAGV_EVG(2,I) = NIM1(NAGV_EVG(1,I))
		  NAGV_EVG(3,I) = NJM2(NAGV_EVG(1,I))
          NAGV_EVG(4,I) = FindCtrd(NAGV_EVG(1,I),-1,-2)
		  NAGV_EVG(5,I) = NJM1(NAGV_EVG(1,I))
          NAGV_EVG(6,I) = FindCtrd(NAGV_EVG(1,I),-1,-1)
          NAGV_EVG(7,I) = NJP1(NAGV_EVG(1,I))
          NAGV_EVG(8,I) = FindCtrd(NAGV_EVG(1,I),-1,1)
        ENDIF
		
		IF (NAGV_EVG(5,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YV(NAGV_EVG(5,I))-YV(NAGV_EVG(1,I))
		    B1 = XV(NAGV_EVG(1,I))-XV(NAGV_EVG(5,I))
		    C1 = XV(N_CTRD_AG+I)*(-A1)+YV(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATV(NAGV_EVG(5,I))-LATV(NAGV_EVG(1,I))
		    B1 = LONV(NAGV_EVG(1,I))-LONV(NAGV_EVG(5,I))
		    C1 = LONV(N_CTRD_AG+I)*(-A1)+LATV(N_CTRD_AG+I)*(-B1)
		  ENDIF
		  
		ELSEIF (NAGV_EVG(7,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YV(NAGV_EVG(7,I))-YV(NAGV_EVG(1,I))
		    B1 = XV(NAGV_EVG(1,I))-XV(NAGV_EVG(7,I))
		    C1 = XV(N_CTRD_AG+I)*(-A1)+YV(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATV(NAGV_EVG(7,I))-LATV(NAGV_EVG(1,I))
		    B1 = LONV(NAGV_EVG(1,I))-LONV(NAGV_EVG(7,I))
		    C1 = LONV(N_CTRD_AG+I)*(-A1)+LATV(N_CTRD_AG+I)*(-B1)
		  ENDIF
		ELSE
		  PRINT*, 'ERROR IN PRE_IE.FOR!'
		  PAUSE
		  STOP
		ENDIF
		
		DO J = 2,8
		  IF (NAGV_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XV(NAGV_EVG(J,I)),YV(NAGV_EVG(J,I))
     * ,XV(1:N_CTRD_AG),YV(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LONV(NAGV_EVG(J,I)),LATV(NAGV_EVG(J,I))
     * ,LONV(1:N_CTRD_AG),LATV(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGV_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		!WEIGHTS FOR LINEAR INTERPOLATION
		IF (COORD.EQ.'XY') THEN
          DO J = 1,4
            DIS1 = 
     * ABS((A1*XV(NAGV_EVG(2*J,I))+B1*YV(NAGV_EVG(2*J,I))+C1)
     * /SQRT(A1**2+B1**2))
				DIS2 = SQRT((XV(NAGV_EVG(2*J,I))-XV(NAGV_EVG(2*J-1,I)))**2
     * +(YV(NAGV_EVG(2*J,I))-YV(NAGV_EVG(2*J-1,I)))**2)
				WTV_EVG(2*J-1,I) = DIS1/(DIS2+EPSILON)
				WTV_EVG(2*J,I) = 1-WTV_EVG(2*J-1,I)
          ENDDO
		  
        ELSEIF (COORD.EQ.'BL') THEN
          DO J = 1,4
            DIS1 = 
     * ABS((A1*LONV(NAGV_EVG(2*J,I))+B1*LATV(NAGV_EVG(2*J,I))+C1)
     * /SQRT(A1**2+B1**2))
				DIS2 = 
     * SQRT((LONV(NAGV_EVG(2*J,I))-LONV(NAGV_EVG(2*J-1,I)))**2
     * +(LATV(NAGV_EVG(2*J,I))-LATV(NAGV_EVG(2*J-1,I)))**2)
				WTV_EVG(2*J-1,I) = DIS1/(DIS2+EPSILON)
				WTV_EVG(2*J,I) = 1-WTV_EVG(2*J-1,I)
          ENDDO
        ENDIF
        
        !DEFINE XI
        RFN = NINT(H2(NAGV_EVG(1,I))/H2(N_CTRD_AG+I)) !REFINING MULTIPLE
		IF (NJM4(N_CTRD_AG+I).EQ.N_CTRDP1) THEN !SOUTH
		  K = 0
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
				II = NJP1(II)
				IF (II.LE.N_CTRD_AG) THEN
					EXIT
				ELSEIF (NAGV_EVG(1,II-N_CTRD_AG).NE.NAGV_EVG(1,I)) THEN
					EXIT
				ENDIF
			  
				K = K+1
			ENDDO
			
			K = RFN-K
			
		ELSE !NORTH
			K = 1
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
				II = NJM1(II)
				IF (II.LE.N_CTRD_AG) THEN
					EXIT
				ELSEIF (NAGV_EVG(1,II-N_CTRD_AG).NE.NAGV_EVG(1,I)) THEN
					EXIT
				ENDIF
			  
				K = K+1
			ENDDO
			
		ENDIF
		  
		!EPSILON
		EPS = (K-1)*H2(N_CTRD_AG+I)/H2(NAGV_EVG(1,I))
        !XIV
        XIV(I) = EPS+0.5
        !WEIGHTS
		WTV_EVG(7,I)  = 1.-EPS
		WTV_EVG(8,I)  = EPS
		WTV_EVG(9,I)  = 1.+EPS
		WTV_EVG(10,I) = -EPS
        
      ENDIF
      
      ENDDO
     
!**********************************************************************
! UPWIND AEI (UPWIND ADVECTION-EQUIVALENT INTERPOLATION)
#elif defined INT_UAEI

      !FIND IN WHICH AG IS THIS EVG
      DO I = 1,N_CTRD_EVG
	  DO J = 1,N_CTRD_AG
		DO K = 1,4
		  XB(K) = XNODE(K,J)
		  YB(K) = YNODE(K,J)
		ENDDO
		  
		IF (COORD.EQ.'XY') THEN
		  CALL INSIDE(XR(N_CTRD_AG+I),YR(N_CTRD_AG+I),XB,YB,NB,IND)
		ELSE
		  CALL INSIDE(LON(N_CTRD_AG+I),LAT(N_CTRD_AG+I),XB,YB,NB,IND)
		ENDIF
		  
		IF (IND.EQ.1) THEN
		  NAGC_EVG(1,I) = J
		  NAGU_EVG(1,I) = J
		  NAGV_EVG(1,I) = J
		  EXIT
		ELSEIF (J.EQ.N_CTRD_AG) THEN
		  PRINT*, 'CANNOT FIND MATCH GRID FOR EVG',N_CTRD_AG+I
		  PAUSE
		  STOP
		ENDIF
	  ENDDO
	ENDDO
	  
	
      DO I = 1,N_CTRD_EVG
C----------------------------------------------------------------------
C NAGC & WTC
	  IF ((NIP2(NAGC_EVG(1,I)).GT.N_CTRD_AG) .OR. 
     * (NIM2(NAGC_EVG(1,I)).GT.N_CTRD_AG)) THEN !WEST OR EAST
		
		IF (NIM1(NAGC_EVG(1,I)).EQ.N_CTRDP1) THEN
		  NAGC_EVG(:,I) = 1
		  NAGU_EVG(:,I) = 1
		  NAGV_EVG(:,I) = 1
		  CYCLE
		ENDIF
		
		IF (COORD.EQ.'XY') THEN
		  AT = YR(NIP1(NAGC_EVG(1,I)))-YR(NAGC_EVG(1,I))
		  BT = XR(NAGC_EVG(1,I))-XR(NIP1(NAGC_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (YR(N_CTRD_AG+I)-YR(NAGC_EVG(1,I)))
     * -(-AT/BT*(XR(N_CTRD_AG+I)-XR(NAGC_EVG(1,I)))-CT/BT)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT = LAT(NIP1(NAGC_EVG(1,I)))-LAT(NAGC_EVG(1,I))
		  BT = LON(NAGC_EVG(1,I))-LON(NIP1(NAGC_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (LAT(N_CTRD_AG+I)-LAT(NAGC_EVG(1,I)))
     * -(-AT/BT*(LON(N_CTRD_AG+I)-LON(NAGC_EVG(1,I)))-CT/BT)
		ENDIF
		
		IF (DIS.GE.0.) THEN
		  NAGC_EVG(2,I) = NJP1(NAGC_EVG(1,I))
		  NAGC_EVG(3,I) = NIP1(NAGC_EVG(1,I))
		  NAGC_EVG(4,I) = NJP1(NIP1(NAGC_EVG(1,I)))
		  NAGC_EVG(5,I) = NIM1(NAGC_EVG(1,I))
		  NAGC_EVG(6,I) = NJP1(NIM1(NAGC_EVG(1,I)))
		ELSE
		  NAGC_EVG(2,I) = NJM1(NAGC_EVG(1,I))
		  NAGC_EVG(3,I) = NIP1(NAGC_EVG(1,I))
		  NAGC_EVG(4,I) = NJM1(NIP1(NAGC_EVG(1,I)))
		  NAGC_EVG(5,I) = NIM1(NAGC_EVG(1,I))
		  NAGC_EVG(6,I) = NJM1(NIM1(NAGC_EVG(1,I)))
		ENDIF
		
		!PERPENDICULAR LINE
		IF (NAGC_EVG(3,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YR(NAGC_EVG(3,I))-YR(NAGC_EVG(1,I))
		    B1 = XR(NAGC_EVG(1,I))-XR(NAGC_EVG(3,I))
		    C1 = XR(N_CTRD_AG+I)*(-A1)+YR(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			A1 = LAT(NAGC_EVG(3,I))-LAT(NAGC_EVG(1,I))
		    B1 = LON(NAGC_EVG(1,I))-LON(NAGC_EVG(3,I))
		    C1 = LON(N_CTRD_AG+I)*(-A1)+LAT(N_CTRD_AG+I)*(-B1)
		  ENDIF
		  
		ELSEIF (NAGC_EVG(5,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YR(NAGC_EVG(5,I))-YR(NAGC_EVG(1,I))
		    B1 = XR(NAGC_EVG(1,I))-XR(NAGC_EVG(5,I))
		    C1 = XR(N_CTRD_AG+I)*(-A1)+YR(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			A1 = LAT(NAGC_EVG(5,I))-LAT(NAGC_EVG(1,I))
		    B1 = LON(NAGC_EVG(1,I))-LON(NAGC_EVG(5,I))
		    C1 = LON(N_CTRD_AG+I)*(-A1)+LAT(N_CTRD_AG+I)*(-B1)
		  ENDIF
		ELSE
		  PRINT*, 'ERROR IN PRE_IE.FOR!'
		  PAUSE
		  STOP
		ENDIF
		
		DO J = 2,6
		  IF (NAGC_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XR(NAGC_EVG(J,I)),YR(NAGC_EVG(J,I))
     * ,XR(1:N_CTRD_AG),YR(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LON(NAGC_EVG(J,I)),LAT(NAGC_EVG(J,I))
     * ,LON(1:N_CTRD_AG),LAT(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGC_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		!WEIGHTS FOR FIRST LINEAR INTERPOLATION
		IF (COORD.EQ.'XY') THEN
		  DIS1 = ABS((A1*XR(NAGC_EVG(2,I))+B1*YR(NAGC_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
            DIS2 = ABS((A1*XR(NAGC_EVG(1,I))+B1*YR(NAGC_EVG(1,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(1,I) = DIS1/(DIS1+DIS2)
		  WTC_EVG(2,I) = 1-WTC_EVG(1,I)
		  
		  DIS1 = ABS((A1*XR(NAGC_EVG(4,I))+B1*YR(NAGC_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = ABS((A1*XR(NAGC_EVG(3,I))+B1*YR(NAGC_EVG(3,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(3,I) = DIS1/(DIS1+DIS2)
		  WTC_EVG(4,I) = 1-WTC_EVG(3,I)
		  
		  DIS1 = ABS((A1*XR(NAGC_EVG(6,I))+B1*YR(NAGC_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = ABS((A1*XR(NAGC_EVG(5,I))+B1*YR(NAGC_EVG(5,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(5,I) = DIS1/(DIS1+DIS2)
		  WTC_EVG(6,I) = 1-WTC_EVG(5,I)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  DIS1 = ABS((A1*LON(NAGC_EVG(2,I))+B1*LAT(NAGC_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = ABS((A1*LON(NAGC_EVG(1,I))+B1*LAT(NAGC_EVG(1,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(1,I) = DIS1/(DIS1+DIS2)
		  WTC_EVG(2,I) = 1-WTC_EVG(1,I)
		  
		  DIS1 = ABS((A1*LON(NAGC_EVG(4,I))+B1*LAT(NAGC_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = ABS((A1*LON(NAGC_EVG(3,I))+B1*LAT(NAGC_EVG(3,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(3,I) = DIS1/(DIS1+DIS2)
		  WTC_EVG(4,I) = 1-WTC_EVG(3,I)
		  
		  DIS1 = ABS((A1*LON(NAGC_EVG(6,I))+B1*LAT(NAGC_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = ABS((A1*LON(NAGC_EVG(5,I))+B1*LAT(NAGC_EVG(5,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(5,I) = DIS1/(DIS1+DIS2)
		  WTC_EVG(6,I) = 1-WTC_EVG(5,I)
		  
		ENDIF
		
		!WEIGHTS FOR SECOND LINEAR INTERPOLATION
		RFN = NINT(H1(NAGC_EVG(1,I))/H1(N_CTRD_AG+I)) !REFINING MULTIPLE
		  IF (NIP2(NAGC_EVG(1,I)).GT.N_CTRD_AG) THEN !WEST
		    K = 0
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NIP1(II)
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGC_EVG(1,II-N_CTRD_AG).NE.NAGC_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
			K = RFN-K
			
		  ELSE !EAST
			K = 1
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NIM1(II)
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGC_EVG(1,II-N_CTRD_AG).NE.NAGC_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
		  ENDIF
		  
		  !EPSILON
		  EPS = (2*K-1)/2.*H1(N_CTRD_AG+I)/H1(NAGC_EVG(1,I))-0.5
		  
		  !WEIGHTS
		  WTC_EVG(7,I)  = 1.-EPS
		  WTC_EVG(8,I)  = EPS
		  WTC_EVG(9,I)  = 1.+EPS
		  WTC_EVG(10,I) = -EPS
		
        ELSE !NORTH OR SOUTH OR OTHERS
          !PERPENDICULAR LINE
          IF (COORD.EQ.'XY') THEN
		  AT = YR(NJP1(NAGC_EVG(1,I)))-YR(NAGC_EVG(1,I))
		  BT = XR(NAGC_EVG(1,I))-XR(NJP1(NAGC_EVG(1,I)))
		  CT = 0
		  
		  DIS = (XR(N_CTRD_AG+I)-XR(NAGC_EVG(1,I)))
     * -(-BT/AT*(YR(N_CTRD_AG+I)-YR(NAGC_EVG(1,I)))-CT/AT)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT = LAT(NJP1(NAGC_EVG(1,I)))-LAT(NAGC_EVG(1,I))
		  BT = LON(NAGC_EVG(1,I))-LON(NJP1(NAGC_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (LON(N_CTRD_AG+I)-LON(NAGC_EVG(1,I)))
     * -(-BT/AT*(LAT(N_CTRD_AG+I)-LAT(NAGC_EVG(1,I)))-CT/AT)
          ENDIF
		  
		IF (DIS.GE.0.) THEN
		  NAGC_EVG(2,I) = NIP1(NAGC_EVG(1,I))
		  NAGC_EVG(3,I) = NJP1(NAGC_EVG(1,I))
		  NAGC_EVG(4,I) = NIP1(NJP1(NAGC_EVG(1,I)))
		  NAGC_EVG(5,I) = NJM1(NAGC_EVG(1,I))
		  NAGC_EVG(6,I) = NIP1(NJM1(NAGC_EVG(1,I)))
		ELSE
		  NAGC_EVG(2,I) = NIM1(NAGC_EVG(1,I))
		  NAGC_EVG(3,I) = NJP1(NAGC_EVG(1,I))
		  NAGC_EVG(4,I) = NIM1(NJP1(NAGC_EVG(1,I)))
		  NAGC_EVG(5,I) = NJM1(NAGC_EVG(1,I))
		  NAGC_EVG(6,I) = NIM1(NJM1(NAGC_EVG(1,I)))
		ENDIF
		
		IF (NAGC_EVG(3,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YR(NAGC_EVG(3,I))-YR(NAGC_EVG(1,I))
		    B1 = XR(NAGC_EVG(1,I))-XR(NAGC_EVG(3,I))
		    C1 = XR(N_CTRD_AG+I)*(-A1)+YR(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			A1 = LAT(NAGC_EVG(3,I))-LAT(NAGC_EVG(1,I))
		    B1 = LON(NAGC_EVG(1,I))-LON(NAGC_EVG(3,I))
		    C1 = LON(N_CTRD_AG+I)*(-A1)+LAT(N_CTRD_AG+I)*(-B1)
		  ENDIF
		  
		ELSEIF (NAGC_EVG(5,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YR(NAGC_EVG(5,I))-YR(NAGC_EVG(1,I))
		    B1 = XR(NAGC_EVG(1,I))-XR(NAGC_EVG(5,I))
		    C1 = XR(N_CTRD_AG+I)*(-A1)+YR(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			A1 = LAT(NAGC_EVG(5,I))-LAT(NAGC_EVG(1,I))
		    B1 = LON(NAGC_EVG(1,I))-LON(NAGC_EVG(5,I))
		    C1 = LON(N_CTRD_AG+I)*(-A1)+LAT(N_CTRD_AG+I)*(-B1)
		  ENDIF
		ELSE
		  PRINT*, 'ERROR IN PRE_IE.FOR!'
		  PAUSE
		  STOP
		ENDIF
		
		DO J = 2,6
		  IF (NAGC_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XR(NAGC_EVG(J,I)),YR(NAGC_EVG(J,I))
     * ,XR(1:N_CTRD_AG),YR(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LON(NAGC_EVG(J,I)),LAT(NAGC_EVG(J,I))
     * ,LON(1:N_CTRD_AG),LAT(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGC_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		!WEIGHTS FOR FIRST LINEAR INTERPOLATION
		IF (COORD.EQ.'XY') THEN
		  DIS1 = ABS((A1*XR(NAGC_EVG(2,I))+B1*YR(NAGC_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
            DIS2 = ABS((A1*XR(NAGC_EVG(1,I))+B1*YR(NAGC_EVG(1,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(1,I) = DIS1/(DIS1+DIS2)
		  WTC_EVG(2,I) = 1-WTC_EVG(1,I)
		  
		  DIS1 = ABS((A1*XR(NAGC_EVG(4,I))+B1*YR(NAGC_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = ABS((A1*XR(NAGC_EVG(3,I))+B1*YR(NAGC_EVG(3,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(3,I) = DIS1/(DIS1+DIS2)
		  WTC_EVG(4,I) = 1-WTC_EVG(3,I)
		  
		  DIS1 = ABS((A1*XR(NAGC_EVG(6,I))+B1*YR(NAGC_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = ABS((A1*XR(NAGC_EVG(5,I))+B1*YR(NAGC_EVG(5,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(5,I) = DIS1/(DIS1+DIS2)
		  WTC_EVG(6,I) = 1-WTC_EVG(5,I)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  DIS1 = ABS((A1*LON(NAGC_EVG(2,I))+B1*LAT(NAGC_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = ABS((A1*LON(NAGC_EVG(1,I))+B1*LAT(NAGC_EVG(1,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(1,I) = DIS1/(DIS1+DIS2)
		  WTC_EVG(2,I) = 1-WTC_EVG(1,I)
		  
		  DIS1 = ABS((A1*LON(NAGC_EVG(4,I))+B1*LAT(NAGC_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = ABS((A1*LON(NAGC_EVG(3,I))+B1*LAT(NAGC_EVG(3,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(3,I) = DIS1/(DIS1+DIS2)
		  WTC_EVG(4,I) = 1-WTC_EVG(3,I)
		  
		  DIS1 = ABS((A1*LON(NAGC_EVG(6,I))+B1*LAT(NAGC_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = ABS((A1*LON(NAGC_EVG(5,I))+B1*LAT(NAGC_EVG(5,I))+C1)
     * /SQRT(A1**2+B1**2))
		  WTC_EVG(5,I) = DIS1/(DIS1+DIS2)
		  WTC_EVG(6,I) = 1-WTC_EVG(5,I)
		  
		ENDIF
		
		!WEIGHTS FOR SECOND LINEAR INTERPOLATION
		  RFN = NINT(H2(NAGC_EVG(1,I))/H2(N_CTRD_AG+I)) !REFINING MULTIPLE
		  IF (NJM4(N_CTRD_AG+I).EQ.N_CTRDP1) THEN !SOUTH
		    K = 0
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NJP1(II)
			  
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGC_EVG(1,II-N_CTRD_AG).NE.NAGC_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
			K = RFN-K
			
		  ELSE !NORTH
			K = 1
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NJM1(II)
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGC_EVG(1,II-N_CTRD_AG).NE.NAGC_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
		  ENDIF
		  
		  !EPSILON
		  EPS = (2*K-1)/2.*H2(N_CTRD_AG+I)/H2(NAGC_EVG(1,I))-0.5
		  
		  !WEIGHTS
		  WTC_EVG(7,I)  = 1.-EPS
		  WTC_EVG(8,I)  = EPS
		  WTC_EVG(9,I)  = 1.+EPS
		  WTC_EVG(10,I) = -EPS
	  
		ENDIF

!----------------------------------------------------------------------
! NAGU & WTU
	  !DETERMINE AG NUMBERS FOR INTERPOLATION & WEIGHTS
	  IF ((NIP2(NAGU_EVG(1,I)).GT.N_CTRD_AG) .OR. 
     * (NIM2(NAGU_EVG(1,I)).GT.N_CTRD_AG)) THEN !WEST OR EAST
		
		IF (COORD.EQ.'XY') THEN
		  AT = YU(NIP1(NAGU_EVG(1,I)))-YU(NAGU_EVG(1,I))
		  BT = XU(NAGU_EVG(1,I))-XU(NIP1(NAGU_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (YU(N_CTRD_AG+I)-YU(NAGU_EVG(1,I)))
     * -(-AT/BT*(XU(N_CTRD_AG+I)-XU(NAGU_EVG(1,I)))-CT/BT)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT = LATU(NIP1(NAGU_EVG(1,I)))-LATU(NAGU_EVG(1,I))
		  BT = LONU(NAGU_EVG(1,I))-LONU(NIP1(NAGU_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (LATU(N_CTRD_AG+I)-LATU(NAGU_EVG(1,I)))
     * -(-AT/BT*(LONU(N_CTRD_AG+I)-LONU(NAGU_EVG(1,I)))-CT/BT)
          ENDIF
		  
		IF (DIS.GE.0.) THEN
		  NAGU_EVG(2,I) = NJP1(NAGU_EVG(1,I))
		  NAGU_EVG(3,I) = NIP1(NAGU_EVG(1,I))
		  NAGU_EVG(4,I) = NJP1(NIP1(NAGU_EVG(1,I)))
		  NAGU_EVG(5,I) = NIM1(NAGU_EVG(1,I))
		  NAGU_EVG(6,I) = NJP1(NIM1(NAGU_EVG(1,I)))
		ELSE
		  NAGU_EVG(2,I) = NJM1(NAGU_EVG(1,I))
		  NAGU_EVG(3,I) = NIP1(NAGU_EVG(1,I))
		  NAGU_EVG(4,I) = NJM1(NIP1(NAGU_EVG(1,I)))
		  NAGU_EVG(5,I) = NIM1(NAGU_EVG(1,I))
		  NAGU_EVG(6,I) = NJM1(NIM1(NAGU_EVG(1,I)))
		ENDIF
		
		!PERPENDICULAR LINE
		IF (NAGU_EVG(3,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YU(NAGU_EVG(3,I))-YU(NAGU_EVG(1,I))
		    B1 = XU(NAGU_EVG(1,I))-XU(NAGU_EVG(3,I))
		    C1 = XU(N_CTRD_AG+I)*(-A1)+YU(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			A1 = LATU(NAGU_EVG(3,I))-LATU(NAGU_EVG(1,I))
		    B1 = LONU(NAGU_EVG(1,I))-LONU(NAGU_EVG(3,I))
		    C1 = LONU(N_CTRD_AG+I)*(-A1)+LATU(N_CTRD_AG+I)*(-B1)
		  ENDIF
		  
		ELSEIF (NAGU_EVG(5,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YU(NAGU_EVG(5,I))-YU(NAGU_EVG(1,I))
		    B1 = XU(NAGU_EVG(1,I))-XU(NAGU_EVG(5,I))
		    C1 = XU(N_CTRD_AG+I)*(-A1)+YU(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			A1 = LATU(NAGU_EVG(5,I))-LATU(NAGU_EVG(1,I))
		    B1 = LONU(NAGU_EVG(1,I))-LONU(NAGU_EVG(5,I))
		    C1 = LONU(N_CTRD_AG+I)*(-A1)+LATU(N_CTRD_AG+I)*(-B1)
		  ENDIF
		ELSE
		  PRINT*, 'ERROR IN PRE_IE.FOR!'
		  PAUSE
		  STOP
		ENDIF
		
		DO J = 2,6
		  IF (NAGU_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XU(NAGU_EVG(J,I)),YU(NAGU_EVG(J,I))
     * ,XU(1:N_CTRD_AG),YU(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LONU(NAGU_EVG(J,I)),LATU(NAGU_EVG(J,I))
     * ,LONU(1:N_CTRD_AG),LATU(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGU_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		!WEIGHTS FOR FIRST LINEAR INTERPOLATION
		IF (COORD.EQ.'XY') THEN
		  DIS1 = ABS((A1*XU(NAGU_EVG(2,I))+B1*YU(NAGU_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XU(NAGU_EVG(2,I))-XU(NAGU_EVG(1,I)))**2
     * +(YU(NAGU_EVG(2,I))-YU(NAGU_EVG(1,I)))**2)
		  WTU_EVG(1,I) = DIS1/DIS2
		  WTU_EVG(2,I) = 1-WTU_EVG(1,I)
		  
		  DIS1 = ABS((A1*XU(NAGU_EVG(4,I))+B1*YU(NAGU_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XU(NAGU_EVG(4,I))-XU(NAGU_EVG(3,I)))**2
     * +(YU(NAGU_EVG(4,I))-YU(NAGU_EVG(3,I)))**2)
		  WTU_EVG(3,I) = DIS1/DIS2
		  WTU_EVG(4,I) = 1-WTU_EVG(3,I)
		  
		  DIS1 = ABS((A1*XU(NAGU_EVG(6,I))+B1*YU(NAGU_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XU(NAGU_EVG(6,I))-XU(NAGU_EVG(5,I)))**2
     * +(YU(NAGU_EVG(6,I))-YU(NAGU_EVG(5,I)))**2)
		  WTU_EVG(5,I) = DIS1/DIS2
		  WTU_EVG(6,I) = 1-WTU_EVG(5,I)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  DIS1 = ABS((A1*LONU(NAGU_EVG(2,I))+B1*LATU(NAGU_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONU(NAGU_EVG(2,I))-LONU(NAGU_EVG(1,I)))**2
     * +(LATU(NAGU_EVG(2,I))-LATU(NAGU_EVG(1,I)))**2)
		  WTU_EVG(1,I) = DIS1/DIS2
		  WTU_EVG(2,I) = 1-WTU_EVG(1,I)
		  
		  DIS1 = ABS((A1*LONU(NAGU_EVG(4,I))+B1*LATU(NAGU_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONU(NAGU_EVG(4,I))-LONU(NAGU_EVG(3,I)))**2
     * +(LATU(NAGU_EVG(4,I))-LATU(NAGU_EVG(3,I)))**2)
		  WTU_EVG(3,I) = DIS1/DIS2
		  WTU_EVG(4,I) = 1-WTU_EVG(3,I)
		  
		  DIS1 = ABS((A1*LONU(NAGU_EVG(6,I))+B1*LATU(NAGU_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONU(NAGU_EVG(6,I))-LONU(NAGU_EVG(5,I)))**2
     * +(LATU(NAGU_EVG(6,I))-LATU(NAGU_EVG(5,I)))**2)
		  WTU_EVG(5,I) = DIS1/DIS2
		  WTU_EVG(6,I) = 1-WTU_EVG(5,I)
		  
		ENDIF
		
		!WEIGHTS FOR SECOND LINEAR INTERPOLATION	
		  RFN = NINT(H1(NAGU_EVG(1,I))/H1(N_CTRD_AG+I)) !REFINING MULTIPLE
		  IF (NIP2(NAGU_EVG(1,I)).GT.N_CTRD_AG) THEN !WEST
		    K = 0
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NIP1(II)
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGU_EVG(1,II-N_CTRD_AG).NE.NAGU_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
			K = RFN-K
			
		  ELSE !EAST
			K = 1
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NIM1(II)
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGU_EVG(1,II-N_CTRD_AG).NE.NAGU_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
		  ENDIF
		  
		  !EPSILON
		  EPS = (K-1)*H1(N_CTRD_AG+I)/H1(NAGU_EVG(1,I))
		  
		  !WEIGHTS
		  WTU_EVG(7,I)  = 1.-EPS
		  WTU_EVG(8,I)  = EPS
		  WTU_EVG(9,I)  = 1.+EPS
		  WTU_EVG(10,I) = -EPS
		  
		  
	  ELSE !NORTH OR SOUTH OR OTHERS
		  
		IF (COORD.EQ.'XY') THEN
		  AT = YU(NJP1(NAGU_EVG(1,I)))-YU(NAGU_EVG(1,I))
		  BT = XU(NAGU_EVG(1,I))-XU(NJP1(NAGU_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (XU(N_CTRD_AG+I)-XU(NAGU_EVG(1,I)))
     * -(-BT/AT*(YU(N_CTRD_AG+I)-YU(NAGU_EVG(1,I)))-CT/AT)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT = LATU(NJP1(NAGU_EVG(1,I)))-LATU(NAGU_EVG(1,I))
		  BT = LONU(NAGU_EVG(1,I))-LONU(NJP1(NAGU_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (LONU(N_CTRD_AG+I)-LONU(NAGU_EVG(1,I)))
     * -(-BT/AT*(LATU(N_CTRD_AG+I)-LATU(NAGU_EVG(1,I)))-CT/AT)
          ENDIF
		  
		IF (DIS.GE.0.) THEN
		  NAGU_EVG(2,I) = NIP1(NAGU_EVG(1,I))
		  NAGU_EVG(3,I) = NJP1(NAGU_EVG(1,I))
		  NAGU_EVG(4,I) = NIP1(NJP1(NAGU_EVG(1,I)))
		  NAGU_EVG(5,I) = NJM1(NAGU_EVG(1,I))
		  NAGU_EVG(6,I) = NIP1(NJM1(NAGU_EVG(1,I)))
		ELSE
		  NAGU_EVG(2,I) = NIM1(NAGU_EVG(1,I))
		  NAGU_EVG(3,I) = NJP1(NAGU_EVG(1,I))
		  NAGU_EVG(4,I) = NIM1(NJP1(NAGU_EVG(1,I)))
		  NAGU_EVG(5,I) = NJM1(NAGU_EVG(1,I))
		  NAGU_EVG(6,I) = NIM1(NJM1(NAGU_EVG(1,I)))
		ENDIF
		
		IF (NAGU_EVG(3,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YU(NAGU_EVG(3,I))-YU(NAGU_EVG(1,I))
		    B1 = XU(NAGU_EVG(1,I))-XU(NAGU_EVG(3,I))
		    C1 = XU(N_CTRD_AG+I)*(-A1)+YU(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			A1 = LATU(NAGU_EVG(3,I))-LATU(NAGU_EVG(1,I))
		    B1 = LONU(NAGU_EVG(1,I))-LONU(NAGU_EVG(3,I))
		    C1 = LONU(N_CTRD_AG+I)*(-A1)+LATU(N_CTRD_AG+I)*(-B1)
		  ENDIF
		  
		ELSEIF (NAGU_EVG(5,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YU(NAGU_EVG(5,I))-YU(NAGU_EVG(1,I))
		    B1 = XU(NAGU_EVG(1,I))-XU(NAGU_EVG(5,I))
		    C1 = XU(N_CTRD_AG+I)*(-A1)+YU(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			A1 = LATU(NAGU_EVG(5,I))-LATU(NAGU_EVG(1,I))
		    B1 = LONU(NAGU_EVG(1,I))-LONU(NAGU_EVG(5,I))
		    C1 = LONU(N_CTRD_AG+I)*(-A1)+LATU(N_CTRD_AG+I)*(-B1)
		  ENDIF
		ELSE
		  PRINT*, 'ERROR IN PRE_IE.FOR!'
		  PAUSE
		  STOP
		ENDIF
		
		DO J = 2,6
		  IF (NAGU_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XU(NAGU_EVG(J,I)),YU(NAGU_EVG(J,I))
     * ,XU(1:N_CTRD_AG),YU(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LONU(NAGU_EVG(J,I)),LATU(NAGU_EVG(J,I))
     * ,LONU(1:N_CTRD_AG),LATU(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGU_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		!WEIGHTS FOR FIRST LINEAR INTERPOLATION
		IF (COORD.EQ.'XY') THEN
		  DIS1 = ABS((A1*XU(NAGU_EVG(2,I))+B1*YU(NAGU_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XU(NAGU_EVG(2,I))-XU(NAGU_EVG(1,I)))**2
     * +(YU(NAGU_EVG(2,I))-YU(NAGU_EVG(1,I)))**2)
		  WTU_EVG(1,I) = DIS1/DIS2
		  WTU_EVG(2,I) = 1-WTU_EVG(1,I)
		  
		  DIS1 = ABS((A1*XU(NAGU_EVG(4,I))+B1*YU(NAGU_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XU(NAGU_EVG(4,I))-XU(NAGU_EVG(3,I)))**2
     * +(YU(NAGU_EVG(4,I))-YU(NAGU_EVG(3,I)))**2)
		  WTU_EVG(3,I) = DIS1/DIS2
		  WTU_EVG(4,I) = 1-WTU_EVG(3,I)
		  
		  DIS1 = ABS((A1*XU(NAGU_EVG(6,I))+B1*YU(NAGU_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XU(NAGU_EVG(6,I))-XU(NAGU_EVG(5,I)))**2
     * +(YU(NAGU_EVG(6,I))-YU(NAGU_EVG(5,I)))**2)
		  WTU_EVG(5,I) = DIS1/DIS2
		  WTU_EVG(6,I) = 1-WTU_EVG(5,I)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  DIS1 = ABS((A1*LONU(NAGU_EVG(2,I))+B1*LATU(NAGU_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONU(NAGU_EVG(2,I))-LONU(NAGU_EVG(1,I)))**2
     * +(LATU(NAGU_EVG(2,I))-LATU(NAGU_EVG(1,I)))**2)
		  WTU_EVG(1,I) = DIS1/DIS2
		  WTU_EVG(2,I) = 1-WTU_EVG(1,I)
		  
		  DIS1 = ABS((A1*LONU(NAGU_EVG(4,I))+B1*LATU(NAGU_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONU(NAGU_EVG(4,I))-LONU(NAGU_EVG(3,I)))**2
     * +(LATU(NAGU_EVG(4,I))-LATU(NAGU_EVG(3,I)))**2)
		  WTU_EVG(3,I) = DIS1/DIS2
		  WTU_EVG(4,I) = 1-WTU_EVG(3,I)
		  
		  DIS1 = ABS((A1*LONU(NAGU_EVG(6,I))+B1*LATU(NAGU_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONU(NAGU_EVG(6,I))-LONU(NAGU_EVG(5,I)))**2
     * +(LATU(NAGU_EVG(6,I))-LATU(NAGU_EVG(5,I)))**2)
		  WTU_EVG(5,I) = DIS1/DIS2
		  WTU_EVG(6,I) = 1-WTU_EVG(5,I)
		  
		ENDIF
		
		!WEIGHTS FOR SECOND LINEAR INTERPOLATION
		  RFN = NINT(H2(NAGU_EVG(1,I))/H2(N_CTRD_AG+I)) !REFINING MULTIPLE
		  IF (NJM4(N_CTRD_AG+I).EQ.N_CTRDP1) THEN !SOUTH
		    K = 0
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NJP1(II)
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGU_EVG(1,II-N_CTRD_AG).NE.NAGU_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
			K = RFN-K
			
		  ELSE !NORTH
			K = 1
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NJM1(II)
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGU_EVG(1,II-N_CTRD_AG).NE.NAGU_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
		  ENDIF
		  
		  !EPSILON
		  EPS = (2*K-1.)/2.*H2(N_CTRD_AG+I)/H2(NAGU_EVG(1,I))-0.5
		  
		  !WEIGHTS
		  WTU_EVG(7,I)  = 1.-EPS
		  WTU_EVG(8,I)  = EPS
		  WTU_EVG(9,I)  = 1.+EPS
		  WTU_EVG(10,I) = -EPS
		
		ENDIF
		
!----------------------------------------------------------------------
! NAGV & WTV
	  !DETERMINE AG NUMBERS FOR INTERPOLATION & WEIGHTS
	  IF ((NIP2(NAGV_EVG(1,I)).GT.N_CTRD_AG) .OR. 
     * (NIM2(NAGV_EVG(1,I)).GT.N_CTRD_AG)) THEN !WEST OR EAS
		
		IF (COORD.EQ.'XY') THEN
		  AT = YV(NIP1(NAGV_EVG(1,I)))-YV(NAGV_EVG(1,I))
		  BT = XV(NAGV_EVG(1,I))-XV(NIP1(NAGV_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (YV(N_CTRD_AG+I)-YV(NAGV_EVG(1,I)))
     * -(-AT/BT*(XV(N_CTRD_AG+I)-XV(NAGV_EVG(1,I)))-CT/BT)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT = LATV(NIP1(NAGV_EVG(1,I)))-LATV(NAGV_EVG(1,I))
		  BT = LONV(NAGV_EVG(1,I))-LONV(NIP1(NAGV_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (LATV(N_CTRD_AG+I)-LATV(NAGV_EVG(1,I)))
     * -(-AT/BT*(LONV(N_CTRD_AG+I)-LONV(NAGV_EVG(1,I)))-CT/BT)
          ENDIF
		  
		IF (DIS.GE.0.) THEN
		  NAGV_EVG(2,I) = NJP1(NAGV_EVG(1,I))
		  NAGV_EVG(3,I) = NIP1(NAGV_EVG(1,I))
		  NAGV_EVG(4,I) = NJP1(NIP1(NAGV_EVG(1,I)))
		  NAGV_EVG(5,I) = NIM1(NAGV_EVG(1,I))
		  NAGV_EVG(6,I) = NJP1(NIM1(NAGV_EVG(1,I)))
		ELSE
		  NAGV_EVG(2,I) = NJM1(NAGV_EVG(1,I))
		  NAGV_EVG(3,I) = NIP1(NAGV_EVG(1,I))
		  NAGV_EVG(4,I) = NJM1(NIP1(NAGV_EVG(1,I)))
		  NAGV_EVG(5,I) = NIM1(NAGV_EVG(1,I))
		  NAGV_EVG(6,I) = NJM1(NIM1(NAGV_EVG(1,I)))
		ENDIF
		
		!PERPENDICULAR LINE
		IF (NAGV_EVG(3,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YV(NAGV_EVG(3,I))-YV(NAGV_EVG(1,I))
		    B1 = XV(NAGV_EVG(1,I))-XV(NAGV_EVG(3,I))
		    C1 = XV(N_CTRD_AG+I)*(-A1)+YV(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			A1 = LATV(NAGV_EVG(3,I))-LATV(NAGV_EVG(1,I))
		    B1 = LONV(NAGV_EVG(1,I))-LONV(NAGV_EVG(3,I))
		    C1 = LONV(N_CTRD_AG+I)*(-A1)+LATV(N_CTRD_AG+I)*(-B1)
		  ENDIF
		  
		ELSEIF (NAGV_EVG(5,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YV(NAGV_EVG(5,I))-YV(NAGV_EVG(1,I))
		    B1 = XV(NAGV_EVG(1,I))-XV(NAGV_EVG(5,I))
		    C1 = XV(N_CTRD_AG+I)*(-A1)+YV(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			A1 = LATV(NAGV_EVG(5,I))-LATV(NAGV_EVG(1,I))
		    B1 = LONV(NAGV_EVG(1,I))-LONV(NAGV_EVG(5,I))
		    C1 = LONV(N_CTRD_AG+I)*(-A1)+LATV(N_CTRD_AG+I)*(-B1)
		  ENDIF
		ELSE
		  PRINT*, 'ERROR IN PRE_IE.FOR!'
		  PAUSE
		  STOP
		ENDIF
		
		DO J = 2,6
		  IF (NAGV_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XV(NAGV_EVG(J,I)),YV(NAGV_EVG(J,I))
     * ,XV(1:N_CTRD_AG),YV(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LONV(NAGV_EVG(J,I)),LATV(NAGV_EVG(J,I))
     * ,LONV(1:N_CTRD_AG),LATV(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGV_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		!WEIGHTS FOR FIRST LINEAR INTERPOLATION
		IF (COORD.EQ.'XY') THEN
		  DIS1 = ABS((A1*XV(NAGV_EVG(2,I))+B1*YV(NAGV_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XV(NAGV_EVG(2,I))-XV(NAGV_EVG(1,I)))**2
     * +(YV(NAGV_EVG(2,I))-YV(NAGV_EVG(1,I)))**2)
		  WTV_EVG(1,I) = DIS1/DIS2
		  WTV_EVG(2,I) = 1-WTV_EVG(1,I)
		  
		  DIS1 = ABS((A1*XV(NAGV_EVG(4,I))+B1*YV(NAGV_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XV(NAGV_EVG(4,I))-XV(NAGV_EVG(3,I)))**2
     * +(YV(NAGV_EVG(4,I))-YV(NAGV_EVG(3,I)))**2)
		  WTV_EVG(3,I) = DIS1/DIS2
		  WTV_EVG(4,I) = 1-WTV_EVG(3,I)
		  
		  DIS1 = ABS((A1*XV(NAGV_EVG(6,I))+B1*YV(NAGV_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XV(NAGV_EVG(6,I))-XV(NAGV_EVG(5,I)))**2
     * +(YV(NAGV_EVG(6,I))-YV(NAGV_EVG(5,I)))**2)
		  WTV_EVG(5,I) = DIS1/DIS2
		  WTV_EVG(6,I) = 1-WTV_EVG(5,I)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  DIS1 = ABS((A1*LONV(NAGV_EVG(2,I))+B1*LATV(NAGV_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONV(NAGV_EVG(2,I))-LONV(NAGV_EVG(1,I)))**2
     * +(LATV(NAGV_EVG(2,I))-LATV(NAGV_EVG(1,I)))**2)
		  WTV_EVG(1,I) = DIS1/DIS2
		  WTV_EVG(2,I) = 1-WTV_EVG(1,I)
		  
		  DIS1 = ABS((A1*LONV(NAGV_EVG(4,I))+B1*LATV(NAGV_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONV(NAGV_EVG(4,I))-LONV(NAGV_EVG(3,I)))**2
     * +(LATV(NAGV_EVG(4,I))-LATV(NAGV_EVG(3,I)))**2)
		  WTV_EVG(3,I) = DIS1/DIS2
		  WTV_EVG(4,I) = 1-WTV_EVG(3,I)
		  
		  DIS1 = ABS((A1*LONV(NAGV_EVG(6,I))+B1*LATV(NAGV_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONV(NAGV_EVG(6,I))-LONV(NAGV_EVG(5,I)))**2
     * +(LATV(NAGV_EVG(6,I))-LATV(NAGV_EVG(5,I)))**2)
		  WTV_EVG(5,I) = DIS1/DIS2
		  WTV_EVG(6,I) = 1-WTV_EVG(5,I)
		  
		ENDIF
		
		!WEIGHTS FOR SECOND LINEAR INTERPOLATION
		  RFN = NINT(H1(NAGV_EVG(1,I))/H1(N_CTRD_AG+I)) !REFINING MULTIPLE
		  IF (NIP2(NAGV_EVG(1,I)).GT.N_CTRD_AG) THEN !WEST
		    K = 0
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NIP1(II)
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGV_EVG(1,II-N_CTRD_AG).NE.NAGV_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
			K = RFN-K
			
		  ELSE !EAST
			K = 1
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NIM1(II)
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGV_EVG(1,II-N_CTRD_AG).NE.NAGV_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
		  ENDIF
		  
		  !EPSILON
		  EPS = (2*K-1.)/2.*H1(N_CTRD_AG+I)/H1(NAGV_EVG(1,I))-0.5
		  
		  !WEIGHTS
		  WTV_EVG(7,I)  = 1.-EPS
		  WTV_EVG(8,I)  = EPS
		  WTV_EVG(9,I)  = 1.+EPS
		  WTV_EVG(10,I) = -EPS

		  
	  ELSE !NORTH OR SOUTH OR OTHERS
		  
		IF (COORD.EQ.'XY') THEN
		  AT = YV(NJP1(NAGV_EVG(1,I)))-YV(NAGV_EVG(1,I))
		  BT = XV(NAGV_EVG(1,I))-XV(NJP1(NAGV_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (XV(N_CTRD_AG+I)-XV(NAGV_EVG(1,I)))
     * -(-BT/AT*(YV(N_CTRD_AG+I)-YV(NAGV_EVG(1,I)))-CT/AT)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT = LATV(NJP1(NAGV_EVG(1,I)))-LATV(NAGV_EVG(1,I))
		  BT = LONV(NAGV_EVG(1,I))-LONV(NJP1(NAGV_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (LONV(N_CTRD_AG+I)-LONV(NAGV_EVG(1,I)))
     * -(-BT/AT*(LATV(N_CTRD_AG+I)-LATV(NAGV_EVG(1,I)))-CT/AT)
          ENDIF
		  
		IF (DIS.GE.0.) THEN
		  NAGV_EVG(2,I) = NIP1(NAGV_EVG(1,I))
		  NAGV_EVG(3,I) = NJP1(NAGV_EVG(1,I))
		  NAGV_EVG(4,I) = NIP1(NJP1(NAGV_EVG(1,I)))
		  NAGV_EVG(5,I) = NJM1(NAGV_EVG(1,I))
		  NAGV_EVG(6,I) = NIP1(NJM1(NAGV_EVG(1,I)))
		ELSE
		  NAGV_EVG(2,I) = NIM1(NAGV_EVG(1,I))
		  NAGV_EVG(3,I) = NJP1(NAGV_EVG(1,I))
		  NAGV_EVG(4,I) = NIM1(NJP1(NAGV_EVG(1,I)))
		  NAGV_EVG(5,I) = NJM1(NAGV_EVG(1,I))
		  NAGV_EVG(6,I) = NIM1(NJM1(NAGV_EVG(1,I)))
		ENDIF
		
		IF (NAGV_EVG(3,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YV(NAGV_EVG(3,I))-YV(NAGV_EVG(1,I))
		    B1 = XV(NAGV_EVG(1,I))-XV(NAGV_EVG(3,I))
		    C1 = XV(N_CTRD_AG+I)*(-A1)+YV(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			A1 = LATV(NAGV_EVG(3,I))-LATV(NAGV_EVG(1,I))
		    B1 = LONV(NAGV_EVG(1,I))-LONV(NAGV_EVG(3,I))
		    C1 = LONV(N_CTRD_AG+I)*(-A1)+LATV(N_CTRD_AG+I)*(-B1)
		  ENDIF
		  
		ELSEIF (NAGV_EVG(5,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YV(NAGV_EVG(5,I))-YV(NAGV_EVG(1,I))
		    B1 = XV(NAGV_EVG(1,I))-XV(NAGV_EVG(5,I))
		    C1 = XV(N_CTRD_AG+I)*(-A1)+YV(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			A1 = LATV(NAGV_EVG(5,I))-LATV(NAGV_EVG(1,I))
		    B1 = LONV(NAGV_EVG(1,I))-LONV(NAGV_EVG(5,I))
		    C1 = LONV(N_CTRD_AG+I)*(-A1)+LATV(N_CTRD_AG+I)*(-B1)
		  ENDIF
		ELSE
		  PRINT*, 'ERROR IN PRE_IE.FOR!'
		  PAUSE
		  STOP
		ENDIF
		
		DO J = 2,6
		  IF (NAGV_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XV(NAGV_EVG(J,I)),YV(NAGV_EVG(J,I))
     * ,XV(1:N_CTRD_AG),YV(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LONV(NAGV_EVG(J,I)),LATV(NAGV_EVG(J,I))
     * ,LONV(1:N_CTRD_AG),LATV(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGV_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		!WEIGHTS FOR FIRST LINEAR INTERPOLATION
		IF (COORD.EQ.'XY') THEN
		  DIS1 = ABS((A1*XV(NAGV_EVG(2,I))+B1*YV(NAGV_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XV(NAGV_EVG(2,I))-XV(NAGV_EVG(1,I)))**2
     * +(YV(NAGV_EVG(2,I))-YV(NAGV_EVG(1,I)))**2)
		  WTV_EVG(1,I) = DIS1/DIS2
		  WTV_EVG(2,I) = 1-WTV_EVG(1,I)
		  
		  DIS1 = ABS((A1*XV(NAGV_EVG(4,I))+B1*YV(NAGV_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XV(NAGV_EVG(4,I))-XV(NAGV_EVG(3,I)))**2
     * +(YV(NAGV_EVG(4,I))-YV(NAGV_EVG(3,I)))**2)
		  WTV_EVG(3,I) = DIS1/DIS2
		  WTV_EVG(4,I) = 1-WTV_EVG(3,I)
		  
		  DIS1 = ABS((A1*XV(NAGV_EVG(6,I))+B1*YV(NAGV_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((XV(NAGV_EVG(6,I))-XV(NAGV_EVG(5,I)))**2
     * +(YV(NAGV_EVG(6,I))-YV(NAGV_EVG(5,I)))**2)
		  WTV_EVG(5,I) = DIS1/DIS2
		  WTV_EVG(6,I) = 1-WTV_EVG(5,I)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  DIS1 = ABS((A1*LONV(NAGV_EVG(2,I))+B1*LATV(NAGV_EVG(2,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONV(NAGV_EVG(2,I))-LONV(NAGV_EVG(1,I)))**2
     * +(LATV(NAGV_EVG(2,I))-LATV(NAGV_EVG(1,I)))**2)
		  WTV_EVG(1,I) = DIS1/DIS2
		  WTV_EVG(2,I) = 1-WTV_EVG(1,I)
		  
		  DIS1 = ABS((A1*LONV(NAGV_EVG(4,I))+B1*LATV(NAGV_EVG(4,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONV(NAGV_EVG(4,I))-LONV(NAGV_EVG(3,I)))**2
     * +(LATV(NAGV_EVG(4,I))-LATV(NAGV_EVG(3,I)))**2)
		  WTV_EVG(3,I) = DIS1/DIS2
		  WTV_EVG(4,I) = 1-WTV_EVG(3,I)
		  
		  DIS1 = ABS((A1*LONV(NAGV_EVG(6,I))+B1*LATV(NAGV_EVG(6,I))+C1)
     * /SQRT(A1**2+B1**2))
		  DIS2 = SQRT((LONV(NAGV_EVG(6,I))-LONV(NAGV_EVG(5,I)))**2
     * +(LATV(NAGV_EVG(6,I))-LATV(NAGV_EVG(5,I)))**2)
		  WTV_EVG(5,I) = DIS1/DIS2
		  WTV_EVG(6,I) = 1-WTV_EVG(5,I)
		  
		ENDIF
		
		!WEIGHTS FOR SECOND LINEAR INTERPOLATION
		  RFN = NINT(H2(NAGV_EVG(1,I))/H2(N_CTRD_AG+I)) !REFINING MULTIPLE
		  IF (NJM4(N_CTRD_AG+I).EQ.N_CTRDP1) THEN !SOUTH
		    K = 0
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NJP1(II)
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGV_EVG(1,II-N_CTRD_AG).NE.NAGV_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
			K = RFN-K
			
		  ELSE !NORTH
			K = 1
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
			  II = NJM1(II)
			  IF (II.LE.N_CTRD_AG) THEN
				EXIT
			  ELSEIF (NAGV_EVG(1,II-N_CTRD_AG).NE.NAGV_EVG(1,I)) THEN
				EXIT
			  ENDIF
			  
			  K = K+1
			ENDDO
			
		  ENDIF
		  
		  !EPSILON
		  EPS = (K-1)*H2(N_CTRD_AG+I)/H2(NAGV_EVG(1,I))
		  
		  !WEIGHTS
		  WTV_EVG(7,I)  = 1.-EPS
		  WTV_EVG(8,I)  = EPS
		  WTV_EVG(9,I)  = 1.+EPS
		  WTV_EVG(10,I) = -EPS
		
		ENDIF
	  
      ENDDO
	  
	!ENDIF




!**********************************************************************
! HAEI (HSIMT ADVECTION-EQUIVALENT INTERPOLATION)
#elif defined INT_HAEI

      !FIND IN WHICH AG IS THIS EVG
      DO I = 1,N_CTRD_EVG
	  DO J = 1,N_CTRD_AG
		DO K = 1,4
		  XB(K) = XNODE(K,J)
		  YB(K) = YNODE(K,J)
		ENDDO
		  
		IF (COORD.EQ.'XY') THEN
		  CALL INSIDE(XR(N_CTRD_AG+I),YR(N_CTRD_AG+I),XB,YB,NB,IND)
		ELSE
		  CALL INSIDE(LON(N_CTRD_AG+I),LAT(N_CTRD_AG+I),XB,YB,NB,IND)
		ENDIF
		  
		IF (IND.EQ.1) THEN
		  NAGC_EVG(1,I) = J
		  NAGU_EVG(1,I) = J
		  NAGV_EVG(1,I) = J
		  EXIT
		ELSEIF (J.EQ.N_CTRD_AG) THEN
		  PRINT*, 'CANNOT FIND MATCH GRID FOR EVG',N_CTRD_AG+I
		  PAUSE
		  STOP
		ENDIF
	  ENDDO
      ENDDO
        
      DO I = 1,N_CTRD_EVG
C----------------------------------------------------------------------
C NAGC & WTC
      !DETERMINE AG NUMBERS FOR INTERPOLATION & WEIGHTS
      IF ((NIP2(NAGC_EVG(1,I)).GT.N_CTRD_AG) .OR. 
     * (NIM2(NAGC_EVG(1,I)).GT.N_CTRD_AG)) THEN !WEST OR EAST
        
        IF (NIM1(NAGC_EVG(1,I)).EQ.N_CTRDP1) THEN
		  NAGC_EVG(:,I) = 1
		  NAGU_EVG(:,I) = 1
		  NAGV_EVG(:,I) = 1
		  CYCLE
        ENDIF
        
        IF (COORD.EQ.'XY') THEN
		  AT = YR(NIP1(NAGC_EVG(1,I)))-YR(NAGC_EVG(1,I))
		  BT = XR(NAGC_EVG(1,I))-XR(NIP1(NAGC_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (YR(N_CTRD_AG+I)-YR(NAGC_EVG(1,I)))
     * -(-AT/BT*(XR(N_CTRD_AG+I)-XR(NAGC_EVG(1,I)))-CT/BT)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT = LAT(NIP1(NAGC_EVG(1,I)))-LAT(NAGC_EVG(1,I))
		  BT = LON(NAGC_EVG(1,I))-LON(NIP1(NAGC_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (LAT(N_CTRD_AG+I)-LAT(NAGC_EVG(1,I)))
     * -(-AT/BT*(LON(N_CTRD_AG+I)-LON(NAGC_EVG(1,I)))-CT/BT)
        ENDIF
        
        IF (DIS.GE.0.) THEN
          NAGC_EVG(2,I) = NJP1(NAGC_EVG(1,I))
          NAGC_EVG(3,I) = NIM2(NAGC_EVG(1,I))
          NAGC_EVG(4,I) = FindCtrd(NAGC_EVG(1,I),-2,1)
          NAGC_EVG(5,I) = NIM1(NAGC_EVG(1,I))
          NAGC_EVG(6,I) = FindCtrd(NAGC_EVG(1,I),-1,1)
          NAGC_EVG(7,I) = NIP1(NAGC_EVG(1,I))
          NAGC_EVG(8,I) = FindCtrd(NAGC_EVG(1,I),1,1)
          NAGC_EVG(9,I) = NIP2(NAGC_EVG(1,I))
          NAGC_EVG(10,I) = FindCtrd(NAGC_EVG(1,I),2,1)
        ELSE
          NAGC_EVG(2,I) = NJM1(NAGC_EVG(1,I))
          NAGC_EVG(3,I) = NIM2(NAGC_EVG(1,I))
          NAGC_EVG(4,I) = FindCtrd(NAGC_EVG(1,I),-2,-1)
          NAGC_EVG(5,I) = NIM1(NAGC_EVG(1,I))
          NAGC_EVG(6,I) = FindCtrd(NAGC_EVG(1,I),-1,-1)
          NAGC_EVG(7,I) = NIP1(NAGC_EVG(1,I))
          NAGC_EVG(8,I) = FindCtrd(NAGC_EVG(1,I),1,-1)
          NAGC_EVG(9,I) = NIP2(NAGC_EVG(1,I))
          NAGC_EVG(10,I) = FindCtrd(NAGC_EVG(1,I),2,-1)
        ENDIF
        
        !PERPENDICULAR LINE
		IF (NAGC_EVG(5,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YR(NAGC_EVG(5,I))-YR(NAGC_EVG(1,I))
		    B1 = XR(NAGC_EVG(1,I))-XR(NAGC_EVG(5,I))
		    C1 = XR(N_CTRD_AG+I)*(-A1)+YR(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LAT(NAGC_EVG(5,I))-LAT(NAGC_EVG(1,I))
		    B1 = LON(NAGC_EVG(1,I))-LON(NAGC_EVG(5,I))
		    C1 = LON(N_CTRD_AG+I)*(-A1)+LAT(N_CTRD_AG+I)*(-B1)
		  ENDIF
		  
		ELSEIF (NAGC_EVG(7,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YR(NAGC_EVG(7,I))-YR(NAGC_EVG(1,I))
		    B1 = XR(NAGC_EVG(1,I))-XR(NAGC_EVG(7,I))
		    C1 = XR(N_CTRD_AG+I)*(-A1)+YR(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LAT(NAGC_EVG(7,I))-LAT(NAGC_EVG(1,I))
		    B1 = LON(NAGC_EVG(1,I))-LON(NAGC_EVG(7,I))
		    C1 = LON(N_CTRD_AG+I)*(-A1)+LAT(N_CTRD_AG+I)*(-B1)
		  ENDIF
		ELSE
		  PRINT*, 'ERROR IN PRE_IE.FOR!'
		  PAUSE
		  STOP
		ENDIF
		
		DO J = 2,10
		  IF (NAGC_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XR(NAGC_EVG(J,I)),YR(NAGC_EVG(J,I))
     * ,XR(1:N_CTRD_AG),YR(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LON(NAGC_EVG(J,I)),LAT(NAGC_EVG(J,I))
     * ,LON(1:N_CTRD_AG),LAT(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGC_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		!WEIGHTS FOR LINEAR INTERPOLATION
		IF (COORD.EQ.'XY') THEN
          DO J = 1,5
            DIS1 = 
     * ABS((A1*XR(NAGC_EVG(2*J,I))+B1*YR(NAGC_EVG(2*J,I))+C1)
     * /SQRT(A1**2+B1**2))
				DIS2 = 
     * ABS((A1*XR(NAGC_EVG(2*J-1,I))+B1*YR(NAGC_EVG(2*J-1,I))+C1)
     * /SQRT(A1**2+B1**2))
				WTC_EVG(2*J-1,I) = DIS1/(DIS1+DIS2+EPSILON)
				WTC_EVG(2*J,I) = 1-WTC_EVG(2*J-1,I)
          ENDDO
		  
        ELSEIF (COORD.EQ.'BL') THEN
          DO J = 1,5
            DIS1 = 
     * ABS((A1*LON(NAGC_EVG(2*J,I))+B1*LAT(NAGC_EVG(2*J,I))+C1)
     * /SQRT(A1**2+B1**2))
				DIS2 = 
     * ABS((A1*LON(NAGC_EVG(2*J-1,I))+B1*LAT(NAGC_EVG(2*J-1,I))+C1)
     * /SQRT(A1**2+B1**2))
				WTC_EVG(2*J-1,I) = DIS1/(DIS1+DIS2+EPSILON)
				WTC_EVG(2*J,I) = 1-WTC_EVG(2*J-1,I)
          ENDDO
          
        ENDIF
        
        !DEFINE EPSC
        RFN = NINT(H1(NAGC_EVG(1,I))/H1(N_CTRD_AG+I)) !REFINING MULTIPLE
		IF (NIP2(NAGC_EVG(1,I)).GT.N_CTRD_AG) THEN !WEST
		  K = 0
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
				II = NIP1(II)
				IF (II.LE.N_CTRD_AG) THEN
					EXIT
				ELSEIF (NAGC_EVG(1,II-N_CTRD_AG).NE.NAGC_EVG(1,I)) THEN
					EXIT
				ENDIF
			  
				K = K+1
			ENDDO
			
			K = RFN-K
			
		ELSE !EAST
			K = 1
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
				II = NIM1(II)
				IF (II.LE.N_CTRD_AG) THEN
					EXIT
				ELSEIF (NAGC_EVG(1,II-N_CTRD_AG).NE.NAGC_EVG(1,I)) THEN
					EXIT
				ENDIF
			  
				K = K+1
			ENDDO
			
		ENDIF
		  
		!EPSILON
		EPSC(I) = (2*K-1.)/2.*H1(N_CTRD_AG+I)/H1(NAGC_EVG(1,I))-0.5
        
        
      ELSE !NORTH OR SOUTH OR OTHERS
     
        IF (COORD.EQ.'XY') THEN
		  AT = YR(NJP1(NAGC_EVG(1,I)))-YR(NAGC_EVG(1,I))
		  BT = XR(NAGC_EVG(1,I))-XR(NJP1(NAGC_EVG(1,I)))
		  CT = 0
		  
		  DIS = (XR(N_CTRD_AG+I)-XR(NAGC_EVG(1,I)))
     * -(-BT/AT*(YR(N_CTRD_AG+I)-YR(NAGC_EVG(1,I)))-CT/AT)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT = LAT(NJP1(NAGC_EVG(1,I)))-LAT(NAGC_EVG(1,I))
		  BT = LON(NAGC_EVG(1,I))-LON(NJP1(NAGC_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (LON(N_CTRD_AG+I)-LON(NAGC_EVG(1,I)))
     * -(-BT/AT*(LAT(N_CTRD_AG+I)-LAT(NAGC_EVG(1,I)))-CT/AT)
		  
        ENDIF
        
        IF (DIS.GE.0.) THEN
		  NAGC_EVG(2,I) = NIP1(NAGC_EVG(1,I))
		  NAGC_EVG(3,I) = NJM2(NAGC_EVG(1,I))
          NAGC_EVG(4,I) = FindCtrd(NAGC_EVG(1,I),1,-2)
		  NAGC_EVG(5,I) = NJM1(NAGC_EVG(1,I))
          NAGC_EVG(6,I) = FindCtrd(NAGC_EVG(1,I),1,-1)
          NAGC_EVG(7,I) = NJP1(NAGC_EVG(1,I))
          NAGC_EVG(8,I) = FindCtrd(NAGC_EVG(1,I),1,1)
          NAGC_EVG(9,I) = NJP2(NAGC_EVG(1,I))
          NAGC_EVG(10,I) = FindCtrd(NAGC_EVG(1,I),1,2)
		ELSE
		  NAGC_EVG(2,I) = NIM1(NAGC_EVG(1,I))
		  NAGC_EVG(3,I) = NJM2(NAGC_EVG(1,I))
          NAGC_EVG(4,I) = FindCtrd(NAGC_EVG(1,I),-1,-2)
		  NAGC_EVG(5,I) = NJM1(NAGC_EVG(1,I))
          NAGC_EVG(6,I) = FindCtrd(NAGC_EVG(1,I),-1,-1)
          NAGC_EVG(7,I) = NJP1(NAGC_EVG(1,I))
          NAGC_EVG(8,I) = FindCtrd(NAGC_EVG(1,I),-1,1)
          NAGC_EVG(9,I) = NJP2(NAGC_EVG(1,I))
          NAGC_EVG(10,I) = FindCtrd(NAGC_EVG(1,I),-1,2)
        ENDIF
        
        IF (NAGC_EVG(5,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YR(NAGC_EVG(5,I))-YR(NAGC_EVG(1,I))
		    B1 = XR(NAGC_EVG(1,I))-XR(NAGC_EVG(5,I))
		    C1 = XR(N_CTRD_AG+I)*(-A1)+YR(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LAT(NAGC_EVG(5,I))-LAT(NAGC_EVG(1,I))
		    B1 = LON(NAGC_EVG(1,I))-LON(NAGC_EVG(5,I))
		    C1 = LON(N_CTRD_AG+I)*(-A1)+LAT(N_CTRD_AG+I)*(-B1)
		  ENDIF
		  
		ELSEIF (NAGC_EVG(7,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YR(NAGC_EVG(7,I))-YR(NAGC_EVG(1,I))
		    B1 = XR(NAGC_EVG(1,I))-XR(NAGC_EVG(7,I))
		    C1 = XR(N_CTRD_AG+I)*(-A1)+YR(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LAT(NAGC_EVG(7,I))-LAT(NAGC_EVG(1,I))
		    B1 = LON(NAGC_EVG(1,I))-LON(NAGC_EVG(7,I))
		    C1 = LON(N_CTRD_AG+I)*(-A1)+LAT(N_CTRD_AG+I)*(-B1)
		  ENDIF
		ELSE
		  PRINT*, 'ERROR IN PRE_IE.FOR!'
		  PAUSE
		  STOP
		ENDIF
		
		DO J = 2,10
		  IF (NAGC_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XR(NAGC_EVG(J,I)),YR(NAGC_EVG(J,I))
     * ,XR(1:N_CTRD_AG),YR(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LON(NAGC_EVG(J,I)),LAT(NAGC_EVG(J,I))
     * ,LON(1:N_CTRD_AG),LAT(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGC_EVG(J,I) = N_AG
		  ENDIF
        ENDDO
        
        !WEIGHTS FOR LINEAR INTERPOLATION
		IF (COORD.EQ.'XY') THEN
          DO J = 1,5
            DIS1 = 
     * ABS((A1*XR(NAGC_EVG(2*J,I))+B1*YR(NAGC_EVG(2*J,I))+C1)
     * /SQRT(A1**2+B1**2))
            DIS2 = 
     * ABS((A1*XR(NAGC_EVG(2*J-1,I))+B1*YR(NAGC_EVG(2*J-1,I))+C1)
     * /SQRT(A1**2+B1**2))
				WTC_EVG(2*J-1,I) = DIS1/(DIS1+DIS2+EPSILON)
				WTC_EVG(2*J,I) = 1-WTC_EVG(2*J-1,I)
          ENDDO
		  
        ELSEIF (COORD.EQ.'BL') THEN
          DO J = 1,5
            DIS1 = 
     * ABS((A1*LON(NAGC_EVG(2*J,I))+B1*LAT(NAGC_EVG(2*J,I))+C1)
     * /SQRT(A1**2+B1**2))
				DIS2 = 
     * ABS((A1*LON(NAGC_EVG(2*J-1,I))+B1*LAT(NAGC_EVG(2*J-1,I))+C1)
     * /SQRT(A1**2+B1**2))
				WTC_EVG(2*J-1,I) = DIS1/(DIS1+DIS2+EPSILON)
				WTC_EVG(2*J,I) = 1-WTC_EVG(2*J-1,I)
          ENDDO
          
        ENDIF
        
        !DEFINE EPSILON
        RFN = NINT(H2(NAGC_EVG(1,I))/H2(N_CTRD_AG+I)) !REFINING MULTIPLE
		IF (NJM4(N_CTRD_AG+I).EQ.N_CTRDP1) THEN !SOUTH
		  K = 0
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
				II = NJP1(II)
				IF (II.LE.N_CTRD_AG) THEN
					EXIT
				ELSEIF (NAGC_EVG(1,II-N_CTRD_AG).NE.NAGC_EVG(1,I)) THEN
					EXIT
				ENDIF
			  
				K = K+1
			ENDDO
			
			K = RFN-K
			
		ELSE !NORTH
			K = 1
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
				II = NJM1(II)
				IF (II.LE.N_CTRD_AG) THEN
				EXIT
				ELSEIF (NAGC_EVG(1,II-N_CTRD_AG).NE.NAGC_EVG(1,I)) THEN
				EXIT
				ENDIF
			  
				K = K+1
			ENDDO
			
		ENDIF
		  
		!EPSILON
		EPSC(I) = (2*K-1)/2.*H2(N_CTRD_AG+I)/H2(NAGC_EVG(1,I))-0.5
        
      ENDIF  
      
C----------------------------------------------------------------------
C NAGU & WTU
      !DETERMINE AG NUMBERS FOR INTERPOLATION & WEIGHTS
	  IF ((NIP2(NAGU_EVG(1,I)).GT.N_CTRD_AG) .OR. 
     * (NIM2(NAGU_EVG(1,I)).GT.N_CTRD_AG)) THEN !WEST OR EAST
		
		IF (COORD.EQ.'XY') THEN
		  AT = YU(NIP1(NAGU_EVG(1,I)))-YU(NAGU_EVG(1,I))
		  BT = XU(NAGU_EVG(1,I))-XU(NIP1(NAGU_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (YU(N_CTRD_AG+I)-YU(NAGU_EVG(1,I)))
     * -(-AT/BT*(XU(N_CTRD_AG+I)-XU(NAGU_EVG(1,I)))-CT/BT)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT = LATU(NIP1(NAGU_EVG(1,I)))-LATU(NAGU_EVG(1,I))
		  BT = LONU(NAGU_EVG(1,I))-LONU(NIP1(NAGU_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (LATU(N_CTRD_AG+I)-LATU(NAGU_EVG(1,I)))
     * -(-AT/BT*(LONU(N_CTRD_AG+I)-LONU(NAGU_EVG(1,I)))-CT/BT)
          ENDIF
        
        IF (DIS.GE.0.) THEN
          NAGU_EVG(2,I) = NJP1(NAGU_EVG(1,I))
          NAGU_EVG(3,I) = NIM2(NAGU_EVG(1,I))
          NAGU_EVG(4,I) = FindCtrd(NAGU_EVG(1,I),-2,1)
          NAGU_EVG(5,I) = NIM1(NAGU_EVG(1,I))
          NAGU_EVG(6,I) = FindCtrd(NAGU_EVG(1,I),-1,1)
          NAGU_EVG(7,I) = NIP1(NAGU_EVG(1,I))
          NAGU_EVG(8,I) = FindCtrd(NAGU_EVG(1,I),1,1)
          NAGU_EVG(9,I) = NIP2(NAGU_EVG(1,I))
          NAGU_EVG(10,I) = FindCtrd(NAGU_EVG(1,I),2,1)
        ELSE
          NAGU_EVG(2,I) = NJM1(NAGU_EVG(1,I))
          NAGU_EVG(3,I) = NIM2(NAGU_EVG(1,I))
          NAGU_EVG(4,I) = FindCtrd(NAGU_EVG(1,I),-2,-1)
          NAGU_EVG(5,I) = NIM1(NAGU_EVG(1,I))
          NAGU_EVG(6,I) = FindCtrd(NAGU_EVG(1,I),-1,-1)
          NAGU_EVG(7,I) = NIP1(NAGU_EVG(1,I))
          NAGU_EVG(8,I) = FindCtrd(NAGU_EVG(1,I),1,-1)
          NAGU_EVG(9,I) = NIP2(NAGU_EVG(1,I))
          NAGU_EVG(10,I) = FindCtrd(NAGU_EVG(1,I),2,-1)
        ENDIF
		
		!PERPENDICULAR LINE
		IF (NAGU_EVG(5,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YU(NAGU_EVG(5,I))-YU(NAGU_EVG(1,I))
		    B1 = XU(NAGU_EVG(1,I))-XU(NAGU_EVG(5,I))
		    C1 = XU(N_CTRD_AG+I)*(-A1)+YU(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATU(NAGU_EVG(5,I))-LATU(NAGU_EVG(1,I))
		    B1 = LONU(NAGU_EVG(1,I))-LONU(NAGU_EVG(5,I))
		    C1 = LONU(N_CTRD_AG+I)*(-A1)+LATU(N_CTRD_AG+I)*(-B1)
		  ENDIF
		  
		ELSEIF (NAGU_EVG(7,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YU(NAGU_EVG(7,I))-YU(NAGU_EVG(1,I))
		    B1 = XU(NAGU_EVG(1,I))-XU(NAGU_EVG(7,I))
		    C1 = XU(N_CTRD_AG+I)*(-A1)+YU(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATU(NAGU_EVG(7,I))-LATU(NAGU_EVG(1,I))
		    B1 = LONU(NAGU_EVG(1,I))-LONU(NAGU_EVG(7,I))
		    C1 = LONU(N_CTRD_AG+I)*(-A1)+LATU(N_CTRD_AG+I)*(-B1)
		  ENDIF
		ELSE
		  PRINT*, 'ERROR IN PRE_IE.FOR!'
		  PAUSE
		  STOP
		ENDIF
		
		DO J = 2,10
		  IF (NAGU_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XU(NAGU_EVG(J,I)),YU(NAGU_EVG(J,I))
     * ,XU(1:N_CTRD_AG),YU(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LONU(NAGU_EVG(J,I)),LATU(NAGU_EVG(J,I))
     * ,LONU(1:N_CTRD_AG),LATU(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGU_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		!WEIGHTS FOR LINEAR INTERPOLATION
		IF (COORD.EQ.'XY') THEN
          DO J = 1,5
            DIS1 = 
     * ABS((A1*XU(NAGU_EVG(2*J,I))+B1*YU(NAGU_EVG(2*J,I))+C1)
     * /SQRT(A1**2+B1**2))
				DIS2 = SQRT((XU(NAGU_EVG(2*J,I))-XU(NAGU_EVG(2*J-1,I)))**2
     * +(YU(NAGU_EVG(2*J,I))-YU(NAGU_EVG(2*J-1,I)))**2)
				WTU_EVG(2*J-1,I) = DIS1/(DIS2+EPSILON)
				WTU_EVG(2*J,I) = 1-WTU_EVG(2*J-1,I)
          ENDDO
          
        ELSEIF (COORD.EQ.'BL') THEN
          DO J = 1,5
            DIS1 = 
     * ABS((A1*LONU(NAGU_EVG(2*J,I))+B1*LATU(NAGU_EVG(2*J,I))+C1)
     * /SQRT(A1**2+B1**2))
				DIS2 = SQRT((LONU(NAGU_EVG(2*J,I))-LONU(NAGU_EVG(2*J-1,I)))**2
     * +(LATU(NAGU_EVG(2*J,I))-LATU(NAGU_EVG(2*J-1,I)))**2)
				WTU_EVG(2*J-1,I) = DIS1/(DIS2+EPSILON)
				WTU_EVG(2*J,I) = 1-WTU_EVG(2*J-1,I)
          ENDDO
        ENDIF
        
        !DEFINE EPSU
        RFN = NINT(H1(NAGU_EVG(1,I))/H1(N_CTRD_AG+I)) !REFINING MULTIPLE
		IF (NIP2(NAGU_EVG(1,I)).GT.N_CTRD_AG) THEN !WEST
		  K = 0
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
				II = NIP1(II)
				IF (II.LE.N_CTRD_AG) THEN
					EXIT
				ELSEIF (NAGU_EVG(1,II-N_CTRD_AG).NE.NAGU_EVG(1,I)) THEN
					EXIT
				ENDIF
			  
				K = K+1
			ENDDO
			
			K = RFN-K
			
		ELSE !EAST
			K = 1
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
				II = NIM1(II)
				IF (II.LE.N_CTRD_AG) THEN
					EXIT
				ELSEIF (NAGU_EVG(1,II-N_CTRD_AG).NE.NAGU_EVG(1,I)) THEN
					EXIT
				ENDIF
			  
				K = K+1
			ENDDO
			
		ENDIF
		  
		!EPSILON
		EPSU(I) = (K-1)*H1(N_CTRD_AG+I)/H1(NAGU_EVG(1,I))
      
      ELSE !NORTH OR SOUTH OR OTHERS
		  
		IF (COORD.EQ.'XY') THEN
		  AT = YU(NJP1(NAGU_EVG(1,I)))-YU(NAGU_EVG(1,I))
		  BT = XU(NAGU_EVG(1,I))-XU(NJP1(NAGU_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (XU(N_CTRD_AG+I)-XU(NAGU_EVG(1,I)))
     * -(-BT/AT*(YU(N_CTRD_AG+I)-YU(NAGU_EVG(1,I)))-CT/AT)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT = LATU(NJP1(NAGU_EVG(1,I)))-LATU(NAGU_EVG(1,I))
		  BT = LONU(NAGU_EVG(1,I))-LONU(NJP1(NAGU_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (LONU(N_CTRD_AG+I)-LONU(NAGU_EVG(1,I)))
     * -(-BT/AT*(LATU(N_CTRD_AG+I)-LATU(NAGU_EVG(1,I)))-CT/AT)
        ENDIF
        
        IF (DIS.GE.0.) THEN
		  NAGU_EVG(2,I) = NIP1(NAGU_EVG(1,I))
		  NAGU_EVG(3,I) = NJM2(NAGU_EVG(1,I))
          NAGU_EVG(4,I) = FindCtrd(NAGU_EVG(1,I),1,-2)
		  NAGU_EVG(5,I) = NJM1(NAGU_EVG(1,I))
          NAGU_EVG(6,I) = FindCtrd(NAGU_EVG(1,I),1,-1)
          NAGU_EVG(7,I) = NJP1(NAGU_EVG(1,I))
          NAGU_EVG(8,I) = FindCtrd(NAGU_EVG(1,I),1,1)
          NAGU_EVG(9,I) = NJP2(NAGU_EVG(1,I))
          NAGU_EVG(10,I) = FindCtrd(NAGU_EVG(1,I),1,2)
		ELSE
		  NAGU_EVG(2,I) = NIM1(NAGU_EVG(1,I))
		  NAGU_EVG(3,I) = NJM2(NAGU_EVG(1,I))
          NAGU_EVG(4,I) = FindCtrd(NAGU_EVG(1,I),-1,-2)
		  NAGU_EVG(5,I) = NJM1(NAGU_EVG(1,I))
          NAGU_EVG(6,I) = FindCtrd(NAGU_EVG(1,I),-1,-1)
          NAGU_EVG(7,I) = NJP1(NAGU_EVG(1,I))
          NAGU_EVG(8,I) = FindCtrd(NAGU_EVG(1,I),-1,1)
          NAGU_EVG(9,I) = NJP2(NAGU_EVG(1,I))
          NAGU_EVG(10,I) = FindCtrd(NAGU_EVG(1,I),-1,2)
        ENDIF
		
		IF (NAGU_EVG(5,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YU(NAGU_EVG(5,I))-YU(NAGU_EVG(1,I))
		    B1 = XU(NAGU_EVG(1,I))-XU(NAGU_EVG(5,I))
		    C1 = XU(N_CTRD_AG+I)*(-A1)+YU(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATU(NAGU_EVG(5,I))-LATU(NAGU_EVG(1,I))
		    B1 = LONU(NAGU_EVG(1,I))-LONU(NAGU_EVG(5,I))
		    C1 = LONU(N_CTRD_AG+I)*(-A1)+LATU(N_CTRD_AG+I)*(-B1)
		  ENDIF
		  
		ELSEIF (NAGU_EVG(7,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YU(NAGU_EVG(7,I))-YU(NAGU_EVG(1,I))
		    B1 = XU(NAGU_EVG(1,I))-XU(NAGU_EVG(7,I))
		    C1 = XU(N_CTRD_AG+I)*(-A1)+YU(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATU(NAGU_EVG(7,I))-LATU(NAGU_EVG(1,I))
		    B1 = LONU(NAGU_EVG(1,I))-LONU(NAGU_EVG(7,I))
		    C1 = LONU(N_CTRD_AG+I)*(-A1)+LATU(N_CTRD_AG+I)*(-B1)
		  ENDIF
		ELSE
		  PRINT*, 'ERROR IN PRE_IE.FOR!'
		  PAUSE
		  STOP
		ENDIF
		
		DO J = 2,10
		  IF (NAGU_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XU(NAGU_EVG(J,I)),YU(NAGU_EVG(J,I))
     * ,XU(1:N_CTRD_AG),YU(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LONU(NAGU_EVG(J,I)),LATU(NAGU_EVG(J,I))
     * ,LONU(1:N_CTRD_AG),LATU(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGU_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		!WEIGHTS FOR LINEAR INTERPOLATION
		IF (COORD.EQ.'XY') THEN
          DO J = 1,5
            DIS1 = 
     * ABS((A1*XU(NAGU_EVG(2*J,I))+B1*YU(NAGU_EVG(2*J,I))+C1)
     * /SQRT(A1**2+B1**2))
				DIS2 = SQRT((XU(NAGU_EVG(2*J,I))-XU(NAGU_EVG(2*J-1,I)))**2
     * +(YU(NAGU_EVG(2*J,I))-YU(NAGU_EVG(2*J-1,I)))**2)
				WTU_EVG(2*J-1,I) = DIS1/(DIS2+EPSILON)
				WTU_EVG(2*J,I) = 1-WTU_EVG(2*J-1,I)
          ENDDO
		  
        ELSEIF (COORD.EQ.'BL') THEN
          DO J = 1,5
            DIS1 = 
     * ABS((A1*LONU(NAGU_EVG(2*J,I))+B1*LATU(NAGU_EVG(2*J,I))+C1)
     * /SQRT(A1**2+B1**2))
				DIS2 = 
     * SQRT((LONU(NAGU_EVG(2*J,I))-LONU(NAGU_EVG(2*J-1,I)))**2
     * +(LATU(NAGU_EVG(2*J,I))-LATU(NAGU_EVG(2*J-1,I)))**2)
				WTU_EVG(2*J-1,I) = DIS1/(DIS2+EPSILON)
				WTU_EVG(2*J,I) = 1-WTU_EVG(2*J-1,I)
          ENDDO
        ENDIF
        
        !DEFINE EPSU
        RFN = NINT(H2(NAGU_EVG(1,I))/H2(N_CTRD_AG+I)) !REFINING MULTIPLE
		IF (NJM4(N_CTRD_AG+I).EQ.N_CTRDP1) THEN !SOUTH
		  K = 0
		II = N_CTRD_AG+I
		DO WHILE (.TRUE.)
			II = NJP1(II)
			IF (II.LE.N_CTRD_AG) THEN
			EXIT
			ELSEIF (NAGU_EVG(1,II-N_CTRD_AG).NE.NAGU_EVG(1,I)) THEN
			EXIT
			ENDIF
			  
			K = K+1
		ENDDO
			
		K = RFN-K
			
		ELSE !NORTH
			K = 1
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
				II = NJM1(II)
				IF (II.LE.N_CTRD_AG) THEN
					EXIT
				ELSEIF (NAGU_EVG(1,II-N_CTRD_AG).NE.NAGU_EVG(1,I)) THEN
					EXIT
				ENDIF
			  
				K = K+1
			ENDDO
			
		ENDIF
		  
		!EPSILON
		EPSU(I) = (2*K-1)/2.*H2(N_CTRD_AG+I)/H2(NAGU_EVG(1,I))-0.5
      
      ENDIF
      
      
C----------------------------------------------------------------------
C NAGV & WTV
      
      !DETERMINE AG NUMBERS FOR INTERPOLATION & WEIGHTS
	  IF ((NIP2(NAGV_EVG(1,I)).GT.N_CTRD_AG) .OR. 
     * (NIM2(NAGV_EVG(1,I)).GT.N_CTRD_AG)) THEN !WEST OR EAST
		
		IF (COORD.EQ.'XY') THEN
		  AT = YV(NIP1(NAGV_EVG(1,I)))-YV(NAGV_EVG(1,I))
		  BT = XV(NAGV_EVG(1,I))-XV(NIP1(NAGV_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (YV(N_CTRD_AG+I)-YV(NAGV_EVG(1,I)))
     * -(-AT/BT*(XV(N_CTRD_AG+I)-XV(NAGV_EVG(1,I)))-CT/BT)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT = LATV(NIP1(NAGV_EVG(1,I)))-LATV(NAGV_EVG(1,I))
		  BT = LONV(NAGV_EVG(1,I))-LONV(NIP1(NAGV_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (LATV(N_CTRD_AG+I)-LATV(NAGV_EVG(1,I)))
     * -(-AT/BT*(LONV(N_CTRD_AG+I)-LONV(NAGV_EVG(1,I)))-CT/BT)
        ENDIF
        
        IF (DIS.GE.0.) THEN
          NAGV_EVG(2,I) = NJP1(NAGV_EVG(1,I))
          NAGV_EVG(3,I) = NIM2(NAGV_EVG(1,I))
          NAGV_EVG(4,I) = FindCtrd(NAGV_EVG(1,I),-2,1)
          NAGV_EVG(5,I) = NIM1(NAGV_EVG(1,I))
          NAGV_EVG(6,I) = FindCtrd(NAGV_EVG(1,I),-1,1)
          NAGV_EVG(7,I) = NIP1(NAGV_EVG(1,I))
          NAGV_EVG(8,I) = FindCtrd(NAGV_EVG(1,I),1,1)
          NAGV_EVG(9,I) = NIP2(NAGV_EVG(1,I))
          NAGV_EVG(10,I) = FindCtrd(NAGV_EVG(1,I),2,1)
        ELSE
          NAGV_EVG(2,I) = NJM1(NAGV_EVG(1,I))
          NAGV_EVG(3,I) = NIM2(NAGV_EVG(1,I))
          NAGV_EVG(4,I) = FindCtrd(NAGV_EVG(1,I),-2,-1)
          NAGV_EVG(5,I) = NIM1(NAGV_EVG(1,I))
          NAGV_EVG(6,I) = FindCtrd(NAGV_EVG(1,I),-1,-1)
          NAGV_EVG(7,I) = NIP1(NAGV_EVG(1,I))
          NAGV_EVG(8,I) = FindCtrd(NAGV_EVG(1,I),1,-1)
          NAGV_EVG(9,I) = NIP2(NAGV_EVG(1,I))
          NAGV_EVG(10,I) = FindCtrd(NAGV_EVG(1,I),2,-1)
        ENDIF
		
		!PERPENDICULAR LINE
		IF (NAGV_EVG(5,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YV(NAGV_EVG(5,I))-YV(NAGV_EVG(1,I))
		    B1 = XV(NAGV_EVG(1,I))-XV(NAGV_EVG(5,I))
		    C1 = XV(N_CTRD_AG+I)*(-A1)+YV(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATV(NAGV_EVG(5,I))-LATV(NAGV_EVG(1,I))
		    B1 = LONV(NAGV_EVG(1,I))-LONV(NAGV_EVG(5,I))
		    C1 = LONV(N_CTRD_AG+I)*(-A1)+LATV(N_CTRD_AG+I)*(-B1)
		  ENDIF
		  
		ELSEIF (NAGV_EVG(7,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YV(NAGV_EVG(7,I))-YV(NAGV_EVG(1,I))
		    B1 = XV(NAGV_EVG(1,I))-XV(NAGV_EVG(7,I))
		    C1 = XV(N_CTRD_AG+I)*(-A1)+YV(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATV(NAGV_EVG(7,I))-LATV(NAGV_EVG(1,I))
		    B1 = LONV(NAGV_EVG(1,I))-LONV(NAGV_EVG(7,I))
		    C1 = LONV(N_CTRD_AG+I)*(-A1)+LATV(N_CTRD_AG+I)*(-B1)
		  ENDIF
		ELSE
		  PRINT*, 'ERROR IN PRE_IE.FOR!'
		  PAUSE
		  STOP
		ENDIF
		
		DO J = 2,10
		  IF (NAGV_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XV(NAGV_EVG(J,I)),YV(NAGV_EVG(J,I))
     * ,XV(1:N_CTRD_AG),YV(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LONV(NAGV_EVG(J,I)),LATV(NAGV_EVG(J,I))
     * ,LONV(1:N_CTRD_AG),LATV(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGV_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		!WEIGHTS FOR LINEAR INTERPOLATION
		IF (COORD.EQ.'XY') THEN
          DO J = 1,5
            DIS1 = 
     * ABS((A1*XV(NAGV_EVG(2*J,I))+B1*YV(NAGV_EVG(2*J,I))+C1)
     * /SQRT(A1**2+B1**2))
				DIS2 = SQRT((XV(NAGV_EVG(2*J,I))-XV(NAGV_EVG(2*J-1,I)))**2
     * +(YV(NAGV_EVG(2*J,I))-YV(NAGV_EVG(2*J-1,I)))**2)
				WTV_EVG(2*J-1,I) = DIS1/(DIS2+EPSILON)
				WTV_EVG(2*J,I) = 1-WTV_EVG(2*J-1,I)
          ENDDO
          
        ELSEIF (COORD.EQ.'BL') THEN
          DO J = 1,5
            DIS1 = 
     * ABS((A1*LONV(NAGV_EVG(2*J,I))+B1*LATV(NAGV_EVG(2*J,I))+C1)
     * /SQRT(A1**2+B1**2))
				DIS2 = 
     * SQRT((LONV(NAGV_EVG(2*J,I))-LONV(NAGV_EVG(2*J-1,I)))**2
     * +(LATV(NAGV_EVG(2*J,I))-LATV(NAGV_EVG(2*J-1,I)))**2)
				WTV_EVG(2*J-1,I) = DIS1/(DIS2+EPSILON)
				WTV_EVG(2*J,I) = 1-WTV_EVG(2*J-1,I)
          ENDDO
        ENDIF
        
        !DEFINE EPSV
        RFN = NINT(H1(NAGV_EVG(1,I))/H1(N_CTRD_AG+I)) !REFINING MULTIPLE
		IF (NIP2(NAGV_EVG(1,I)).GT.N_CTRD_AG) THEN !WEST
		  K = 0
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
				II = NIP1(II)
				IF (II.LE.N_CTRD_AG) THEN
					EXIT
				ELSEIF (NAGV_EVG(1,II-N_CTRD_AG).NE.NAGV_EVG(1,I)) THEN
					EXIT
				ENDIF
			  
				K = K+1
			ENDDO
			
			K = RFN-K
			
		ELSE !EAST
			K = 1
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
				II = NIM1(II)
				IF (II.LE.N_CTRD_AG) THEN
					EXIT
				ELSEIF (NAGV_EVG(1,II-N_CTRD_AG).NE.NAGV_EVG(1,I)) THEN
					EXIT
				ENDIF
			  
				K = K+1
			ENDDO
			
		ENDIF
		  
		!EPSILON
		EPSV(I) = (2*K-1)/2.*H1(N_CTRD_AG+I)/H1(NAGV_EVG(1,I))-0.5
      
      ELSE !NORTH OR SOUTH OR OTHERS
		  
		IF (COORD.EQ.'XY') THEN
		  AT = YV(NJP1(NAGV_EVG(1,I)))-YV(NAGV_EVG(1,I))
		  BT = XV(NAGV_EVG(1,I))-XV(NJP1(NAGV_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (XV(N_CTRD_AG+I)-XV(NAGV_EVG(1,I)))
     * -(-BT/AT*(YV(N_CTRD_AG+I)-YV(NAGV_EVG(1,I)))-CT/AT)
		  
		ELSEIF (COORD.EQ.'BL') THEN
		  AT = LATV(NJP1(NAGV_EVG(1,I)))-LATV(NAGV_EVG(1,I))
		  BT = LONV(NAGV_EVG(1,I))-LONV(NJP1(NAGV_EVG(1,I)))
		  CT = 0.
		  
		  DIS = (LONV(N_CTRD_AG+I)-LONV(NAGV_EVG(1,I)))
     * -(-BT/AT*(LATV(N_CTRD_AG+I)-LATV(NAGV_EVG(1,I)))-CT/AT)
        ENDIF
        
        IF (DIS.GE.0.) THEN
		  NAGV_EVG(2,I) = NIP1(NAGV_EVG(1,I))
		  NAGV_EVG(3,I) = NJM2(NAGV_EVG(1,I))
          NAGV_EVG(4,I) = FindCtrd(NAGV_EVG(1,I),1,-2)
		  NAGV_EVG(5,I) = NJM1(NAGV_EVG(1,I))
          NAGV_EVG(6,I) = FindCtrd(NAGV_EVG(1,I),1,-1)
          NAGV_EVG(7,I) = NJP1(NAGV_EVG(1,I))
          NAGV_EVG(8,I) = FindCtrd(NAGV_EVG(1,I),1,1)
          NAGV_EVG(9,I) = NJP2(NAGV_EVG(1,I))
          NAGV_EVG(10,I) = FindCtrd(NAGV_EVG(1,I),1,2)
		ELSE
		  NAGV_EVG(2,I) = NIM1(NAGV_EVG(1,I))
		  NAGV_EVG(3,I) = NJM2(NAGV_EVG(1,I))
          NAGV_EVG(4,I) = FindCtrd(NAGV_EVG(1,I),-1,-2)
		  NAGV_EVG(5,I) = NJM1(NAGV_EVG(1,I))
          NAGV_EVG(6,I) = FindCtrd(NAGV_EVG(1,I),-1,-1)
          NAGV_EVG(7,I) = NJP1(NAGV_EVG(1,I))
          NAGV_EVG(8,I) = FindCtrd(NAGV_EVG(1,I),-1,1)
          NAGV_EVG(9,I) = NJP2(NAGV_EVG(1,I))
          NAGV_EVG(10,I) = FindCtrd(NAGV_EVG(1,I),-1,2)
        ENDIF
		
		IF (NAGV_EVG(5,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YV(NAGV_EVG(5,I))-YV(NAGV_EVG(1,I))
		    B1 = XV(NAGV_EVG(1,I))-XV(NAGV_EVG(5,I))
		    C1 = XV(N_CTRD_AG+I)*(-A1)+YV(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATV(NAGV_EVG(5,I))-LATV(NAGV_EVG(1,I))
		    B1 = LONV(NAGV_EVG(1,I))-LONV(NAGV_EVG(5,I))
		    C1 = LONV(N_CTRD_AG+I)*(-A1)+LATV(N_CTRD_AG+I)*(-B1)
		  ENDIF
		  
		ELSEIF (NAGV_EVG(7,I).LE.N_CTRD_AG) THEN
		  IF (COORD.EQ.'XY') THEN
		    A1 = YV(NAGV_EVG(7,I))-YV(NAGV_EVG(1,I))
		    B1 = XV(NAGV_EVG(1,I))-XV(NAGV_EVG(7,I))
		    C1 = XV(N_CTRD_AG+I)*(-A1)+YV(N_CTRD_AG+I)*(-B1)
		  ELSEIF (COORD.EQ.'BL') THEN
			  A1 = LATV(NAGV_EVG(7,I))-LATV(NAGV_EVG(1,I))
		    B1 = LONV(NAGV_EVG(1,I))-LONV(NAGV_EVG(7,I))
		    C1 = LONV(N_CTRD_AG+I)*(-A1)+LATV(N_CTRD_AG+I)*(-B1)
		  ENDIF
		ELSE
		  PRINT*, 'ERROR IN PRE_IE.FOR!'
		  PAUSE
		  STOP
		ENDIF
		
		DO J = 2,10
		  IF (NAGV_EVG(J,I).GT.N_CTRD_AG) THEN
		    IF (COORD.EQ.'XY') THEN
			  CALL FIND_NEAREST_AG(XV(NAGV_EVG(J,I)),YV(NAGV_EVG(J,I))
     * ,XV(1:N_CTRD_AG),YV(1:N_CTRD_AG),N_AG)
		    ELSEIF (COORD.EQ.'BL') THEN
			  CALL FIND_NEAREST_AG(LONV(NAGV_EVG(J,I)),LATV(NAGV_EVG(J,I))
     * ,LONV(1:N_CTRD_AG),LATV(1:N_CTRD_AG),N_AG)
		    ENDIF
			
		    NAGV_EVG(J,I) = N_AG
		  ENDIF
		ENDDO
		
		!WEIGHTS FOR LINEAR INTERPOLATION
		IF (COORD.EQ.'XY') THEN
          DO J = 1,5
            DIS1 = 
     * ABS((A1*XV(NAGV_EVG(2*J,I))+B1*YV(NAGV_EVG(2*J,I))+C1)
     * /SQRT(A1**2+B1**2))
				DIS2 = SQRT((XV(NAGV_EVG(2*J,I))-XV(NAGV_EVG(2*J-1,I)))**2
     * +(YV(NAGV_EVG(2*J,I))-YV(NAGV_EVG(2*J-1,I)))**2)
				WTV_EVG(2*J-1,I) = DIS1/(DIS2+EPSILON)
				WTV_EVG(2*J,I) = 1-WTV_EVG(2*J-1,I)
          ENDDO
		  
        ELSEIF (COORD.EQ.'BL') THEN
          DO J = 1,5
            DIS1 = 
     * ABS((A1*LONV(NAGV_EVG(2*J,I))+B1*LATV(NAGV_EVG(2*J,I))+C1)
     * /SQRT(A1**2+B1**2))
				DIS2 = 
     * SQRT((LONV(NAGV_EVG(2*J,I))-LONV(NAGV_EVG(2*J-1,I)))**2
     * +(LATV(NAGV_EVG(2*J,I))-LATV(NAGV_EVG(2*J-1,I)))**2)
				WTV_EVG(2*J-1,I) = DIS1/(DIS2+EPSILON)
				WTV_EVG(2*J,I) = 1-WTV_EVG(2*J-1,I)
          ENDDO
        ENDIF
        
        !DEFINE EPSILON
        RFN = NINT(H2(NAGV_EVG(1,I))/H2(N_CTRD_AG+I)) !REFINING MULTIPLE
		IF (NJM4(N_CTRD_AG+I).EQ.N_CTRDP1) THEN !SOUTH
		  K = 0
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
				II = NJP1(II)
				IF (II.LE.N_CTRD_AG) THEN
					EXIT
				ELSEIF (NAGV_EVG(1,II-N_CTRD_AG).NE.NAGV_EVG(1,I)) THEN
					EXIT
				ENDIF
			  
				K = K+1
			ENDDO
			
			K = RFN-K
			
		ELSE !NORTH
			K = 1
			II = N_CTRD_AG+I
			DO WHILE (.TRUE.)
				II = NJM1(II)
				IF (II.LE.N_CTRD_AG) THEN
					EXIT
				ELSEIF (NAGV_EVG(1,II-N_CTRD_AG).NE.NAGV_EVG(1,I)) THEN
					EXIT
				ENDIF
			  
				K = K+1
			ENDDO
			
		ENDIF
		  
		!EPSILON
		EPSV(I) = (K-1)*H2(N_CTRD_AG+I)/H2(NAGV_EVG(1,I))
        
      ENDIF
      
      ENDDO




!**********************************************************************
#else
      PRINT*, 'Interpolation scheme undefined!'
	PAUSE
	STOP
#endif
	
!======================================================================
!                             END: INTERPOLATION
!======================================================================



	  
!======================================================================
!                                 UPDATE
!======================================================================
	
!**********************************************************************
! DIRECTLY-REPLACING
#ifdef UD_DR	  
	  !FIND THE NEAREST GRID
	  DO I = 1,N_CTRD_IVG
		DO J = 1,N_CTRD_AG
		  
		  IF (COORD .EQ. 'XY') THEN
			DISC(J) = SQRT((XR(N_CTRD_AG+N_CTRD_EVG+I)-XR(J))**2
     * +(YR(N_CTRD_AG+N_CTRD_EVG+I)-YR(J))**2)
			DISU(J) = SQRT((XU(N_CTRD_AG+N_CTRD_EVG+I)-XU(J))**2
     * +(YU(N_CTRD_AG+N_CTRD_EVG+I)-YU(J))**2)
			DISV(J) = SQRT((XV(N_CTRD_AG+N_CTRD_EVG+I)-XV(J))**2
     * +(YV(N_CTRD_AG+N_CTRD_EVG+I)-YV(J))**2)
		  ELSEIF (COORD .EQ. 'BL') THEN
			!DISC
			DX = (LON(N_CTRD_AG+N_CTRD_EVG+I)-LON(J))*PI*ER
     * *COS(0.5*(LAT(N_CTRD_AG+N_CTRD_EVG+I)+LAT(J))*PI)
              DY = (LAT(N_CTRD_AG+N_CTRD_EVG+I)-LAT(J))*ER*PI
	        DISC(J) = SQRT(DX**2+DY**2)
			!DISU
			DX = (LONU(N_CTRD_AG+N_CTRD_EVG+I)-LONU(J))*PI*ER
     * *COS(0.5*(LATU(N_CTRD_AG+N_CTRD_EVG+I)+LATU(J))*PI)
              DY = (LATU(N_CTRD_AG+N_CTRD_EVG+I)-LATU(J))*ER*PI
	        DISU(J) = SQRT(DX**2+DY**2)
			!DISV
			DX = (LONV(N_CTRD_AG+N_CTRD_EVG+I)-LONV(J))*PI*ER
     * *COS(0.5*(LATV(N_CTRD_AG+N_CTRD_EVG+I)+LATV(J))*PI)
              DY = (LATV(N_CTRD_AG+N_CTRD_EVG+I)-LATV(J))*ER*PI
	        DISV(J) = SQRT(DX**2+DY**2)
			
		  ELSE
                
		  ENDIF
		  
		ENDDO
		
		ITEMP = MINLOC(DISC)
		NAGC_IVG(1,I) = ITEMP(1)
		ITEMP = MINLOC(DISU)
		NAGU_IVG(1,I) = ITEMP(1)
		ITEMP = MINLOC(DISV)
		NAGV_IVG(1,I) = ITEMP(1)
		
		
        ENDDO
		
        
!**********************************************************************
! INVERSE DISTANCE WEIGHTING INTERPOLATION
#elif defined UD_IDW
	!IF (TRIM(UDSCH).EQ.'IDW') THEN
	  
	  IF (COORD .EQ. 'XY') THEN
      
      DO I = 1,N_CTRD_IVG
        CALL WEIGHT_IDW(XR(N_CTRD_AG+N_CTRD_EVG+I)
     * ,YR(N_CTRD_AG+N_CTRD_EVG+I),XR,YR,WTC_IVG(:,I),NAGC_IVG(:,I))
        CALL WEIGHT_IDW(XU(N_CTRD_AG+N_CTRD_EVG+I)
     * ,YU(N_CTRD_AG+N_CTRD_EVG+I),XU,YU,WTU_IVG(:,I),NAGU_IVG(:,I))
        CALL WEIGHT_IDW(XV(N_CTRD_AG+N_CTRD_EVG+I)
     * ,YV(N_CTRD_AG+N_CTRD_EVG+I),XV,YV,WTV_IVG(:,I),NAGV_IVG(:,I))
      ENDDO
      
        ELSE
        
      DO I = 1,N_CTRD_IVG
        CALL WEIGHT_IDW(LON(N_CTRD_AG+N_CTRD_EVG+I)
     * ,LAT(N_CTRD_AG+N_CTRD_EVG+I),LON,LAT,WTC_IVG(:,I),NAGC_IVG(:,I))
        CALL WEIGHT_IDW(LONU(N_CTRD_AG+N_CTRD_EVG+I)
     * ,LATU(N_CTRD_AG+N_CTRD_EVG+I),LONU,LATU
     * ,WTU_IVG(:,I),NAGU_IVG(:,I))
        CALL WEIGHT_IDW(LONV(N_CTRD_AG+N_CTRD_EVG+I)
     * ,LATV(N_CTRD_AG+N_CTRD_EVG+I),LONV,LATV
     * ,WTV_IVG(:,I),NAGV_IVG(:,I))
      ENDDO
        
        ENDIF
        
	
!**********************************************************************
! AREA-AVERAGING
#elif defined UD_AVE
	
	DO I = 1,N_CTRD_IVG
		
	  DO K = 1,4
	    XB(K) = XNODE(K,N_CTRD_AG+N_CTRD_EVG+I)
	    YB(K) = YNODE(K,N_CTRD_AG+N_CTRD_EVG+I)
	  ENDDO
	  
	  !NAGC
	  II = 0
	  DO J = 1,N_CTRD_AG
	    IF (COORD.EQ.'XY') THEN
	  	  CALL INSIDE(XR(J),YR(J),XB,YB,NB,IND)
	    ELSE
	  	  CALL INSIDE(LON(J),LAT(J),XB,YB,NB,IND)
	    ENDIF
	    
	    IF (IND.EQ.1) THEN
	  		II = II+1
	  		NAGC_IVG(II,I) = J
				WTC_IVG(II,I) = DJ(J)
	    ENDIF
	  ENDDO
	  
	  !NAGU
	  II = 1
	  NREF = NAGC_IVG(1,I)
	  NAGU_IVG(II,I) = NREF
        WTU_IVG(II,I) = H2(NREF)
	  
	  DO WHILE (.TRUE.)
			NREF = NJP1(NREF)
			IF (COORD.EQ.'XY') THEN
	  			CALL INSIDE(XR(NREF),YR(NREF),XB,YB,NB,IND)
				ELSE
	  			CALL INSIDE(LON(NREF),LAT(NREF),XB,YB,NB,IND)
			ENDIF
		
			IF (IND.EQ.1) THEN
				II = II+1
				NAGU_IVG(II,I) = NREF
				WTU_IVG(II,I) = H2(NREF)
			ELSE
				EXIT
			ENDIF
	  
	  ENDDO
	  
	  !NAGV
	  II = 1
	  NREF = NAGC_IVG(1,I)
	  NAGV_IVG(II,I) = NREF
        WTV_IVG(II,I) = H1(NREF)
	  
	  DO WHILE (.TRUE.)
		NREF = NIP1(NREF)
		IF (COORD.EQ.'XY') THEN
	  	  CALL INSIDE(XR(NREF),YR(NREF),XB,YB,NB,IND)
	    ELSE
	  	  CALL INSIDE(LON(NREF),LAT(NREF),XB,YB,NB,IND)
		ENDIF
		
		IF (IND.EQ.1) THEN
		  II = II+1
		  NAGV_IVG(II,I) = NREF
		  WTV_IVG(II,I) = H1(NREF)
		ELSE
		  EXIT
		ENDIF
	  
	  ENDDO
		
      ENDDO

      
!**********************************************************************
!9-POINT SHAPIRO FILTERING
#elif defined UD_SF
	DO I = 1,N_CTRD_IVG
	  DO J = 1,N_CTRD_AG
		  
		  IF (COORD .EQ. 'XY') THEN
			DISC(J) = SQRT((XR(N_CTRD_AG+N_CTRD_EVG+I)-XR(J))**2
     * +(YR(N_CTRD_AG+N_CTRD_EVG+I)-YR(J))**2)
			DISU(J) = SQRT((XU(N_CTRD_AG+N_CTRD_EVG+I)-XU(J))**2
     * +(YU(N_CTRD_AG+N_CTRD_EVG+I)-YU(J))**2)
			DISV(J) = SQRT((XV(N_CTRD_AG+N_CTRD_EVG+I)-XV(J))**2
     * +(YV(N_CTRD_AG+N_CTRD_EVG+I)-YV(J))**2)
		  ELSEIF (COORD .EQ. 'BL') THEN
			!DISC
			DX = (LON(N_CTRD_AG+N_CTRD_EVG+I)-LON(J))*PI*ER
     * *COS(0.5*(LAT(N_CTRD_AG+N_CTRD_EVG+I)+LAT(J))*PI)
              DY = (LAT(N_CTRD_AG+N_CTRD_EVG+I)-LAT(J))*ER*PI
	        DISC(J) = SQRT(DX**2+DY**2)
			!DISU
			DX = (LONU(N_CTRD_AG+N_CTRD_EVG+I)-LONU(J))*PI*ER
     * *COS(0.5*(LATU(N_CTRD_AG+N_CTRD_EVG+I)+LATU(J))*PI)
              DY = (LATU(N_CTRD_AG+N_CTRD_EVG+I)-LATU(J))*ER*PI
	        DISU(J) = SQRT(DX**2+DY**2)
			!DISV
			DX = (LONV(N_CTRD_AG+N_CTRD_EVG+I)-LONV(J))*PI*ER
     * *COS(0.5*(LATV(N_CTRD_AG+N_CTRD_EVG+I)+LATV(J))*PI)
              DY = (LATV(N_CTRD_AG+N_CTRD_EVG+I)-LATV(J))*ER*PI
	        DISV(J) = SQRT(DX**2+DY**2)
			
		  ELSE
                
		  ENDIF
		  
		ENDDO
			
		!NAGC/WTC
		ITEMP = MINLOC(DISC)
		NAGC_IVG(1,I) = ITEMP(1)
		WTC_IVG(1,I) = (SE-1)**2
		NAGC_IVG(2,I) = NIM1(NAGC_IVG(1,I))
		WTC_IVG(2,I) = SE/2.*(1-SE)
		NAGC_IVG(3,I) = NIP1(NAGC_IVG(1,I))
		WTC_IVG(3,I) = SE/2.*(1-SE)
		NAGC_IVG(4,I) = NJM1(NAGC_IVG(1,I))
		WTC_IVG(4,I) = SE/2.*(1-SE)
		NAGC_IVG(5,I) = NJP1(NAGC_IVG(1,I))
		WTC_IVG(5,I) = SE/2.*(1-SE)
		NAGC_IVG(6,I) = NIM1JP1(NAGC_IVG(1,I))
		WTC_IVG(6,I) = SE**2/4.
		NAGC_IVG(7,I) = NIP1JP1(NAGC_IVG(1,I))
		WTC_IVG(7,I) = SE**2/4.
		NAGC_IVG(8,I) = NIP1JM1(NAGC_IVG(1,I))
		WTC_IVG(8,I) = SE**2/4.
		NAGC_IVG(9,I) = NIM1JM1(NAGC_IVG(1,I))
		WTC_IVG(9,I) = SE**2/4.
		
		!NAGU/WTU
		ITEMP = MINLOC(DISU)
		NAGU_IVG(1,I) = ITEMP(1)
		WTU_IVG(1,I) = 1.-SE
		NAGU_IVG(2,I) = NJM1(NAGU_IVG(1,I))
		WTU_IVG(2,I) = SE/2.
		NAGU_IVG(3,I) = NJP1(NAGU_IVG(1,I))
		WTU_IVG(3,I) = SE/2.
		
		!NAGV/WTV
		ITEMP = MINLOC(DISV)
		NAGV_IVG(1,I) = ITEMP(1)
		WTV_IVG(1,I) = 1.-SE
		NAGV_IVG(2,I) = NIM1(NAGV_IVG(1,I))
		WTV_IVG(2,I) = SE/2.
		NAGV_IVG(3,I) = NIP1(NAGV_IVG(1,I))
		WTV_IVG(3,I) = SE/2.
		
      ENDDO
	
      
!**********************************************************************
!FULL-WEIGHTING OPERATOR (25 POINT)
#elif defined UD_FWO
      !FIND THE NEAREST GRID
	 DO I = 1,N_CTRD_IVG
		DO J = 1,N_CTRD_AG
		  
		  IF (COORD .EQ. 'XY') THEN
			DISC(J) = SQRT((XR(N_CTRD_AG+N_CTRD_EVG+I)-XR(J))**2
     * +(YR(N_CTRD_AG+N_CTRD_EVG+I)-YR(J))**2)
			DISU(J) = SQRT((XU(N_CTRD_AG+N_CTRD_EVG+I)-XU(J))**2
     * +(YU(N_CTRD_AG+N_CTRD_EVG+I)-YU(J))**2)
			DISV(J) = SQRT((XV(N_CTRD_AG+N_CTRD_EVG+I)-XV(J))**2
     * +(YV(N_CTRD_AG+N_CTRD_EVG+I)-YV(J))**2)
		  ELSEIF (COORD .EQ. 'BL') THEN
			!DISC
			DX = (LON(N_CTRD_AG+N_CTRD_EVG+I)-LON(J))*PI*ER
     * *COS(0.5*(LAT(N_CTRD_AG+N_CTRD_EVG+I)+LAT(J))*PI)
              DY = (LAT(N_CTRD_AG+N_CTRD_EVG+I)-LAT(J))*ER*PI
	        DISC(J) = SQRT(DX**2+DY**2)
			!DISU
			DX = (LONU(N_CTRD_AG+N_CTRD_EVG+I)-LONU(J))*PI*ER
     * *COS(0.5*(LATU(N_CTRD_AG+N_CTRD_EVG+I)+LATU(J))*PI)
              DY = (LATU(N_CTRD_AG+N_CTRD_EVG+I)-LATU(J))*ER*PI
	        DISU(J) = SQRT(DX**2+DY**2)
			!DISV
			DX = (LONV(N_CTRD_AG+N_CTRD_EVG+I)-LONV(J))*PI*ER
     * *COS(0.5*(LATV(N_CTRD_AG+N_CTRD_EVG+I)+LATV(J))*PI)
              DY = (LATV(N_CTRD_AG+N_CTRD_EVG+I)-LATV(J))*ER*PI
	        DISV(J) = SQRT(DX**2+DY**2)
			
		  ELSE
                
		  ENDIF
		  
		ENDDO
			
		!NAGC/WTC
		ITEMP = MINLOC(DISC)
		NAGC_IVG(1,I) = ITEMP(1)
		WTC_IVG(1,I) = (1-SE)**2
		NAGC_IVG(2,I) = NIM1(NAGC_IVG(1,I))
		WTC_IVG(2,I) = SE/3.*(1-SE)
		NAGC_IVG(3,I) = NIP1(NAGC_IVG(1,I))
		WTC_IVG(3,I) = SE/3.*(1-SE)
		NAGC_IVG(4,I) = NJM1(NAGC_IVG(1,I))
		WTC_IVG(4,I) = SE/3.*(1-SE)
		NAGC_IVG(5,I) = NJP1(NAGC_IVG(1,I))
		WTC_IVG(5,I) = SE/3.*(1-SE)
		NAGC_IVG(6,I) = NIM1JP1(NAGC_IVG(1,I))
		WTC_IVG(6,I) = SE**2/9.
		NAGC_IVG(7,I) = NIP1JP1(NAGC_IVG(1,I))
		WTC_IVG(7,I) = SE**2/9.
		NAGC_IVG(8,I) = NIP1JM1(NAGC_IVG(1,I))
		WTC_IVG(8,I) = SE**2/9.
		NAGC_IVG(9,I) = NIM1JM1(NAGC_IVG(1,I))
		WTC_IVG(9,I) = SE**2/9.
		NAGC_IVG(10,I) = NJM2(NAGC_IVG(1,I))
		WTC_IVG(10,I) = SE/6.*(1-SE)
		NAGC_IVG(11,I) = NJP2(NAGC_IVG(1,I))
		WTC_IVG(11,I) = SE/6.*(1-SE)
		NAGC_IVG(12,I) = NIM2(NAGC_IVG(1,I))
		WTC_IVG(12,I) = SE/6.*(1-SE)
		NAGC_IVG(13,I) = NIP2(NAGC_IVG(1,I))
		WTC_IVG(13,I) = SE/6.*(1-SE)
		NAGC_IVG(14,I) = NJM2(NIM1(NAGC_IVG(1,I)))
		WTC_IVG(14,I) = SE**2/18.
		NAGC_IVG(15,I) = NJP2(NIM1(NAGC_IVG(1,I)))
		WTC_IVG(15,I) = SE**2/18.
		NAGC_IVG(16,I) = NJM2(NIP1(NAGC_IVG(1,I)))
		WTC_IVG(16,I) = SE**2/18.
		NAGC_IVG(17,I) = NJP2(NIP1(NAGC_IVG(1,I)))
		WTC_IVG(17,I) = SE**2/18.
		NAGC_IVG(18,I) = NJM1(NIM2(NAGC_IVG(1,I)))
		WTC_IVG(18,I) = SE**2/18.
		NAGC_IVG(19,I) = NJP1(NIM2(NAGC_IVG(1,I)))
		WTC_IVG(19,I) = SE**2/18.
		NAGC_IVG(20,I) = NJM1(NIP2(NAGC_IVG(1,I)))
		WTC_IVG(20,I) = SE**2/18.
		NAGC_IVG(21,I) = NJP1(NIP2(NAGC_IVG(1,I)))
		WTC_IVG(21,I) = SE**2/18.
		NAGC_IVG(22,I) = NJM2(NIM2(NAGC_IVG(1,I)))
		WTC_IVG(22,I) = SE**2/36.
		NAGC_IVG(23,I) = NJP2(NIM2(NAGC_IVG(1,I)))
		WTC_IVG(23,I) = SE**2/36.
		NAGC_IVG(24,I) = NJM2(NIP2(NAGC_IVG(1,I)))
		WTC_IVG(24,I) = SE**2/36.
		NAGC_IVG(25,I) = NJP2(NIP2(NAGC_IVG(1,I)))
		WTC_IVG(25,I) = SE**2/36.
		
		!NAGU/WTU
		ITEMP = MINLOC(DISU)
		NAGU_IVG(1,I) = ITEMP(1)
		WTU_IVG(1,I) = 1.-SE
		NAGU_IVG(2,I) = NJM1(NAGU_IVG(1,I))
		WTU_IVG(2,I) = SE/3.
		NAGU_IVG(3,I) = NJP1(NAGU_IVG(1,I))
		WTU_IVG(3,I) = SE/3.
		NAGU_IVG(4,I) = NJM2(NAGU_IVG(1,I))
		WTU_IVG(4,I) = SE/6.
		NAGU_IVG(5,I) = NJP2(NAGU_IVG(1,I))
		WTU_IVG(5,I) = SE/6.
		
		!NAGV/WTV
		ITEMP = MINLOC(DISV)
		NAGV_IVG(1,I) = ITEMP(1)
		WTV_IVG(1,I) = 1.-SE
		NAGV_IVG(2,I) = NIM1(NAGV_IVG(1,I))
		WTV_IVG(2,I) = SE/3.
		NAGV_IVG(3,I) = NIP1(NAGV_IVG(1,I))
		WTV_IVG(3,I) = SE/3.
		NAGV_IVG(4,I) = NIM2(NAGV_IVG(1,I))
		WTV_IVG(4,I) = SE/6.
		NAGV_IVG(5,I) = NIP2(NAGV_IVG(1,I))
		WTV_IVG(5,I) = SE/6.
		
	ENDDO

!**********************************************************************
	
	
#else
      PRINT*, 'Update scheme undefined!'
	PAUSE
	STOP
#endif
	
!======================================================================
!							 END: UPDATE
!======================================================================
	
	PRINT*, 'DONE!'
	
!AVE
!      OPEN (100,FILE='NAGC_CHK_AVE.DAT')
!	DO I = 1,N_CTRD_IVG
!        WRITE (100,1113) N_CTRD_AG+N_CTRD_EVG+I,NAGC_IVG(:,I)
!        WRITE (100,1114) N_CTRD_AG+N_CTRD_EVG+I,WTC_IVG(:,I)
!      ENDDO
!      CLOSE (100)
!	
!	OPEN (100,FILE='NAGU_CHK_AVE.DAT')
!	DO I = 1,N_CTRD_IVG
!        WRITE (100,1113) N_CTRD_AG+N_CTRD_EVG+I,NAGU_IVG(:,I)
!        WRITE (100,1114) N_CTRD_AG+N_CTRD_EVG+I,WTU_IVG(:,I)
!      ENDDO
!      CLOSE (100)
!	
!	OPEN (100,FILE='NAGV_CHK_AVE.DAT')
!	DO I = 1,N_CTRD_IVG
!        WRITE (100,1113) N_CTRD_AG+N_CTRD_EVG+I,NAGV_IVG(:,I)
!        WRITE (100,1114) N_CTRD_AG+N_CTRD_EVG+I,WTV_IVG(:,I)
!      ENDDO
!      CLOSE (100)
!	
!      PRINT*, "AVE NAG CHECK DONE!"
!      PAUSE 
!      STOP
!	
!1113  FORMAT (I6,100I7) 
!1114  FORMAT (I6,100F12.5)
      
1111	FORMAT (I6,100I7)   
1112  FORMAT (I6,100F14.6) 
		
	RETURN
	
	END SUBROUTINE PRE_IE
	
	
	SUBROUTINE FIND_NEAREST_AG (XINT,YINT,X,Y,N_AG)
	
	USE MOD_GLOBAL
	
	IMPLICIT NONE
	
	INTEGER I
	INTEGER N_AG
	INTEGER ITEMP(1)
	REAL, PARAMETER :: PI = 3.1415926/180.
	REAL XINT,YINT,DX,DY
	REAL X(N_CTRD_AG),Y(N_CTRD_AG),DIS(N_CTRD_AG)
	
	DO I = 1,N_CTRD_AG
		  
	  IF (COORD .EQ. 'XY') THEN
		DIS(I) = SQRT((XINT-X(I))**2+(YINT-Y(I))**2)
	  ELSEIF (COORD .EQ. 'BL') THEN
		DX = (XINT-X(I))*PI*ER*COS(0.5*(YINT+Y(I))*PI)
		DY = (YINT-Y(I))*ER*PI
		DIS(I) = SQRT(DX**2+DY**2)
	  ENDIF
		  
	ENDDO
	
	ITEMP = MINLOC(DIS)
	N_AG = ITEMP(1)
	
	RETURN
	END SUBROUTINE FIND_NEAREST_AG
	

!      !CHECK
!      
!      OPEN (100,FILE='WEIGHT_CHK_VG.DAT')
!      DO I = 1,N_CTRD_VG
!        WRITE (100,1111) N_CTRD_AG+I,NAGREF_C(:,I),WEIREF_C(:,I)
!      ENDDO
!      CLOSE (100)
!      
!      OPEN (100,FILE='WEIGHT_CHK_VGU.DAT')
!      DO I = 1,N_CTRD_VG
!        WRITE (100,1111) N_CTRD_AG+I,NAGREF_U(:,I),WEIREF_U(:,I)
!      ENDDO
!      CLOSE (100)
!      
!      OPEN (100,FILE='WEIGHT_CHK_VGV.DAT')
!      DO I = 1,N_CTRD_VG
!        WRITE (100,1111) N_CTRD_AG+I,NAGREF_V(:,I),WEIREF_V(:,I)
!      ENDDO
!      CLOSE (100)
!      
!      OPEN (100,FILE='XYUV.DAT')
!      DO I = N_CTRD_AG+1,N_CTRD
!        WRITE (100,*) I,XU(I),YU(I),XV(I),YV(I)
!      ENDDO
!      CLOSE (100)
!      
!      PRINT*, "WEIGHT CHECK DONE!"
!      PAUSE 
!      STOP
!      
!1111  FORMAT (I6,4I7,4F14.6)   
	
! DR
!      OPEN (100,FILE='NAG_CHK_DR.DAT')
!      DO I = 1,N_CTRD_EVG
!        WRITE (100,1112) N_CTRD_AG+I,NAGC_EVG(1,I)
!	ENDDO
!	DO I = 1,N_CTRD_IVG
!        WRITE (100,1112) N_CTRD_AG+N_CTRD_EVG+I,NAGC_IVG(1,I)
!      ENDDO
!      CLOSE (100)
!	
!      PRINT*, "DR NAG CHECK DONE!"
!      PAUSE 
!      STOP
!1112	FORMAT (I6,I7) 
	
!AVE
!      OPEN (100,FILE='NAGC_CHK_AVE.DAT')
!	DO I = 1,N_CTRD_IVG
!        WRITE (100,1113) N_CTRD_AG+N_CTRD_EVG+I,NAGC_IVG(:,I)
!      ENDDO
!      CLOSE (100)
!	
!	OPEN (100,FILE='NAGU_CHK_AVE.DAT')
!	DO I = 1,N_CTRD_IVG
!        WRITE (100,1113) N_CTRD_AG+N_CTRD_EVG+I,NAGU_IVG(:,I)
!      ENDDO
!      CLOSE (100)
!	
!	OPEN (100,FILE='NAGV_CHK_AVE.DAT')
!	DO I = 1,N_CTRD_IVG
!        WRITE (100,1113) N_CTRD_AG+N_CTRD_EVG+I,NAGV_IVG(:,I)
!      ENDDO
!      CLOSE (100)
!	
!      PRINT*, "AVE NAG CHECK DONE!"
!      PAUSE 
!      STOP
!	
!1113	FORMAT (I6,100I7) 
	
!      !SHAPIRO/FULL-WEIGHTING
!      OPEN (100,FILE='WEIGHT_CHK_IVGC.DAT')
!      DO I = 1,N_CTRD_IVG
!        WRITE (100,1111) N_CTRD_AG+N_CTRD_EVG+I,NAGC_IVG(:,I)
!	  WRITE (100,1112) N_CTRD_AG+N_CTRD_EVG+I,WTC_IVG(:,I)
!      ENDDO
!      CLOSE (100)
!      
!      OPEN (100,FILE='WEIGHT_CHK_IVGU.DAT')
!      DO I = 1,N_CTRD_IVG
!        WRITE (100,1111) N_CTRD_AG+N_CTRD_EVG+I,NAGU_IVG(:,I)
!	  WRITE (100,1112) N_CTRD_AG+N_CTRD_EVG+I,WTU_IVG(:,I)
!      ENDDO
!      CLOSE (100)
!      
!      OPEN (100,FILE='WEIGHT_CHK_IVGV.DAT')
!      DO I = 1,N_CTRD_IVG
!        WRITE (100,1111) N_CTRD_AG+N_CTRD_EVG+I,NAGV_IVG(:,I)
!	  WRITE (100,1112) N_CTRD_AG+N_CTRD_EVG+I,WTV_IVG(:,I)
!      ENDDO
!      CLOSE (100)
!      
!      
!      PRINT*, "WEIGHT CHECK DONE!"
!      PAUSE 
!      STOP
!      
!1111	FORMAT (I6,100I7) 
!1112  FORMAT (I6,100F14.6)
	
!!      !QUADRATIC/LINEAR/UAEI/HPI
!      OPEN (100,FILE='WEIGHT_CHK_EVGC.DAT')
!      DO I = 1,N_CTRD_EVG
!        WRITE (100,1111) N_CTRD_AG+I,NAGC_EVG(:,I)
!	  WRITE (100,1112) N_CTRD_AG+I,WTC_EVG(:,I)
!      ENDDO
!      CLOSE (100)
!      
!      OPEN (100,FILE='WEIGHT_CHK_EVGU.DAT')
!      DO I = 1,N_CTRD_EVG
!        WRITE (100,1111) N_CTRD_AG+I,NAGU_EVG(:,I)
!	  WRITE (100,1112) N_CTRD_AG+I,WTU_EVG(:,I)
!      ENDDO
!      CLOSE (100)
!      
!      OPEN (100,FILE='WEIGHT_CHK_EVGV.DAT')
!      DO I = 1,N_CTRD_EVG
!        WRITE (100,1111) N_CTRD_AG+I,NAGV_EVG(:,I)
!	  WRITE (100,1112) N_CTRD_AG+I,WTV_EVG(:,I)
!      ENDDO
!      CLOSE (100)
!      
!      
!      PRINT*, "WEIGHT CHECK DONE!"
!      PAUSE 
!      STOP
!      
!1111	FORMAT (I6,100I7)   
!1112  FORMAT (I6,100F14.6) 