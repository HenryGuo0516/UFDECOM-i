#include "DEFS.h"
	
	MODULE MOD_OUTPUT
	SAVE
	CONTAINS

	SUBROUTINE OUTPUT_FIELD
      
	USE MOD_GLOBAL
      
	CHARACTER*80 FNAME,F1
	INTEGER IUFPT,IUVPT,IUEPT,IUMMT,IUBED	
      
      REAL, ALLOCATABLE :: ACCEU(:,:)
	REAL, ALLOCATABLE :: ACCEV(:,:)
	REAL, ALLOCATABLE :: PGFU(:,:)
	REAL, ALLOCATABLE :: PGFV(:,:)
	REAL, ALLOCATABLE :: PGBPU(:)
	REAL, ALLOCATABLE :: PGBPV(:)
	REAL, ALLOCATABLE :: PGBCU(:,:)
	REAL, ALLOCATABLE :: PGBCV(:,:)
	REAL, ALLOCATABLE :: FRICU(:,:)
	REAL, ALLOCATABLE :: FRICV(:,:)
	REAL, ALLOCATABLE :: CORIU(:,:)
	REAL, ALLOCATABLE :: CORIV(:,:)
	REAL, ALLOCATABLE :: ADV_U(:,:)
	REAL, ALLOCATABLE :: ADV_V(:,:)
	REAL, ALLOCATABLE :: ADVW_U(:,:)
	REAL, ALLOCATABLE :: ADVW_V(:,:)	

	REAL, ALLOCATABLE :: ACCEUR(:,:)
	REAL, ALLOCATABLE :: ACCEVR(:,:)
	REAL, ALLOCATABLE :: PGFUR(:,:)
	REAL, ALLOCATABLE :: PGFVR(:,:)
	REAL, ALLOCATABLE :: PGBPUR(:)
	REAL, ALLOCATABLE :: PGBPVR(:)
	REAL, ALLOCATABLE :: PGBCUR(:,:)
	REAL, ALLOCATABLE :: PGBCVR(:,:)
	REAL, ALLOCATABLE :: FRICUR(:,:)
	REAL, ALLOCATABLE :: FRICVR(:,:)
	REAL, ALLOCATABLE :: CORIUR(:,:)
	REAL, ALLOCATABLE :: CORIVR(:,:)
	REAL, ALLOCATABLE :: ADV_UR(:,:)
	REAL, ALLOCATABLE :: ADV_VR(:,:)
	REAL, ALLOCATABLE :: ADVW_UR(:,:)
	REAL, ALLOCATABLE :: ADVW_VR(:,:)
      
	IUSPT = 394
	IUVPT = 395
	IUEPT = 396
      IUBED = 397
	IUMMT = 398

!*=======================================================================
!*               OUTPUT [on/off] list for [GPU to CPU]
!*=======================================================================
!* Turn on/off the field output
![on]#define FPT
![on]#ifdef FPT
!* # define FSMPT //turn on/off fsm field output
![on]# define VPT //turn on/off velocity field output
![on]# define EPT //turn on/off elevation field output
![on]# define SPT //turn on/off salinity field output
!* # define MPT //turn on/off KM and KH field output
![on]# define TPT //turn on/off temperature field output
!* # define SEDPT //turn on/off sediment field output
!*#define BEDPT //turn on/off sediment bed field output
![on]#define MMTPT //turn on/off momentum field output
!#endif      
c!$ACC UPDATE HOST (T,UR,VR,WR,EL) !GPU to CPU fo Sub OUTPUT_FIELD-SEC
!$ACC UPDATE HOST (T,S,UR,VR,WR,EL) !GPU to CPU fo Sub OUTPUT_FIELD
!$ACC UPDATE HOST (FSM,DUM,DVM,U,V,UNN,VNN) !GPU to CPU fo Sub OUTPUT_FIELD
!$ACC UPDATE HOST (DRHOX,DRHOY,DU,DV) !GPU to CPU fo Sub OUTPUT_FIELD
c!$ACC UPDATE HOST (DU,DV) !GPU to CPU fo Sub OUTPUT_FIELD-BAROPG5B
!$ACC UPDATE HOST (KM,WUSURF,WVSURF,WUBOT,WVBOT) !GPU to CPU fo Sub OUTPUT_FIELD
c!$ACC UPDATE HOST (KM,WUBOT,WVBOT) !GPU to CPU fo Sub OUTPUT_FIELD-BCOND(9)
      
      
	WRITE (F1,1000) NN_FPT
1000	FORMAT (I6.6)

#ifdef FSMPT
      !IF (LOG_FSMPT) THEN
        FNAME="./"//TRIM(OUT_DIRE)//TRIM(XG)//"/FIELD_DISTRI/FSM/
     *FSM_FIELD_"//TRIM(F1)
	  OPEN (IUSPT,FILE=FNAME,FORM='UNFORMATTED',CONVERT='BIG_ENDIAN')
	  WRITE (IUSPT) THOUR,FSM
	  CLOSE (IUSPT)
	!ENDIF
#endif
  
#ifdef SPT
	!IF (LOG_SPT) THEN
	  FNAME="./"//TRIM(OUT_DIRE)//TRIM(XG)
     *//"/FIELD_DISTRI/SALINITY/S_FIELD_"//TRIM(F1)
c	  OPEN (IUSPT,FILE=FNAME) !c V2410
c	  WRITE (IUSPT,*) THOUR
c        DO I=1,N_CTRD
c            WRITE (IUSPT,9876) (S(I,K),K=1,10)
c        ENDDO    
	  OPEN (IUSPT,FILE=FNAME,FORM='UNFORMATTED',CONVERT='BIG_ENDIAN') !V2410
	  WRITE (IUSPT) THOUR,S
	  CLOSE (IUSPT)
	!ENDIF
#endif
9876  FORMAT (10F12.4)
#ifdef TPT
	!IF (LOG_SPT) THEN
	  FNAME="./"//TRIM(OUT_DIRE)//TRIM(XG)
     *//"/FIELD_DISTRI/TEMPERATURE/T_FIELD_"//TRIM(F1)
	  OPEN (IUSPT,FILE=FNAME,FORM='UNFORMATTED',CONVERT='BIG_ENDIAN')
	  WRITE (IUSPT) THOUR,T
	  CLOSE (IUSPT)
	!ENDIF
#endif

#ifdef MPT
!$ACC UPDATE HOST (CONCENTRATION) !GPU to CPU V2410
	  FNAME="./"//TRIM(OUT_DIRE)//TRIM(XG)
     *//"/FIELD_DISTRI/MATRERIAL/M_FIELD_"//TRIM(F1)
	  OPEN (IUSPT,FILE=FNAME,FORM='UNFORMATTED',CONVERT='BIG_ENDIAN')
c	  WRITE (IUSPT) THOUR,KM,KH !c V2410
	  WRITE (IUSPT) THOUR,CONCENTRATION !V2410
	  CLOSE (IUSPT)
#endif

#ifdef SEDPT
	!IF (LOG_SEDPT) THEN
	  FNAME="./"//TRIM(OUT_DIRE)//"/"//"/FIELD_DISTRI/SEDIMENTS/SED_FIELD_"
     *//TRIM(F1)
	  OPEN (IUSPT,FILE=FNAME,FORM='UNFORMATTED',CONVERT='BIG_ENDIAN')
	  WRITE (IUSPT) THOUR,SED
	  CLOSE (IUSPT)
	!ENDIF
#endif

#ifdef VPT
	!IF (LOG_VPT) THEN
	  FNAME="./"//TRIM(OUT_DIRE)//TRIM(XG)
     *//"/FIELD_DISTRI/CURRENT/V_FIELD_"//TRIM(F1)
	  OPEN (IUVPT,FILE=FNAME,FORM='UNFORMATTED',CONVERT='BIG_ENDIAN')
	  WRITE (IUVPT) THOUR,UR,VR,WR 
	  CLOSE (IUVPT)
	!ENDIF
#endif

#ifdef EPT
	!IF (LOG_EPT) THEN
      FNAME="./"//TRIM(OUT_DIRE)//TRIM(XG)
     *//"/FIELD_DISTRI/ELEVATION/EL_FIELD_"//TRIM(F1)
c	OPEN (IUEPT,FILE=FNAME) !c V2410
c	WRITE (IUEPT,*) THOUR
c      DO I=1,N_CTRD
c            WRITE (IUEPT,9876) EL(I)
c      ENDDO    
	OPEN (IUEPT,FILE=FNAME,FORM='UNFORMATTED',CONVERT='BIG_ENDIAN') !V2410
	WRITE (IUEPT) THOUR,EL
	CLOSE (IUEPT)
	!ENDIF
#endif

#ifdef BEDPT
	!IF (LOG_BEDPT) THEN                                                   
        FNAME=
     * "./"//TRIM(OUT_DIRE)//"/"//"/FIELD_DISTRI/BED/ZB_FIELD_"         
     *//TRIM(F1)
        OPEN (IUBED,FILE=FNAME,FORM='UNFORMATTED',CONVERT='BIG_ENDIAN')
        WRITE (IUBED) THOUR,TAU,TAU_WAVE,TAU_TIDE,QERO,QDEP,ZBED,ZBEDD
	  CLOSE (IUBED)
	!ENDIF
#endif     
      

	
C================================MOMENTUM==============================
	
! CALCULATING EACH PART OF MOMENTUM EQUATION BELOW (AT GRID CENTER)

#ifdef MMTPT
      !IF (LOG_MMTPT) THEN
        
      FNAME=
     * "./"//TRIM(OUT_DIRE)//"/"//"/FIELD_DISTRI/MOMENTUM/MMT_FIELD_"
     *//TRIM(F1)
	OPEN (IUMMT,FILE=FNAME,FORM='UNFORMATTED',CONVERT='BIG_ENDIAN')
	WRITE (IUMMT) THOUR
!----------------------------------------------------------------------
! ±º‰œÓ
      ALLOCATE (ACCEU(N_CTRDP1,KB),ACCEV(N_CTRDP1,KB))
      ALLOCATE (ACCEUR(N_CTRD_AG,KB),ACCEVR(N_CTRD_AG,KB))     
      ACCEU = 0
      ACCEV = 0
      ACCEUR = 0
      ACCEVR = 0
      DO I = 1,N_CTRD
      IF (DUM(I).EQ.1) THEN
      DO K = 1,KBM1
            ACCEU(I,K) = (U(I,K)-UNN(I,K))/DTI
      ENDDO											   
      ENDIF
      ENDDO
      DO I = 1,N_CTRD
      IF (DVM(I).EQ.1) THEN
	    DO K = 1,KBM1
              ACCEV(I,K) = (V(I,K)-VNN(I,K))/DTI
          ENDDO
      ENDIF
      ENDDO
      DO K = 1,KBM1
      DO I = 1,N_CTRD_AG
      IF (FSM(I)*FSM(NIM1(I))*FSM(NIP1(I))
     * *FSM(NJM1(I))*FSM(NJP1(I)).GT.0.5) THEN
		! CALCULATE VARIABLE AT THE CENTER OF GRID
          UIJK = (ACCEU(I,K)+ACCEU(NIP1(I),K))/2.
          VIJK = (ACCEV(I,K)+ACCEV(NJP1(I),K))/2.
     	    IF (FSM(NIP1(I)).EQ.0..AND.UIJK.GT.0) UIJK=0.
	    IF (FSM(NIM1(I)).EQ.0..AND.UIJK.LT.0) UIJK=0.
     	    IF (FSM(NJP1(I)).EQ.0..AND.VIJK.GT.0) VIJK=0.
	    IF (FSM(NJM1(I)).EQ.0..AND.VIJK.LT.0) VIJK=0.
		! TRANSFER VARIABLE TO CARTESIAN COORDINATES
          ACCEUR(I,K) = YETA(I)*UIJK/H2(I)-YXI(I)*VIJK/H1(I)
          ACCEVR(I,K) = -XETA(I)*UIJK/H2(I)+XXI(I)*VIJK/H1(I)
      ELSE 
	    ACCEUR(I,K) = 0.
	    ACCEVR(I,K) = 0.
      ENDIF
      ENDDO
      ENDDO
	WRITE (IUMMT) ACCEUR,ACCEVR
      DEALLOCATE (ACCEU,ACCEV)
      DEALLOCATE (ACCEUR,ACCEVR)
!----------------------------------------------------------------------


!----------------------------------------------------------------------
! ADV_UR, ADV_VR
! ADVECTION TERMS
      ALLOCATE (ADV_U(N_CTRDP1,KB),ADV_V(N_CTRDP1,KB))
      ALLOCATE (ADV_UR(N_CTRD_AG,KB),ADV_VR(N_CTRD_AG,KB))
      
      ADV_U = 0
      ADV_V = 0
      ADV_UR = 0
      ADV_VR = 0
      
      DO I = 1,N_CTRD
      IF (DUM(I).EQ.1) THEN
	DO K = 1,KBM2
      ADV_U(I,K) = 
     * U(I,K)*(U(NIP1(I),K)-U(NIM1(I),K))/(H1(I)+H1(NIM1(I)))
     * +0.25*(V(NIM1(I),K)+V(I,K)+V(NIM1JP1(I),K)+V(NJP1(I),K))
     * *(U(NJP1(I),K)-U(NJM1(I),K))
     * /(H2(I)+0.5*H2(NJM1(I))+0.5*H2(NJP1(I)))
     * +0.25*(WR(I,K)+WR(NIM1(I),K)+WR(I,K+1)+WR(NIM1(I),K+1))
     * *(U(I,K)-U(I,K+1))/(DZZ(K)*DU(I))
      ENDDO
      K = KBM1
      ADV_U(I,K) = 
     * U(I,K)*(U(NIP1(I),K)-U(NIM1(I),K))/(H1(I)+H1(NIM1(I)))
     * +0.25*(V(NIM1(I),K)+V(I,K)+V(NIM1JP1(I),K)+V(NJP1(I),K))
     * *(U(NJP1(I),K)-U(NJM1(I),K))
     * /(H2(I)+0.5*H2(NJM1(I))+0.5*H2(NJP1(I)))
      ENDIF
      ENDDO
      DO I = 1,N_CTRD
      IF (DVM(I).EQ.1) THEN 
      DO K = 1,KBM2
          ADV_V(I,K) =
     *  0.25*(U(I,K)+U(NIP1(I),K)+U(NJM1(I),K)+U(NIP1JM1(I),K))
     *  *(V(NIP1(I),K)-V(NIM1(I),K))
     *  /(H1(I)+0.5*H1(NIM1(I))+0.5*H1(NIP1(I)))
     *  +V(I,K)*(V(NJP1(I),K)-V(NJM1(I),K))/(H2(I)+H2(NJM1(I)))
     *  +0.25*(WR(I,K)+WR(NJM1(I),K)+WR(I,K+1)+WR(NJM1(I),K+1))
     *  *(V(I,K)-V(I,K+1))/(DZZ(K)*DV(I))
      ENDDO
      K = KBM1
      ADV_V(I,K) =
     * 0.25*(U(I,K)+U(NIP1(I),K)+U(NJM1(I),K)+U(NIP1JM1(I),K))
     * *(V(NIP1(I),K)-V(NIM1(I),K))
     * /(H1(I)+0.5*H1(NIM1(I))+0.5*H1(NIP1(I)))
     * +V(I,K)*(V(NJP1(I),K)-V(NJM1(I),K))/(H2(I)+H2(NJM1(I)))
      ENDIF
      ENDDO
      
      DO K = 1,KBM1
      DO I = 1,N_CTRD_AG
      IF (FSM(I)*FSM(NIM1(I))*FSM(NIP1(I))
     * *FSM(NJM1(I))*FSM(NJP1(I)).GT.0.5) THEN
        UIJK = (ADV_U(I,K)+ADV_U(NIP1(I),K))/2.
        VIJK = (ADV_V(I,K)+ADV_V(NJP1(I),K))/2.
     	  IF (FSM(NIP1(I)).EQ.0..AND.UIJK.GT.0) UIJK = 0.
	  IF (FSM(NIM1(I)).EQ.0..AND.UIJK.LT.0) UIJK = 0.
     	  IF (FSM(NJP1(I)).EQ.0..AND.VIJK.GT.0) VIJK = 0.
	  IF (FSM(NJM1(I)).EQ.0..AND.VIJK.LT.0) VIJK = 0.
        ADV_UR(I,K) = YETA(I)*UIJK/H2(I)-YXI(I)*VIJK/H1(I)
        ADV_VR(I,K) = -XETA(I)*UIJK/H2(I)+XXI(I)*VIJK/H1(I)
	ELSE 
	  ADV_UR(I,K) = 0.
	  ADV_VR(I,K) = 0.
      ENDIF
      ENDDO
      ENDDO
	WRITE (IUMMT) -1.*ADV_UR,-1.*ADV_VR
      DEALLOCATE (ADV_U,ADV_V)
      DEALLOCATE (ADV_UR,ADV_VR)
!----------------------------------------------------------------------

      
!----------------------------------------------------------------------
! CORIUR, CORIVR
! CORIOLIS FORCE
      ALLOCATE (CORIU(N_CTRDP1,KB),CORIV(N_CTRDP1,KB))
      ALLOCATE (CORIUR(N_CTRD_AG,KB),CORIVR(N_CTRD_AG,KB))
      
      CORIU = 0
      CORIV = 0
      CORIUR = 0
      CORIVR = 0
      DO I = 1,N_CTRD
      IF (DUM(I).EQ.1) THEN
      DO K = 1,KBM1
            CORIU(I,K) = 
     *-0.25*COR(I)*(V(I,K)+V(NJP1(I),K)+V(NIM1JP1(I),K)+V(NIM1(I),K))
      ENDDO
      ENDIF
      ENDDO
       
      DO I = 1,N_CTRD
      IF (DVM(I).EQ.1) THEN
      DO K = 1,KBM1
      CORIV(I,K) = 
     *0.25*COR(I)*(U(I,K)+U(NJM1(I),K)+U(NIP1(I),K)+U(NIP1JM1(I),K))
      ENDDO
      ENDIF
      ENDDO
      
      
      DO K = 1,KBM1
      DO I = 1,N_CTRD_AG
      IF (FSM(I)*FSM(NIM1(I))*FSM(NIP1(I))
     * *FSM(NJM1(I))*FSM(NJP1(I)).GT.0.5) THEN
          UIJK = (CORIU(I,K)+CORIU(NIP1(I),K))/2.
          VIJK = (CORIV(I,K)+CORIV(NJP1(I),K))/2.
     		IF (FSM(NIP1(I)).EQ.0..AND.UIJK.GT.0) UIJK = 0.
		IF (FSM(NIM1(I)).EQ.0..AND.UIJK.LT.0) UIJK = 0.
     		IF (FSM(NJP1(I)).EQ.0..AND.VIJK.GT.0) VIJK = 0.
		IF (FSM(NJM1(I)).EQ.0..AND.VIJK.LT.0) VIJK = 0.
          CORIUR(I,K) = YETA(I)*UIJK/H2(I)-YXI(I)*VIJK/H1(I)
          CORIVR(I,K) = -XETA(I)*UIJK/H2(I)+XXI(I)*VIJK/H1(I)
      ELSE 
		  CORIUR(I,K) = 0.
		  CORIVR(I,K) = 0.
      ENDIF
      ENDDO
      ENDDO
      
	WRITE (IUMMT) -1.*CORIUR,-1.*CORIVR
      
      DEALLOCATE (CORIU,CORIV)
      DEALLOCATE (CORIUR,CORIVR)
!----------------------------------------------------------------------


!----------------------------------------------------------------------
! PGBPUR, PGBPVR
! BAROTROPIC PRESSURE GRADIENT FORCE
      
      ALLOCATE (PGBPU(N_CTRDP1),PGBPV(N_CTRDP1))
      ALLOCATE (PGBPUR(N_CTRD_AG),PGBPVR(N_CTRD_AG))
      
      PGBPU = 0
      PGBPV = 0
      PGBPUR = 0
      PGBPVR = 0
      
      DO I = 1,N_CTRD
      IF (DUM(I).EQ.1) THEN
      PGBPU(I) = -1.*GRAV*(EL(I)-EL(NIM1(I)))/(0.5*(H1(I)+H1(NIM1(I))))
      ENDIF
      ENDDO
      
      DO I = 1,N_CTRD
      IF (DVM(I).EQ.1) THEN
      PGBPV(I) = -1.*GRAV*(EL(I)-EL(NJM1(I)))/(0.5*(H2(I)+H2(NJM1(I))))
      ENDIF
      ENDDO
      
      DO I = 1,N_CTRD_AG
      IF (FSM(I)*FSM(NIM1(I))*FSM(NIP1(I))
     * *FSM(NJM1(I))*FSM(NJP1(I)).GT.0.5) THEN
      UIJK = (PGBPU(I)+PGBPU(NIP1(I)))/2.
      VIJK = (PGBPV(I)+PGBPV(NJP1(I)))/2.
      IF (FSM(NIP1(I)).EQ.0..AND.UIJK.GT.0) UIJK = 0.
      IF (FSM(NIM1(I)).EQ.0..AND.UIJK.LT.0) UIJK = 0.
      IF (FSM(NJP1(I)).EQ.0..AND.VIJK.GT.0) VIJK = 0.
      IF (FSM(NJM1(I)).EQ.0..AND.VIJK.LT.0) VIJK = 0.
        PGBPUR(I) = YETA(I)*UIJK/H2(I)-YXI(I)*VIJK/H1(I)
        PGBPVR(I) = -XETA(I)*UIJK/H2(I)+XXI(I)*VIJK/H1(I)
      ELSE                                               
	  PGBPUR(I) = 0.
	  PGBPVR(I) = 0.
      ENDIF
	ENDDO
	
      WRITE (IUMMT) PGBPUR,PGBPVR
	
      DEALLOCATE (PGBPU,PGBPV)
      DEALLOCATE (PGBPUR,PGBPVR)
!----------------------------------------------------------------------


!----------------------------------------------------------------------
! PGBCUR, PGBCVR
! BAROCLINIC PRESSURE GRADIENT FORCE

      ALLOCATE (PGBCU(N_CTRDP1,KB),PGBCV(N_CTRDP1,KB))
      ALLOCATE (PGBCUR(N_CTRD_AG,KB),PGBCVR(N_CTRD_AG,KB))
      
      PGBCU = 0
      PGBCV = 0
      PGBCUR = 0
      PGBCVR = 0
      
      DO I = 1,N_CTRD
        IF (DUM(I).EQ.1) THEN
	    DO K = 1,KBM1
            PGBCU(I,K) = DRHOX(I,K)/(DU(I)*0.5*(DJ(I)+DJ(NIM1(I))))
          ENDDO
        ENDIF
      ENDDO
      
      DO I = 1,N_CTRD
        IF (DVM(I).EQ.1) THEN
	    DO K = 1,KBM1
            PGBCV(I,K) = DRHOY(I,K)/(DV(I)*0.5*(DJ(I)+DJ(NJM1(I))))
          ENDDO
        ENDIF
      ENDDO
      
      
      DO K = 1,KBM1
      DO I = 1,N_CTRD_AG
      IF (FSM(I)*FSM(NIM1(I))*FSM(NIP1(I))
     * *FSM(NJM1(I))*FSM(NJP1(I)).GT.0.5) THEN
        UIJK = (PGBCU(I,K)+PGBCU(NIP1(I),K))/2.
        VIJK = (PGBCV(I,K)+PGBCV(NJP1(I),K))/2.
     	  IF (FSM(NIP1(I)).EQ.0..AND.UIJK.GT.0) UIJK = 0.
	  IF (FSM(NIM1(I)).EQ.0..AND.UIJK.LT.0) UIJK = 0.
     	  IF (FSM(NJP1(I)).EQ.0..AND.VIJK.GT.0) VIJK = 0.
	  IF (FSM(NJM1(I)).EQ.0..AND.VIJK.LT.0) VIJK = 0.
        PGBCUR(I,K) = YETA(I)*UIJK/H2(I)-YXI(I)*VIJK/H1(I)
        PGBCVR(I,K) = -XETA(I)*UIJK/H2(I)+XXI(I)*VIJK/H1(I)
	ELSE 
	  PGBCUR(I,K) = 0.
	  PGBCVR(I,K) = 0.
	ENDIF
	ENDDO
      ENDDO
      
	WRITE (IUMMT) -1.*PGBCUR,-1.*PGBCVR
	
      DEALLOCATE (PGBCU,PGBCV)
      DEALLOCATE (PGBCUR,PGBCVR)
!----------------------------------------------------------------------


!----------------------------------------------------------------------
! FRICUR, FRICVR
! FRICTION
      ALLOCATE (FRICU(N_CTRDP1,KB),FRICV(N_CTRDP1,KB))
      ALLOCATE (FRICUR(N_CTRD_AG,KB),FRICVR(N_CTRD_AG,KB))
      
      FRICU = 0
      FRICV = 0
      FRICUR = 0
      FRICVR = 0
          
      DO I = 1,N_CTRD
      IF (DUM(I).EQ.1) THEN
	DO K = 1,KBM1
      IF (K.EQ.1) THEN
      FRICU(I,K) = 
     * (-1.*WUSURF(I)-KM(I,K+1)*(U(I,K)-U(I,K+1))/(DZZ(K)*DU(I)))
     * /(DZ(K)*DU(I))
      
      ELSEIF (K.EQ.KBM1) THEN
      FRICU(I,K) = 
     * (KM(I,K)*(U(I,K-1)-U(I,K))/(DZZ(K-1)*DU(I))-(-1.)*WUBOT(I))
     * /(DZ(K)*DU(I))

      ELSEIF (K.GT.1.AND.K.LT.KBM1) THEN
      FRICU(I,K) = (KM(I,K)*(U(I,K-1)-U(I,K))/DZZ(K-1)
     * -KM(I,K+1)*(U(I,K)-U(I,K+1))/DZZ(K))/DZ(K)/DU(I)**2
      ENDIF
      ENDDO
      ENDIF
      ENDDO
      
      DO I = 1,N_CTRD
      IF (DVM(I).EQ.1) THEN
	DO K = 1,KBM1
      IF (K.EQ.1) THEN
      FRICV(I,K)=
     * (-1.*WVSURF(I)-KM(I,K+1)*(V(I,K)-V(I,K+1))/(DZZ(K)*DV(I)))
     * /(DZ(K)*DV(I))
   
      ELSEIF (K.EQ.KBM1) THEN
      FRICV(I,K)=
     * (KM(I,K)*(V(I,K-1)-V(I,K))/(DZZ(K-1)*DV(I))-(-1.)*WVBOT(I))
     * /(DZ(K)*DV(I))

      ELSEIF (K.GT.1.AND.K.LT.KBM1) THEN
      FRICV(I,K)=
     * (KM(I,K)*(V(I,K-1)-V(I,K))/DZZ(K-1)-KM(I,K+1)*(V(I,K)-V(I,K+1))
     * /DZZ(K))/DZ(K)/DV(I)**2
      ENDIF
      ENDDO
      ENDIF
      ENDDO
       
      DO K = 1,KBM1
      DO I = 1,N_CTRD_AG
      IF (FSM(I)*FSM(NIM1(I))*FSM(NIP1(I))
     * *FSM(NJM1(I))*FSM(NJP1(I)).GT.0.5) THEN
        UIJK = (FRICU(I,K)+FRICU(NIP1(I),K))/2.
        VIJK = (FRICV(I,K)+FRICV(NJP1(I),K))/2.
     	  IF (FSM(NIP1(I)).EQ.0..AND.UIJK.GT.0) UIJK = 0.
	  IF (FSM(NIM1(I)).EQ.0..AND.UIJK.LT.0) UIJK = 0.
     	  IF (FSM(NJP1(I)).EQ.0..AND.VIJK.GT.0) VIJK = 0.
	  IF (FSM(NJM1(I)).EQ.0..AND.VIJK.LT.0) VIJK = 0.
        FRICUR(I,K) = YETA(I)*UIJK/H2(I)-YXI(I)*VIJK/H1(I)
        FRICVR(I,K) = -XETA(I)*UIJK/H2(I)+XXI(I)*VIJK/H1(I)
      ELSE 
        FRICUR(I,K) = 0.
        FRICVR(I,K) = 0.
      ENDIF
      ENDDO
      ENDDO
      
	WRITE (IUMMT) FRICUR,FRICVR
      
      DEALLOCATE (FRICU,FRICV)
      DEALLOCATE (FRICUR,FRICVR)
!----------------------------------------------------------------------

	CLOSE (IUMMT)
	!ENDIF
#endif

	
	RETURN
	END SUBROUTINE OUTPUT_FIELD
	
      
      
	
	SUBROUTINE OUTPUT_POINT
	USE MOD_GLOBAL
      
c      REAL U_T(N_CTRD_AG),V_T(N_CTRD_AG),S_T(N_CTRD_AG),SED_T(N_CTRD_AG)
c      REAL T_T(N_CTRD_AG)
c      REAL XMFL_T(N_CTRD_AG),YMFL_T(N_CTRD_AG) !LXY
	REAL UR_E(KBM1),VR_E(KBM1),XMFL_E(KBM1),YMFL_E(KBM1)  !LXY
	REAL UV(KBM1),DIRE(KBM1)
	REAL S_E(KBM1),T_E(KBM1)
	REAL SED_E(KBM1)
	REAL RE,REEL,REDEP,TIMES,RETAU_WAVE,RETAU_TIDE,RETAU,RQERO,RQDEP
	REAL RETAUD,RETAUE,REHSIG,RETSIG,REWAVEDIR,REUBM,REWINDU,REWINDV
      REAL RE8,REZBEDD
      REAL REKM(KBM1),REKH(KBM1)
	REAL NUM
	CHARACTER*100 VELO(50),DIRE1(50),SALI(50),TEMP(50),FNN,SEDI(50)
      CHARACTER*100 XMFL(50),YMFL(50)    !LXY

      
      !INITIAL VALUE
	IF ( (RESTAR.EQ.'cold'.AND.NSTEP.EQ.ISTART) .OR.
     *(RESTAR.EQ.'hot'.AND.NSTEP.EQ.ISTART.AND.THOUR_HOT.EQ.0.) ) THEN
      
cc      S_T = 0.0
cc      T_T = 0.0
cc      SED_T = 0.0
c      UR_E = 0.0
c      VR_E = 0.0
c      S_E = 0.0
c      T_E = 0.0
      SED_E = 0.0
      REEL = 0.0
      REDEP = 0.0
      RETAU_WAVE = 0.0
      RETAU_TIDE = 0.0
      RETAU = 0.0
      RQERO = 0.0
      RQDEP = 0.0
      RETAUD = 0.0
      RETAUE = 0.0
      REHSIG = 0.0
      RETSIG = 0.0
      REWAVEDIR = 0.0
      REUBM = 0.0
      REWINDU = 0.0
      REWINDV = 0.0
      REZBEDD = 0.0
c      REKM = 0.0
c      REKH = 0.0

      ENDIF !INITIAL VALUE
      SED_E = 0.0
      
c	PI = 180./3.141592653589
	NIT = 400
      
c!$ACC UPDATE HOST (UR,VR,T,KM,KH,EL) !GPU to CPU fo Sub OUTPUT_POINT-SEC
!$ACC UPDATE HOST (UR,VR,T,S,KM,KH,EL) !GPU to CPU fo Sub OUTPUT_POINT

	DO II = 1,N_TSR
	X0 = X_TSR(II)
	Y0 = Y_TSR(II)
      N0 = NN_TSR(II)
	NIT = NIT+1
      
	DO K = 1,KBM1
          
c ---------- CBR: cancelled heavy loops:          
!	DO I = 1,N_CTRD_AG
      
!      U_T(I) = UR(I,K)
!	V_T(I) = VR(I,K)
      
!      XMFL_T(I) = XMFLUX(I,K)
!      YMFL_T(I) = YMFLUX(I,K)
      
!#ifdef MODULE_SAL
!	  S_T(I) = S(I,K)
!#endif
!#ifdef MODULE_TMP
!      T_T(I) = T(I,K)
!#endif
!#ifdef MODULE_SED
!	  SED_T(I) = SED(I,K)
!#endif
!      ENDDO
c      CALL INTERP_POINT(UR_E(K),U_T,FSM(1:N_CTRD_AG),NAGC_TSR(:,II)
c     * ,WTC_TSR(:,II))
      CALL INTERP_POINT3D(UR_E(K),UR,K,FSM,NAGC_TSR(:,II)
     * ,WTC_TSR(:,II))
c      CALL INTERP_POINT(VR_E(K),V_T,FSM(1:N_CTRD_AG),NAGC_TSR(:,II)
c     * ,WTC_TSR(:,II))
      CALL INTERP_POINT3D(VR_E(K),VR,K,FSM,NAGC_TSR(:,II)
     * ,WTC_TSR(:,II))
C************ U AND DIR****************
	UV(K) = SQRT(UR_E(K)**2+VR_E(K)**2)
c ------ CBR: Calculate Dir Ver 1:
      DIRE(k)=ATAN2D(UR_E(k),VR_E(k))
      if (DIRE(k)<0.0) DIRE(k)=360.+DIRE(k)
c ------ CBR: Calculate Dir Ver 0:
c	IF(UR_E(K).GT.0.0 .AND. VR_E(K).GT.0.0) 
c     *DIRE(K) = ATAN(UR_E(K)/VR_E(K))*PI
	
c	IF(UR_E(K).GT.0.0 .AND. VR_E(K).LT.0.0) 
c     *DIRE(K) = 180.-ATAN(UR_E(K)/ABS(VR_E(K)))*PI

c	IF(UR_E(K).LT.0.0 .AND. VR_E(K).LT.0.0) 
c     *DIRE(K) = 180.+ATAN(ABS(UR_E(K))/ABS(VR_E(K)))*PI

c	IF(UR_E(K).LT.0.0 .AND. VR_E(K).GT.0.0) 
c     *DIRE(K) = 360.-ATAN(ABS(UR_E(K))/VR_E(K))*PI

c	IF(UR_E(K).EQ.0.0 .AND. VR_E(K).GE.0.0) DIRE(K) = 0.0
c	IF(UR_E(K).EQ.0.0 .AND. VR_E(K).LT.0.0) DIRE(K) = 180.0
c	IF(UR_E(K).LT.0.0 .AND. VR_E(K).EQ.0.0) DIRE(K) = 270.0
c	IF(UR_E(K).GT.0.0 .AND. VR_E(K).EQ.0.0) DIRE(K) = 90.0

      !SAL
#ifdef MODULE_SAL
	  !CALL INTER_POINT(N0,X0,Y0,S_T,RE)
c        CALL INTERP_POINT(S_E(K),S_T,FSM(1:N_CTRD_AG),NAGC_TSR(:,II)
c     * ,WTC_TSR(:,II))
        CALL INTERP_POINT3D(S_E(K),S,K,FSM,NAGC_TSR(:,II)
     * ,WTC_TSR(:,II))
#endif

#ifdef MODULE_TMP
c      CALL INTERP_POINT(T_E(K),T_T,FSM(1:N_CTRD_AG),NAGC_TSR(:,II)
c     * ,WTC_TSR(:,II))
      CALL INTERP_POINT3D(T_E(K),T,K,FSM,NAGC_TSR(:,II)
     * ,WTC_TSR(:,II))
#endif

      !SED
#ifdef MODULE_SED
	  !CALL INTER_POINT(N0,X0,Y0,SED_T,RE)
   !     SED_E(K)=RE
c        CALL INTERP_POINT(SED_E(K),SED_T,FSM(1:N_CTRD_AG),NAGC_TSR(:,II)
c     * ,WTC_TSR(:,II))
        CALL INTERP_POINT3D(SED_E(K),SED,K,FSM,NAGC_TSR(:,II)
     * ,WTC_TSR(:,II))
#endif
      ENDDO
      DO K = 1,KBM1
c        CALL INTERP_POINT(REKM(K),KM(1:N_CTRD_AG,K),FSM(1:N_CTRD_AG)
c     * ,NAGC_TSR(:,II),WTC_TSR(:,II))
        CALL INTERP_POINT3D(REKM(K),KM,K,FSM
     * ,NAGC_TSR(:,II),WTC_TSR(:,II))
        REKM(K) = 1.E4*REKM(K)
c        CALL INTERP_POINT(REKH(K),KH(1:N_CTRD_AG,K),FSM
c     * ,NAGC_TSR(:,II),WTC_TSR(:,II))
        CALL INTERP_POINT3D(REKH(K),KH,K,FSM
     * ,NAGC_TSR(:,II),WTC_TSR(:,II))
        REKH(K) = 1.E4*REKH(K)
      ENDDO
      CALL INTERP_POINT(REEL,EL,FSM
     * ,NAGC_TSR(:,II),WTC_TSR(:,II))
#ifdef MODULE_SED
      CALL INTERP_POINT(RETAUE,TAUE,FSM
     * ,NAGC_TSR(:,II),WTC_TSR(:,II))
	
	CALL INTERP_POINT(RETAUD,TAUD,FSM
     * ,NAGC_TSR(:,II),WTC_TSR(:,II))

	CALL INTERP_POINT(RETAU,TAU,FSM
     * ,NAGC_TSR(:,II),WTC_TSR(:,II))

	CALL INTERP_POINT(RETAU_WAVE,TAU_WAVE,FSM
     * ,NAGC_TSR(:,II),WTC_TSR(:,II))

	CALL INTERP_POINT(RETAU_TIDE,TAU_TIDE,FSM
     * ,NAGC_TSR(:,II),WTC_TSR(:,II))

	CALL INTERP_POINT(RQERO,QERO,FSM
     * ,NAGC_TSR(:,II),WTC_TSR(:,II))

	CALL INTERP_POINT(RQDEP,QDEP,FSM
     * ,NAGC_TSR(:,II),WTC_TSR(:,II))

	CALL INTERP_POINT(REHSIG,HSIG,FSM
     * ,NAGC_TSR(:,II),WTC_TSR(:,II))

	CALL INTERP_POINT(RETSIG,TSIG,FSM
     * ,NAGC_TSR(:,II),WTC_TSR(:,II))

	CALL INTERP_POINT(REWAVEDIR,WAVEDIR,FSM
     * ,NAGC_TSR(:,II),WTC_TSR(:,II))

	CALL INTERP_POINT(REUBM,UBM,FSM
     * ,NAGC_TSR(:,II),WTC_TSR(:,II))

	CALL INTERP_POINT(REZBEDD,ZBEDD,FSM
     * ,NAGC_TSR(:,II),WTC_TSR(:,II))

	CALL INTERP_POINT(REWINDU,WINDU,FSM
     * ,NAGC_TSR(:,II),WTC_TSR(:,II))

	CALL INTERP_POINT(REWINDV,WINDV,FSM
     * ,NAGC_TSR(:,II),WTC_TSR(:,II))

	CALL INTERP_POINT(REUBM_WB,UBM_WB,FSM
     * ,NAGC_TSR(:,II),WTC_TSR(:,II))

!      ENDIF  
#endif
      
	IF ( (RESTAR.EQ.'cold'.AND.NSTEP.EQ.ISTART) .OR.
     *(RESTAR.EQ.'hot'.AND.NSTEP.EQ.ISTART.AND.THOUR_HOT.EQ.0.) ) THEN
      
C      CALL INTER_POINT(N0,X0,Y0,H,RE)
C	REDEP = RE
      CALL INTERP_POINT(REDEP,H,FSMADD
     * ,NAGC_TSR(:,II),WTC_TSR(:,II))
      
	WRITE (NIT,6001) REDEP
	DO K = 1,KBM1
	WRITE (VELO(K),5001) K
	WRITE (DIRE1(K),5002) K
	WRITE (SALI(K),5003) K
      WRITE (TEMP(K),5004) K
	WRITE (SEDI(K),5005) K
      WRITE (XMFL(K),5008) K
	WRITE (YMFL(K),5009) K	
	ENDDO
	WRITE (NIT,'(200A16)') 'TIME(HOUR)','ELEVATION',(TRIM(VELO(K))
     *,TRIM(DIRE1(K)),K=1,KBM1),(TRIM(SALI(K)),K=1,KBM1)
     *,(TRIM(TEMP(K)),K=1,KBM1),(TRIM(SEDI(K)),K=1,KBM1),'TAU',
     *'TAU_WAVE','TAU_TIDE','TAUE','TAUD','QERO','QDEP',
     *'HSIG','TSIG','WAVEDIR','ZBEDD', !LXY
     *'WIND_U10','WIND_V10','UBM_CALC',
     *'KM02E4','KH02E4','KM09E4', 'KH09E4' !LXY      
C     *,(TRIM(XMFL(K)),TRIM(YMFL(K)),K=1,KBM1)  !LXY
	
      ENDIF
	
      
5001	FORMAT('VELOCITY',I2.2)
5002	FORMAT ('DIRECTION',I2.2)
5003  FORMAT('SALINITY',I2.2)
5004  FORMAT('TEMPERATURE',I2.2)
5005  FORMAT('SEDIMENT',I2.2)
5008  FORMAT('XMFLUX',I2.2)
5009	FORMAT('YMFLUX',I2.2)     
      
C      CALL INTER_POINT(N0,X0,Y0,EL,RE)
C	REEL = RE
      	
	WRITE (NIT,6000)  THOUR,REEL,(UV(K),DIRE(K),K=1,KBM1),
     * (S_E(K),K=1,KBM1),(T_E(K),K=1,KBM1),
     * (SED_E(K),K=1,KBM1),RETAU,RETAU_WAVE,
     * RETAU_TIDE,RETAUE,RETAUD,RQERO,RQDEP,
     * REHSIG,RETSIG,REWAVEDIR,REZBEDD, !LXY
     * REWINDU,REWINDV,REUBM_WB,
     * REKM(2),REKH(2),REKM(KBM2),REKH(KBM2)  !LXY
C     * ,(XMFL_E(K),YMFL_E(K),K=1,KBM1)  !LXY
	NUMPOINT = NUMPOINT+1
	ENDDO

6000	FORMAT (200F16.6)
6001	FORMAT (' STEADY_DEPTH  ',F8.2)
	RETURN
	END SUBROUTINE OUTPUT_POINT
	
	
      
      SUBROUTINE OUTPUT_SECFLUX
	USE MOD_GLOBAL

	NNN = 400+N_TSR+10
      
	DO II=1,N_SEC_N
      NIT=NNN+II
      WRITE(NIT,6001) THOUR,FLUXSUM(II),FLUX_ACM(II),SFLUXSUM(II),   !LXY
     *    SFLUX_ACM(II),SEDFLUXSUM(II),SEDFLUX_ACM(II)
      ENDDO

6001  FORMAT (F14.2,F16.0,F26.0,F16.0,F26.0,F16.0,F26.0)
!7004	FORMAT (F14.2,F16.0,F16.0,F16.0)

	RETURN
	END SUBROUTINE OUTPUT_SECFLUX
      
      

      SUBROUTINE RESIDUAL
	USE MOD_GLOBAL
      DIMENSION UA(N_CTRD_AG),CA(N_CTRD_AG),VA(N_CTRD_AG)
      DIMENSION UHH(N_CTRD_AG),VHH(N_CTRD_AG),CHH(N_CTRD_AG)
      DIMENSION UCHH(N_CTRD_AG),VCHH(N_CTRD_AG)
      DIMENSION UCDH(N_CTRD_AG),VCDH(N_CTRD_AG)
      CHARACTER*20 RE_U(100),RE_V(100),REX_U(100),REX_V(100),RESX_U(100)
     *,RESX_V(100),RES(100)
      REAL DTT
      
      NIT1=2000
	NIT2=2200
	NIT3=2400
	NIT4=2600
	NIT5=2800
      DTT=DTI/3600.
      DO 320 JHIST=1,JHM
      IF (THOUR.GE.HIST(JHIST,1).AND.THOUR.LE.HIST(JHIST,2)) THEN
!$ACC UPDATE HOST (UR,VR,S,D,FSM,EL) !GPU to CPU for RESIDUAL
c!$ACC UPDATE HOST (UR,VR,D,FSM,EL) !GPU to CPU for RESIDUAL-SEC
      DO I=1,N_CTRD_AG
      IF (FSM(I).GT.0.) THEN
          UA_TEMP=0.
		CA_TEMP=0.
          VA_TEMP=0.
          
		UCDH_TEMP=0.
		VCDH_TEMP=0.
          DO K=1,KBM1
              UA_TEMP=UA_TEMP+UR(I,K)*DZ(K)
		    CA_TEMP=CA_TEMP+S(I,K)*DZ(K)
		    VA_TEMP=VA_TEMP+VR(I,K)*DZ(K)
          ENDDO
          UA(I)=UA_TEMP
          VA(I)=VA_TEMP
          CA(I)=CA_TEMP
          UHH(I)=UA_TEMP*D(I)
          VHH(I)=VA_TEMP*D(I)
          CHH(I)=CA_TEMP*D(I)
          UCHH(I)=UA_TEMP*D(I)*CA_TEMP
          VCHH(I)=VA_TEMP*D(I)*CA_TEMP
          DO K=1,KBM1
              UD_TEMP=UR(I,K)-UA(I)
		    VD_TEMP=VR(I,K)-VA(I)
		    CD_TEMP=S(I,K)-CA(I)
		    UCDH_TEMP=UCDH_TEMP+UD_TEMP*CD_TEMP*DZ(K)
		    VCDH_TEMP=VCDH_TEMP+VD_TEMP*CD_TEMP*DZ(K)
          ENDDO
          UCDH(I)=UCDH_TEMP*D(I)
		VCDH(I)=VCDH_TEMP*D(I)
      ELSE
          UA(I)=0.
		CA(I)=0.
		VA(I)=0.
		UHH(I)=0.
		CHH(I)=0.
		UCHH(I)=0.
		VHH(I)=0.
		VCHH(I)=0.
		UCDH(I)=0.
		VCDH(I)=0.
      ENDIF
      ENDDO
      DO I=1,N_CTRD_AG
          ARC_UA(I,JHIST)=ARC_UA(I,JHIST)+UA(I)*DTT*DEI(JHIST)
	    ARC_CA(I,JHIST)=ARC_CA(I,JHIST)+CA(I)*DTT*DEI(JHIST)
	    ARC_VA(I,JHIST)=ARC_VA(I,JHIST)+VA(I)*DTT*DEI(JHIST)
	    ARC_UHH(I,JHIST)=ARC_UHH(I,JHIST)+UHH(I)*DTT*DEI(JHIST)
	    ARC_CHH(I,JHIST)=ARC_CHH(I,JHIST)+CHH(I)*DTT*DEI(JHIST)
	    ARC_UCHH(I,JHIST)=ARC_UCHH(I,JHIST)+UCHH(I)*DTT*DEI(JHIST)
	    ARC_VHH(I,JHIST)=ARC_VHH(I,JHIST)+VHH(I)*DTT*DEI(JHIST)
	    ARC_VCHH(I,JHIST)=ARC_VCHH(I,JHIST)+VCHH(I)*DTT*DEI(JHIST)
	    ARC_UCDH(I,JHIST)=ARC_UCDH(I,JHIST)+UCDH(I)*DTT*DEI(JHIST)
	    ARC_VCDH(I,JHIST)=ARC_VCDH(I,JHIST)+VCDH(I)*DTT*DEI(JHIST)
	    ARC_H(I,JHIST)=ARC_H(I,JHIST)+D(I)*DTT*DEI(JHIST)
      ENDDO
      DO K=1,KBM1
      DO I=1,N_CTRD_AG
      IF (FSM(I).GT.0.0) THEN
          ARCUR(I,K,JHIST)=ARCUR(I,K,JHIST)+UR(I,K)*DTT
	    ARCVR(I,K,JHIST)=ARCVR(I,K,JHIST)+VR(I,K)*DTT
	    ARCUX(I,K,JHIST)=ARCUX(I,K,JHIST)+UR(I,K)*DTT*D(I)*DZ(K)
	    ARCVX(I,K,JHIST)=ARCVX(I,K,JHIST)+VR(I,K)*DTT*D(I)*DZ(K)
	    ARCSUX(I,K,JHIST)=ARCSUX(I,K,JHIST)+S(I,K)*UR(I,K)*D(I)*DZ(K)*DTT
	    ARCSVX(I,K,JHIST)=ARCSVX(I,K,JHIST)+S(I,K)*VR(I,K)*D(I)*DZ(K)*DTT
          ARCS(I,K,JHIST)=ARCS(I,K,JHIST)+S(I,K)*DTT
      ENDIF
      ENDDO
      ENDDO
      DO I=1,N_CTRD_AG
          IF (FSM(I).GT.0.0) THEN
              DTT22(I,JHIST)=DTT22(I,JHIST)+DTT
              ARCET(I,JHIST)=ARCET(I,JHIST)+EL(I)*DTT*DEI(JHIST)
          ENDIF
      ENDDO
      ELSEIF(THOUR.GT.HIST(JHIST,2)) THEN
      IF (LOG_ARC(JHIST)) THEN
          
      LOG_ARC(JHIST)=.FALSE.    
      DO I=1,N_CTRD_AG
          FU1_WATER(I,JHIST)=ARC_UA(I,JHIST)*ARC_H(I,JHIST)
          FU2_WATER(I,JHIST)=ARC_UHH(I,JHIST)
     *    -ARC_UA(I,JHIST)*ARC_H(I,JHIST)
		FU1(I,JHIST)=ARC_UA(I,JHIST)*ARC_H(I,JHIST)
     *    *ARC_CA(I,JHIST)
		FU2(I,JHIST)=ARC_UHH(I,JHIST)*ARC_CA(I,JHIST)
     *	-ARC_UA(I,JHIST)*ARC_H(I,JHIST)*ARC_CA(I,JHIST)
		FU3(I,JHIST)=ARC_CHH(I,JHIST)*ARC_UA(I,JHIST)-
     *    ARC_UA(I,JHIST)*ARC_H(I,JHIST)*ARC_CA(I,JHIST)
		FU4(I,JHIST)=ARC_UCHH(I,JHIST)-ARC_UHH(I,JHIST)
     *    *ARC_CA(I,JHIST)-ARC_CHH(I,JHIST)*ARC_UA(I,JHIST)
     *    +ARC_H(I,JHIST)*ARC_UA(I,JHIST)*ARC_CA(I,JHIST)
		FU5(I,JHIST)=ARC_UCDH(I,JHIST)
          
          FV1_WATER(I,JHIST)=ARC_VA(I,JHIST)*ARC_H(I,JHIST)
		FV2_WATER(I,JHIST)=ARC_VHH(I,JHIST)
     *    -ARC_VA(I,JHIST)*ARC_H(I,JHIST)
		FV1(I,JHIST)=ARC_VA(I,JHIST)*ARC_H(I,JHIST)
     *    *ARC_CA(I,JHIST)
		FV2(I,JHIST)=ARC_VHH(I,JHIST)*ARC_CA(I,JHIST)
     *    -ARC_VA(I,JHIST)*ARC_H(I,JHIST)*ARC_CA(I,JHIST)
		FV3(I,JHIST)=ARC_CHH(I,JHIST)*ARC_VA(I,JHIST)
     *    -ARC_VA(I,JHIST)*ARC_H(I,JHIST)*ARC_CA(I,JHIST)
		FV4(I,JHIST)=ARC_VCHH(I,JHIST)-ARC_VHH(I,JHIST)
     *    *ARC_CA(I,JHIST)- ARC_CHH(I,JHIST)*ARC_VA(I,JHIST)+
     *    ARC_H(I,JHIST)*ARC_VA(I,JHIST)*ARC_CA(I,JHIST)
		FV5(I,JHIST)=ARC_VCDH(I,JHIST)
      ENDDO
      DO K=1,KBM1
      DO I=1,N_CTRD_AG
      IF(FSMADD(I).EQ.1..AND.DTT22(I,JHIST).GT.0.) THEN
          DEI22(I,JHIST)=1./DTT22(I,JHIST)
	    ARCUR(I,K,JHIST)=ARCUR(I,K,JHIST)*DEI22(I,JHIST)
	    ARCVR(I,K,JHIST)=ARCVR(I,K,JHIST)*DEI22(I,JHIST)
	    ARCUX(I,K,JHIST)=ARCUX(I,K,JHIST)*DEI22(I,JHIST)
	    ARCVX(I,K,JHIST)=ARCVX(I,K,JHIST)*DEI22(I,JHIST)
	    ARCSUX(I,K,JHIST)=ARCSUX(I,K,JHIST)*DEI22(I,JHIST)
	    ARCSVX(I,K,JHIST)=ARCSVX(I,K,JHIST)*DEI22(I,JHIST)
          ARCS(I,K,JHIST)=ARCS(I,K,JHIST)*DEI22(I,JHIST)
      ENDIF
      ENDDO
      ENDDO
      NIT=NIT1+JHIST
      WRITE (NIT,2000) HIST(JHIST,1),HIST(JHIST,2)
2000  FORMAT ('3D EULERIAN RESIDUAL CURRENT , AVERAGED
     *FROM: HOUR',F10.4,' TO: HOUR',F10.4)
      DO K=1,KBM1
	    WRITE (RE_U(K),"('   RESI_U_',I2.2)") K
	    WRITE (RE_V(K),"('   RESI_V_',I2.2)") K
	ENDDO
      WRITE (NIT,3000) (TRIM(RE_U(K)),TRIM(RE_V(K)),K=1,KBM1)
3000  FORMAT ('     I     J   X-COORD   Y-COORD   RESI_EL',200A12)
      DO I=1,N_CTRD_AG
      IF (FSMADD(I).EQ.1..AND.DTT22(I,JHIST).GT.0.) THEN
	WRITE (NIT,2001) I,J,XR(I),YR(I),ARCET(I,JHIST),
     *(ARCUR(I,K,JHIST),ARCVR(I,K,JHIST),K=1,KBM1)
	ENDIF    
      ENDDO
2001  FORMAT (2I6,2F10.1,F10.2,200F12.2)
      CLOSE (NIT)
      NIT=NIT5+JHIST
	WRITE (NIT,2007) HIST(JHIST,1),HIST(JHIST,2)
2007	FORMAT ('3D RESIDUAL UNIT FLUX, AVERAGED
     *FROM: HOUR',F10.4,' TO: HOUR',F10.4)
	DO K=1,KBM1
	    WRITE (REX_U(K),"('  RESI_UF_',I2.2)") K
	    WRITE (REX_V(K),"('  RESI_VF_',I2.2)") K
	ENDDO        
      WRITE (NIT,3007) (TRIM(REX_U(K)),TRIM(REX_V(K)),K=1,KBM1)        
3007	FORMAT ('     I     J   X-COORD   Y-COORD',200A12)
      DO I=1,N_CTRD_AG
      IF (FSMADD(I).EQ.1..AND.DTT22(I,JHIST).GT.0.) THEN     !QC
	WRITE (NIT,2001) I,J,XR(I),YR(I),
     *(ARCUX(I,K,JHIST),ARCVX(I,K,JHIST),K=1,KBM1)
	ENDIF    
      ENDDO
      CLOSE (NIT)
#if defined MODULE_SAL
      NIT=NIT2+JHIST
	WRITE (NIT,2003) HIST(JHIST,1),HIST(JHIST,2)
2003	FORMAT ('3D RESIDUAL SALT FLUX, AVERAGED
     *FROM: HOUR',F10.4,' TO: HOUR',F10.4)
	DO K=1,KBM1
	    WRITE (RES(K),"(' RESI_SAL_',I2.2)") K	
	    WRITE (RESX_U(K),"(' RESI_USF_',I2.2)") K
	    WRITE (RESX_V(K),"(' RESI_VSF_',I2.2)") K
	ENDDO
	WRITE (NIT,3001) (TRIM(RES(K)),K=1,KBM1),(TRIM(RESX_U(K))
     *,TRIM(RESX_V(K)),K=1,KBM1)
3001	FORMAT ('     I     J   X-COORD   Y-COORD',200A12)
	DO I=1,N_CTRD_AG
      IF (FSMADD(I).EQ.1..AND.DTT22(I,JHIST).GT.0.) THEN
	WRITE (NIT,2001) I,J,XR(I),YR(I),
     *(ARCS(I,K,JHIST),K=1,KBM1)
     *,(ARCSUX(I,K,JHIST),ARCSVX(I,K,JHIST),K=1,KBM1),
     *DTT22(I,JHIST)                                        
	ENDIF
      ENDDO
      CLOSE (NIT)
      
      NIT=NIT4+JHIST
	WRITE (NIT,2005) HIST(JHIST,1),HIST(JHIST,2)
2005	FORMAT ('RESIDUAL SALT FLUX AND ITS MECHANISM ANALYSIS 
     *CONSITUTENTS,AVERAGED FROM: HOUR',F10.4,' TO: HOUR',F10.4)
	WRITE (NIT,'(A224)') '     I     J   X-COORD   Y-COORD
     *  SALTFLUX_U  SALTFLUX_V       LAG_U       LAG_V      ELER_U
     *      ELER_V    SOTKES_U    SOTKES_V  TIDEPUMP_U  TIDEPUMP_V
     * TIDEPUMP1_U TIDEPUMP1_V TIDEPUMP2_U TIDEPUMP2_V VERTSHEAR_U
     * VERTSHEAR_V'
	DO I=1,N_CTRD_AG
	IF (FSMADD(I).EQ.1.) THEN
	AA=FU1(I,JHIST)+FU2(I,JHIST)+FU3(I,JHIST)+FU4(I,JHIST)
     *+FU5(I,JHIST)
	BB=FV1(I,JHIST)+FV2(I,JHIST)+FV3(I,JHIST)+FV4(I,JHIST)
     *+FV5(I,JHIST)
	WRITE (NIT,2001) I,J,XR(I),YR(I),AA,BB,FU1(I,JHIST)+FU2(I,JHIST)
     *,FV1(I,JHIST)+FV2(I,JHIST),FU1(I,JHIST),FV1(I,JHIST),
     *FU2(I,JHIST),FV2(I,JHIST),FU3(I,JHIST)+FU4(I,JHIST)
     *,FV3(I,JHIST)+FV4(I,JHIST),FU3(I,JHIST),FV3(I,JHIST),
     *FU4(I,JHIST),FV4(I,JHIST),FU5(I,JHIST),FV5(I,JHIST)
	ENDIF
	ENDDO
	CLOSE (NIT)
#endif
      ENDIF
      ENDIF
320   CONTINUE
      RETURN
      
      END SUBROUTINE RESIDUAL


	
      END MODULE
