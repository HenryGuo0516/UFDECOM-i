#include "DEFS.h"
	SUBROUTINE	  OPENFILE(IDX)
	USE MOD_GLOBAL
      INTEGER II
	INTEGER IDX
      INTEGER NNN,I0,K0,F0 !V2410
      REAL X0,Y0,Z0,TT,R1,R2,D1,D2,C1,C2,DD1,CC1 !V2410
      REAL DHMIN !V2410
	CHARACTER*80 FNAME,FN1,COMMENT
	CHARACTER*100 MDOUT,MDTSR,MDSEC,MDFPT,MDFPTS,MDFPTUV,MDFPTEL,MDFPTMMT
	CHARACTER*100 MDMDG,MDLAG,MDTSP,MDRSC,MDRSF,MDPLT,MDHAM,MDIMG,MDFPTSED
	CHARACTER*100 MDFPTBED,MDFPTFSM,MDFCK,MDCHK,MDRST,MDFPTT,MDPOI,MDFPTM
	CHARACTER*100 VELO(50),DIRE(50),SALI(50),FNN
	CHARACTER*800 MAKEDIR
	CHARACTER*10 TMP
	LOGICAL ALIVE
	PI = 180./3.1415926535
      SELECT CASE (IDX)
***********************************************************************
*                                                                     *
*                               CASE 1                                *
*                                                                     *
***********************************************************************
      CASE(1)
      IF (TRIM(OSTYPE).EQ.'windows') THEN
          TMP = 'mkdir'
      ELSEIF (TRIM(OSTYPE).EQ.'linux') THEN
          TMP = 'mkdir -p'
      ENDIF

	MDOUT = TRIM(TMP)//" "//TRIM(OUT_DIRE)
	MDTSR = TRIM(OUT_DIRE)//TRIM(XG)//"TIMESERIES"
	MDSEC = TRIM(OUT_DIRE)//TRIM(XG)//"SECFLUX"
	MDLAG = TRIM(OUT_DIRE)//TRIM(XG)//"LAG_TRACKING"
	MDFPT = TRIM(OUT_DIRE)//TRIM(XG)//"FIELD_DISTRI"
      MDFPTFSM = TRIM(OUT_DIRE)//TRIM(XG)//
     *	"FIELD_DISTRI"//TRIM(XG)//"FSM"
	MDFPTS = TRIM(OUT_DIRE)//TRIM(XG)//
     *  "FIELD_DISTRI"//TRIM(XG)//"SALINITY"
      MDFPTM = TRIM(OUT_DIRE)//TRIM(XG)//
     *  "FIELD_DISTRI"//TRIM(XG)//"MATRERIAL"
      MDFPTT = TRIM(OUT_DIRE)//TRIM(XG)//
     *  "FIELD_DISTRI"//TRIM(XG)//"TEMPERATURE"
	MDFPTUV = TRIM(OUT_DIRE)//TRIM(XG)//
     *	"FIELD_DISTRI"//TRIM(XG)//"CURRENT"
	MDFPTEL = TRIM(OUT_DIRE)//TRIM(XG)//
     *	"FIELD_DISTRI"//TRIM(XG)//"ELEVATION"
	MDFPTSED = TRIM(OUT_DIRE)//TRIM(XG)//
     *	"FIELD_DISTRI"//TRIM(XG)//"SEDIMENTS"
	MDFPTBED = TRIM(OUT_DIRE)//TRIM(XG)//
     *	"FIELD_DISTRI"//TRIM(XG)//"BED"
 	MDFPTMMT = TRIM(OUT_DIRE)//TRIM(XG)//
     *	"FIELD_DISTRI"//TRIM(XG)//"MOMENTUM"
	MDTSP = TRIM(OUT_DIRE)//TRIM(XG)//"TIMESTEP"
	MDRSC = TRIM(OUT_DIRE)//TRIM(XG)//"RESI_EULER"
	MDRSF = TRIM(OUT_DIRE)//TRIM(XG)//"RESI_FLUX"
	MDHAM = TRIM(OUT_DIRE)//TRIM(XG)//"HARMONIC"
      MDIMG = TRIM(OUT_DIRE)//TRIM(XG)//"SAVEIMG"
      MDPOI=TRIM(OUT_DIRE)//TRIM(XG)// "POI_ERROR"
	MDMDG = "MODGEN"
      MDFCK = TRIM(MDMDG)//TRIM(XG)//"FIELD_CHECK"
      MDCHK = TRIM(MDMDG)//TRIM(XG)//"VGA_CHECK"
	MDRST = "RSTFILES"
      
      MAKEDIR=TRIM(MDOUT)//" "//TRIM(MDTSR)//" "//
     *TRIM(MDSEC)//" "//TRIM(MDLAG)//" "//TRIM(MDFPTM)//" "//
     *TRIM(MDFPT)//" "//TRIM(MDFPTFSM)//" "//TRIM(MDFPTS)//" "//
     *TRIM(MDFPTT)//" "//TRIM(MDFPTUV)//" "//TRIM(MDFPTEL)//" "//
     *TRIM(MDTSP)//" "//TRIM(MDRSC)//" "//
     *TRIM(MDRSF)//" "//TRIM(MDHAM)//" "//TRIM(MDIMG)//" "//
     *TRIM(MDMDG)//" "//TRIM(MDFPTSED)//" "//
     *TRIM(MDFPTBED)//" "//TRIM(MDFPTMMT)//" "//TRIM(MDFCK)//" "//
     *TRIM(MDCHK)//" "//TRIM(MDRST)//" "//TRIM(MDPOI) 

      IF (RESTAR.EQ.'cold') CALL SYSTEM(TRIM(MAKEDIR))
	FN="./"//TRIM(IN_DIRE)//TRIM(XG)//TRIM(CASENAME)//".grd"
	OPEN (IUGRD,FILE=FN,STATUS='OLD')
      READ (IUGRD,*) 
      READ (IUGRD,*)
	READ (IUGRD,*) KB
      KBM1 = KB-1
      KBM2 = KB-2
      ALLOCATE (DZR(KB))
      ALLOCATE (Z(KB))
      ALLOCATE (ZZ(KB))
      ALLOCATE (DZ(KB))
      ALLOCATE (DZZ(KB))
      DO 10 K = 1, KB
          READ (IUGRD,*) Z(K)
10    CONTINUE
      DO 20 K = 1, KBM1
          DZ(K) = Z(K) - Z(K+1)
          DZR(K) = 1. / DZ(K)
          ZZ(K) = .5 * (Z(K)+Z(K+1))
20    CONTINUE
      DO 30 K = 1, KBM2
          DZZ(K) = ZZ(K) - ZZ(K+1)
30    CONTINUE
      DZZ(KBM1) = 0.0
      DZ(KB) = 0.0 
**********DEFINE THE METRICS OF THE COORDINATE TRANSFORMATION**********
      READ (IUGRD,*) 
      READ (IUGRD,*) N_CTRD_AG,N_CTRD_EVG,N_CTRD_IVG
	N_CTRD_VG=N_CTRD_EVG+N_CTRD_IVG
      N_CTRD=N_CTRD_AG+N_CTRD_VG
      N_CTRDP1=N_CTRD+1
	FN="./"//TRIM(IN_DIRE)//TRIM(XG)//TRIM(CASENAME)//"_dep.dat"
	OPEN (IUDEP,FILE=FN,STATUS='OLD')
      FN="./"//TRIM(IN_DIRE)//TRIM(XG)//TRIM(CASENAME)//"_dchg.dat" 
      INQUIRE (FILE=FN,EXIST=LOG_DCHG)
      IF (LOG_DCHG) THEN
          OPEN (IUDCHG,FILE=FN,STATUS='OLD')
      ENDIF
      FN="./"//TRIM(IN_DIRE)//TRIM(XG)//TRIM(CASENAME)//"_cbcadj.dat" 
      INQUIRE (FILE=FN,EXIST=LOG_CBCADJ)
      IF (LOG_CBCADJ) THEN
	    OPEN (IUCBCADJ,FILE=FN,STATUS='OLD')
      ENDIF
	FN="./"//TRIM(IN_DIRE)//TRIM(XG)//TRIM(CASENAME)//"_tsr.dat"
	OPEN (IUTSR,FILE=FN,STATUS='OLD')
	FN="./"//TRIM(IN_DIRE)//TRIM(XG)//TRIM(CASENAME)//"_sec.dat"
	OPEN (IUSEC,FILE=FN,STATUS='OLD')
     
	FN="./"//TRIM(IN_DIRE)//TRIM(XG)//TRIM(CASENAME)//"_el_obc.dat"
	OPEN (IUBCS,FILE=FN,STATUS='OLD')
#ifdef TIDE_FLATHER
      FN="./"//TRIM(IN_DIRE)//TRIM(XG)//TRIM(CASENAME)//"_flux_obc.dat"
	OPEN (IUFOBC,FILE=FN,STATUS='OLD')
#endif
#ifdef HEATFLUX_BULK
      FN="./"//TRIM(IN_DIRE)//"/"//TRIM(CASENAME)//"_atp.dat"  
      OPEN (IUATP,FILE=FN,FORM='UNFORMATTED',STATUS='OLD')
      FN="./"//TRIM(IN_DIRE)//"/"//TRIM(CASENAME)//"_rhm.dat"  
      OPEN (IURHM,FILE=FN,FORM='UNFORMATTED',STATUS='OLD')
      FN="./"//TRIM(IN_DIRE)//"/"//TRIM(CASENAME)//"_cld.dat"  
      OPEN (IUCLD,FILE=FN,FORM='UNFORMATTED',STATUS='OLD')
#endif
#if defined HEATFLUX_BULK  ||  defined AIRPRESSURE 
      FN="./"//TRIM(IN_DIRE)//"/"//TRIM(CASENAME)//"_apr.dat"  
      OPEN (IUAPR,FILE=FN,FORM='UNFORMATTED',STATUS='OLD')
#endif
 
#ifdef ELB
	FN="./"//TRIM(IN_DIRE)//TRIM(XG)//TRIM(CASENAME)//"_rel_obc.dat"  !RESIDUAL ELEVATION
	OPEN (IUREBC,FILE=FN,STATUS='OLD')
#endif 
      FN="./"//TRIM(IN_DIRE)//TRIM(XG)//TRIM(CASENAME)//"_fluxbond.dat"
	OPEN (IUFBC,FILE=FN,STATUS='OLD')
	FN="./"//TRIM(IN_DIRE)//TRIM(XG)//TRIM(CASENAME)//"_wds.dat"
#ifdef WDTYP_UNIFORM
	  OPEN (IUWDS,FILE=FN,STATUS='OLD')
#elif defined WDTYP_FIELD
	  OPEN (IUWDS,FILE=FN,FORM='UNFORMATTED',STATUS='OLD')
#else
      PRINT*, 'Wind data type undefined!'
	PAUSE
	STOP
#endif
#if defined MODULE_SAL ||  defined MODULE_TMP
	  FN="./"//TRIM(IN_DIRE)//TRIM(XG)//TRIM(CASENAME)//"_its.dat"
	  OPEN (IUITS,FILE=FN,STATUS='OLD')
	  FN="./"//TRIM(IN_DIRE)//TRIM(XG)//TRIM(CASENAME)//"_ts_obc.dat"  
	  OPEN (IUSBC,FILE=FN,STATUS='OLD')
#endif

#ifdef MODULE_SED  
      FN="./"//TRIM(IN_DIRE)//TRIM(XG)//TRIM(CASENAME)//"_ised.dat"  
      INQUIRE(FILE=FN,EXIST=LOG_ISED)
      IF (LOG_ISED)THEN
	  OPEN (IUISED,FILE=FN,STATUS='OLD') 
      ENDIF
	FN="./"//TRIM(IN_DIRE)//TRIM(XG)//TRIM(CASENAME)//"_sed_obc.dat"  
	OPEN (IUSDBC,FILE=FN,STATUS='OLD')
	
	FN="./"//TRIM(IN_DIRE)//TRIM(XG)//TRIM(CASENAME)//"_sed_flux.dat"  
	OPEN (IUSDFLX,FILE=FN,STATUS='OLD')
#ifdef TAUTYP_USERDEF
      FN="./"//TRIM(IN_DIRE)//TRIM(XG)//TRIM(CASENAME)//"_d50_tau.dat"  
      OPEN (IUD50,FILE=FN,STATUS='OLD')
#  endif
#endif
#ifdef MODULE_WAVE
      FN="./"//TRIM(IN_DIRE)//TRIM(XG)//TRIM(CASENAME)//"_wave.dat"  
      OPEN (IUWAVE,FILE=FN,FORM='UNFORMATTED',CONVERT='BIG_ENDIAN',
     *    STATUS='OLD')	
#endif
#ifdef MODULE_MATERIAL
      FN="./"//TRIM(IN_DIRE)//TRIM(XG)//TRIM(CASENAME)//"_mat.dat"  
      OPEN (IUMAT,FILE=FN,STATUS='OLD')
#endif
! V2410:
#ifdef MODULE_LAG
      FN="./"//TRIM(IN_DIRE)//TRIM(XG)//TRIM(CASENAME)//"_lag.dat"  
      OPEN (IULAG,FILE=FN,STATUS='OLD')
#endif
      OPEN (IUPOI,FILE='./'//TRIM(OUT_DIRE)//TRIM(XG)//
     *'POI_ERROR/POI_ERROR.dat',position='append')
! OPEN FN TSR ! rem V2410
	NIT=400
	READ (IUTSR,*) 
	N_TSR=0
	DO WHILE (.TRUE.)
	    READ (IUTSR,*,END=4999)
	    N_TSR=N_TSR+1
	ENDDO
4999  REWIND(IUTSR)
	READ(IUTSR,*)
	ALLOCATE(FN_TSR(N_TSR))
	IF (N_TSR.NE.0) THEN
	    ALLOCATE (X_TSR(N_TSR))
	    ALLOCATE (Y_TSR(N_TSR))
          ALLOCATE (NN_TSR(N_TSR))
          ALLOCATE (WTC_TSR(NMAX,N_TSR))   ;WTC_TSR = 0.0
          ALLOCATE (WTU_TSR(NMAX,N_TSR))   ;WTU_TSR = 0.0
          ALLOCATE (WTV_TSR(NMAX,N_TSR))   ;WTV_TSR = 0.0
          ALLOCATE (NAGC_TSR(NMAX,N_TSR))  ;NAGC_TSR = 1
          ALLOCATE (NAGU_TSR(NMAX,N_TSR))  ;NAGU_TSR = 1
          ALLOCATE (NAGV_TSR(NMAX,N_TSR))  ;NAGV_TSR = 1
      ENDIF
	IF (N_TSR.NE.0) THEN
      DO NUM=1,N_TSR
	    READ(IUTSR,*) X_TSR(NUM),Y_TSR(NUM),FN
	    FN_TSR(NUM)=FN
	    FN1=FN
	    NIT=NIT+1
	    FN="./"//TRIM(OUT_DIRE)//TRIM(XG)//"TIMESERIES"//TRIM(XG)//
     *    TRIM(CASENAME)//"_"//TRIM(FN)//".OUT"
	    INQUIRE (FILE=FN,EXIST=ALIVE)
	    OPEN (NIT,FILE=FN,POSITION='APPEND')
	    IF (RESTAR.EQ.'hot') THEN
          IF (ALIVE) THEN
	        BACKSPACE(NIT)
	        READ (NIT,*) THOUR_HOT
          ELSE 
	        THOUR_HOT=0
	        WRITE (NIT,*) 'STATION:  ',TRIM(FN1)
	        WRITE (NIT,*) 'LOCATION:  ',X_TSR(NUM),Y_TSR(NUM)
          ENDIF
	    ENDIF
	    IF (RESTAR.EQ.'cold') THEN
	        WRITE (NIT,*) 'STATION:  ',TRIM(FN1)
	        WRITE (NIT,*) 'LOCATION:  ',X_TSR(NUM),Y_TSR(NUM)
	    ENDIF
      ENDDO
      ENDIF
! OPEN FN SEC !V2410
      NIT=400+N_TSR+10
      READ (IUSEC,*)
      READ (IUSEC,*) COMMENT,N_SEC_N
      K=1
      IF (N_SEC_N.NE.0) THEN
          ALLOCATE (FLUXSUM     (N_SEC_N))
          ALLOCATE (SFLUXSUM    (N_SEC_N))
          ALLOCATE (SEDFLUXSUM  (N_SEC_N))
          ALLOCATE (FLUX_ACM    (N_SEC_N))
          ALLOCATE (SFLUX_ACM   (N_SEC_N))
          ALLOCATE (SEDFLUX_ACM (N_SEC_N)) 
          ALLOCATE (TOLNUM_SEC  (N_SEC_N))
          FLUX_ACM=0.
          SFLUX_ACM=0.
          SEDFLUX_ACM=0.
      ENDIF
      
      DO N=1,N_SEC_N
          READ (IUSEC,*) COMMENT,TOLNUM_SEC(N)
          IF (TOLNUM_SEC(N).GT.K) K=TOLNUM_SEC(N)
          DO I=1,TOLNUM_SEC(N)
              READ (IUSEC,*)
          ENDDO
      ENDDO
      REWIND(IUSEC)
      IF (N_SEC_N.NE.0) THEN
          ALLOCATE (NUM_SEC(K,N_SEC_N))
      ENDIF
      READ (IUSEC,*)
      READ (IUSEC,*)
      DO N=1,N_SEC_N
          READ(IUSEC,*) FN,TOLNUM_SEC(N)
          DO I=1,TOLNUM_SEC(N)
              READ (IUSEC,*) NUM_SEC(I,N)
          ENDDO
          NIT=NIT+1
	    FN="./"//TRIM(OUT_DIRE)//TRIM(XG)//"SECFLUX"//TRIM(XG)//
     *    TRIM(CASENAME)//"_"//TRIM(FN)//".OUT"
  	    INQUIRE(FILE=FN,EXIST=ALIVE)
	    OPEN (NIT,FILE=FN,POSITION='APPEND')
          IF (RESTAR.EQ.'hot') THEN
          IF (ALIVE) THEN
          ELSE
	      WRITE (NIT,7000) 'SECTION:',TRIM(FN1)
            WRITE (NIT,5998) 'SECTION TYPE:  COORDINATE'
            WRITE (NIT,7002) 'TIME(HOUR)','FLUX(M3/S)',
     *      'ACCUMULATED_FLUX(M3)','SALT_FLUX(KG/S)',
     *      'ACCUMULATED_SALTFLUX(KG)','SED_FLUX(KG/S)',
     *      'ACCUMULATED_SEDFLUX(KG)'
          ENDIF
          ENDIF
	    IF (RESTAR.EQ.'cold') THEN
	    WRITE (NIT,7000) 'SECTION:',TRIM(FN1)
          WRITE (NIT,5998) 'SECTION TYPE:  COORDINATE'      
          WRITE (NIT,7002) 'TIME(HOUR)','FLUX(M3/S)',
     *    'ACCUMULATED_FLUX(M3)','SALT_FLUX(KG/S)',
     *    'ACCUMULATED_SALTFLUX(KG)','SED_FLUX(KG/S)',
     *    'ACCUMULATED_SEDFLUX(KG)'
          ENDIF
      ENDDO 
7000  FORMAT(A8,A10)
5998  FORMAT (A25)
7001  FORMAT(A18,2F10.0,A4,2F10.0,A13,F5.0,A2)
7002  FORMAT(A14,A16,A26,A16,A26,A16,A26) !LXY  
! --- V2410:      
!======================================================================
! LAGRANGE TRACKING      
#ifdef MODULE_LAG      
      NIT=7488+N_TSR+N_SEC_XY+20
      EARLIEST_PARTICLE=9999999.
      READ (IULAG,*)
	READ (IULAG,*) N_LAG
      IF (N_LAG.EQ.0) THEN
	    PRINT*, 'NO PARTICLE FOR TRACKING!!! PLEASE CHECK.'
	    PAUSE
	    STOP
      ELSE
          ALLOCATE (PIX_LAG(N_LAG))
          ALLOCATE (PIY_LAG(N_LAG))
          ALLOCATE (PZ_LAG(N_LAG))
          ALLOCATE (PI_LAG(N_LAG))
          ALLOCATE (PB_LAG(N_LAG))
          ALLOCATE (PS_LAG(N_LAG))
          ALLOCATE (PF_LAG(N_LAG))
          ALLOCATE (PT_LAG(N_LAG))
          ALLOCATE (PENDX_LAG(N_LAG))
          ALLOCATE (PENDY_LAG(N_LAG))
      ENDIF      
      DO N=1,N_LAG
          READ (IULAG,*) NNN,X0,Y0,Z0,TT,F0
          IF (N.NE.NNN) THEN
              PRINT*,'IN CASENAME_LAG.DAT N NOT EQ NNN'
              PAUSE
              STOP
          ENDIF
          PB_LAG(N)=0
          PENDX_LAG(N)=X0
	    PENDY_LAG(N)=Y0
	    PZ_LAG(N)=Z0
	    PT_LAG(N)=TT
          PF_LAG(N)=F0
	    WRITE (FN,500) N
500	    FORMAT ('LAG',I3.3)
	    NIT=NIT+1
          FN="./"//TRIM(OUT_DIRE)//TRIM(XG)//"LAG_TRACKING/"//
     *	TRIM(CASENAME)//"_"//TRIM(FN)//".OUT"
          INQUIRE(FILE=FN,EXIST=ALIVE)
          
	    OPEN (NIT,FILE=FN,POSITION='APPEND')
          IF (RESTAR.EQ.'hot') THEN
          IF (ALIVE) THEN
          ELSE 
	        WRITE (NIT,501) 'LAGRANGIAN PARTICLE TRACKING RESULT:'
              WRITE (NIT,502) 'INITIAL POSITION:  ',X0,Y0,Z0,TT
              WRITE (NIT,503) 'TIME(HOUR)  X-COORDINATE  Y-COORDINATE
     *  Z-COORDINATE'
          ENDIF
	    ENDIF
	    IF (RESTAR.EQ.'cold') THEN
	        WRITE (NIT,501) 'LAGRANGIAN PARTICLE TRACKING RESULT:'
              WRITE (NIT,502) 'INITIAL POSITION:  ',X0,Y0,Z0,TT
              WRITE (NIT,503) 'TIME(HOUR)  X-COORDINATE  Y-COORDINATE
     *  Z-COORDINATE'
	    ENDIF
501	    FORMAT(A36)
502	    FORMAT(A19,2F10.0,2F8.2)
503       FORMAT (A108)
          IF (TT.LT.EARLIEST_PARTICLE) EARLIEST_PARTICLE=TT

      ENDDO

#endif
! END: LAGRANGE TRACKING    
!======================================================================   
!======================================================================
! SLUICE   
#ifdef MODULE_SLUICE 
      FN="./"//trim(in_dire)//trim(xg)//trim(casename)
     *//"_sluice.dat"
      OPEN (IUSLUICE,FILE=FN,STATUS='OLD')
      READ (IUSLUICE,*)
      READ (IUSLUICE,*) N_SLUICE
      IF (N_SLUICE.NE.0) THEN
          ALLOCATE (SLUICE_NUM(N_SLUICE))
          ALLOCATE (LOG_SLUICE_ON(N_SLUICE))
	    ALLOCATE (SLUICE_TYP(N_SLUICE))
	    ALLOCATE (SLUICE_GRD(N_SLUICE,500))
          ALLOCATE (SLUICE_FSM(N_SLUICE,500))
      ENDIF
      LOG_SLUICE_ON=.FALSE.
#endif
! END: SLUICE   
!======================================================================   
! --- V2410.      

************** RESIDUAL CURRENT AND RESIDUAL FLUX FILES ***************
	NIT1=2000
	NIT2=2200
	NIT3=2400
	NIT4=2600
	NIT5=2800
	IF (JHM.NE.0) THEN   !BEGIN
	  DO J=1,JHM
	    NIT=NIT1+J
	    WRITE (FNN,'(I2.2)') J
	    FN="./"//TRIM(OUT_DIRE)//TRIM(XG)//"RESI_EULER"//TRIM(XG)//
     *"RESI_EULER_CURRENT_"//TRIM(FNN)//".OUT"
	    OPEN (NIT,FILE=FN,STATUS='UNKNOWN')
	  ENDDO
#ifdef MODULE_SAL
	  DO J=1,JHM
	    NIT=NIT2+J
	    WRITE (FNN,'(I2.2)') J
	    FN="./"//TRIM(OUT_DIRE)//TRIM(XG)//"RESI_FLUX"//TRIM(XG)//
     *"RESI_SALT_FLUX_3D_"//TRIM(FNN)//".OUT"
	    OPEN (NIT,FILE=FN,STATUS='UNKNOWN')
	  ENDDO
#endif
	  DO J=1,JHM
	    NIT=NIT3+J
	    WRITE (FNN,'(I2.2)') J
	FN="./"//TRIM(OUT_DIRE)//TRIM(XG)//"RESI_FLUX"//TRIM(XG)
     *//"RESI_UNITFLUX_WATERMASS_"//TRIM(FNN)//".OUT"
	    OPEN (NIT,FILE=FN,STATUS='UNKNOWN')
	  ENDDO
#ifdef MODULE_SAL
	  DO J=1,JHM
	    NIT=NIT4+J
	    WRITE (FNN,'(I2.2)') J
	    FN="./"//TRIM(OUT_DIRE)//TRIM(XG)//"RESI_FLUX"//TRIM(XG)//
     *"RESI_UNITFLUX_SALT_"//TRIM(FNN)//".OUT"
	    OPEN (NIT,FILE=FN,STATUS='UNKNOWN')
        ENDDO
#endif
      ENDIF
	DO J=1,JHM
	  NIT=NIT5+J
	  WRITE (FNN,'(I2.2)') J
	  FN="./"//TRIM(OUT_DIRE)//TRIM(XG)//"RESI_FLUX"//TRIM(XG)//
     *"RESI_FLUX_3D_"//TRIM(FNN)//".OUT"
	  OPEN (NIT,FILE=FN,STATUS='UNKNOWN')
      ENDDO
      FN = "./"//TRIM(OUT_DIRE)//TRIM(XG)//"TIMESTEP"//TRIM(XG)
     *//"TIMESTEP.DAT"
      OPEN (IUTSP,FILE=FN,POSITION='APPEND')
#ifdef HEATFLUX_BULK
      FN="./MODGEN"//TRIM(XG)//TRIM(CASENAME)//".T86"
	OPEN (IUT86,FILE=FN,FORM='UNFORMATTED')
      FN="./MODGEN"//TRIM(XG)//TRIM(CASENAME)//".T87"
	OPEN (IUT87,FILE=FN,FORM='UNFORMATTED')
      FN="./MODGEN"//TRIM(XG)//TRIM(CASENAME)//".T88"
	OPEN (IUT88,FILE=FN,FORM='UNFORMATTED')
#endif
#if defined HEATFLUX_BULK  ||  defined AIRPRESSURE 
      FN="./MODGEN"//TRIM(XG)//TRIM(CASENAME)//".T89"
	OPEN (IUT89,FILE=FN,FORM='UNFORMATTED')
#endif
      
	FN="./MODGEN"//TRIM(XG)//TRIM(CASENAME)//".T90" 
	OPEN (IUT90,FILE=FN)
	FN="./MODGEN"//TRIM(XG)//TRIM(CASENAME)//".T91"
	OPEN (IUT91,FILE=FN)
	FN="./MODGEN"//TRIM(XG)//TRIM(CASENAME)//".T92"
	OPEN (IUT92,FILE=FN,FORM='UNFORMATTED')	
	FN="./MODGEN"//TRIM(XG)//TRIM(CASENAME)//".T93"
	OPEN (IUT93,FILE=FN,FORM='UNFORMATTED')
#if defined MODULE_SAL ||  defined MODULE_TMP
	FN="./MODGEN"//TRIM(XG)//TRIM(CASENAME)//".T94"
	OPEN (IUT94,FILE=FN)
     
#endif
#ifdef ELB
	FN="./MODGEN"//TRIM(XG)//TRIM(CASENAME)//".T95"
	OPEN (IUT95,FILE=FN)
#endif
#if defined TIDE_FLUX  || defined TIDE_FLATHER
      FN="./MODGEN"//TRIM(XG)//TRIM(CASENAME)//".T96"
	OPEN (IUT96,FILE=FN)
#endif
#ifdef MODULE_SED
	  FN="./MODGEN"//TRIM(XG)//TRIM(CASENAME)//".T21"
	  OPEN (IUT21,FILE=FN)
	  FN="./MODGEN"//TRIM(XG)//TRIM(CASENAME)//".T20"
	  OPEN (IUT20,FILE=FN)	
#endif
#ifdef MODULE_MATERIAL
	FN="./MODGEN"//TRIM(XG)//TRIM(CASENAME)//".T98"
	  OPEN (IUT98,FILE=FN)
#endif		
	FN="./MODGEN"//TRIM(XG)//TRIM(CASENAME)//"_GRIDDEPTH" 
      OPEN (IUGDH,FILE=FN,FORM='UNFORMATTED')
	FN="./MODGEN"//TRIM(XG)//TRIM(CASENAME)//".PRT"
	OPEN (IUPRT,FILE=FN)
	FN="./MODGEN"//TRIM(XG)//"UV_CHECK.DAT"
	OPEN (IUUVK,FILE=FN) 
***********************************************************************
*                                                                     *
*                               CASE 2                                *
*                                                                     *
***********************************************************************
      CASE(2)    
#ifdef WEIR
      FN="./"//TRIM(IN_DIRE)//TRIM(XG)//TRIM(CASENAME)//"_weir.dat"
      OPEN (IUWEIR,FILE=FN,STATUS='OLD')
      READ (IUWEIR,*)
      READ (IUWEIR,*) ND
      IF (ND.NE.0) THEN
        ALLOCATE (WEIR_N(ND))  
        ALLOCATE (WEIR_TYP(ND))
        ALLOCATE (WEIR_I(ND,500))
        ALLOCATE (WEIR_H(ND,500))     
        ALLOCATE (DELTA0(ND,500))
        ALLOCATE (DELTA1(ND,500))
        ALLOCATE (OCSL(ND,500))
        ALLOCATE (KB_WEIR_0(ND,500))
        ALLOCATE (KB_WEIR_1(ND,500))
        ALLOCATE (KB_WEIR_5(ND,500))
      ENDIF
      DO I=1,ND
        READ (IUWEIR,*) WEIR_N(I)  
        DO J=1,WEIR_N(I)
          READ (IUWEIR,*) WEIR_I(I,J),WEIR_H(I,J) !LTY      
        ENDDO
        IF (WEIR_I(I,1).EQ.NJM1(WEIR_I(I,2))
     *   .OR.WEIR_I(I,1).EQ.NJP1(WEIR_I(I,2))) THEN
          WEIR_TYP(I)='J_DIRECT'
        ELSEIF (WEIR_I(I,1).EQ.NIM1(WEIR_I(I,2))
     *  .OR.WEIR_I(I,1).EQ.NIP1(WEIR_I(I,2))) THEN
          WEIR_TYP(I)='I_DIRECT'
        ELSE
          PRINT*, 'WRONG WEIR SETTING UP!'
          PAUSE
          STOP
        ENDIF
      ENDDO
#endif
      END SELECT
	END SUBROUTINE OPENFILE