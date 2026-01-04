#include "DEFS.h"
	
	SUBROUTINE READ_NML     
	USE MOD_GLOBAL

	CHARACTER*200 COMMENT
	INTEGER OMP_GET_NUM_THREADS,OMP_GET_MAX_THREADS
	INTEGER OMP_GET_NUM_PROCS
	
      FN = "./namelist.txt"
	OPEN (IUNML,FILE=FN,STATUS='OLD')
	FN="INFO.TXT"
	OPEN (IUIFO,FILE=FN)
***********************************************************************
*                                                                     *
*                       READ THE *_RUN DATA                           *
*                                                                     *
***********************************************************************

**************************READ CASE NAME*******************************
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'CASENAME')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,CASENAME
      WRITE (IUIFO,*) COMMENT,CASENAME
*************************READ OS PARAMETER*****************************
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'OSTYPE')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,OSTYPE
      WRITE (IUIFO,*) COMMENT,OSTYPE
      IF (OSTYPE.NE.'windows'.AND.OSTYPE.NE.'linux') THEN
          WRITE (*,*) 'ERROR OS PARAMETER.'
           PAUSE
          STOP
      ENDIF
      IF (TRIM(OSTYPE).EQ.'windows') THEN
          XG = '\\'
      ELSEIF (TRIM(OSTYPE).EQ.'linux') THEN
          XG = '/'
      ENDIF
************************READ PARALLEL SETTINGS*************************
c#ifdef OMP
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'NTHREADS')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,NTHREADS
      WRITE (IUIFO,*) COMMENT,NTHREADS
      IF (NTHREADS.LE.0) THEN
          PRINT*, 'INVALID NTHREADS!'
          PAUSE
          STOP
      ENDIF
      CALL OMP_SET_NUM_THREADS(NTHREADS)
      CALL ACC_SET_NUM_CORES(NTHREADS)
c#endif
************************READ COORDINATE SELCETION**********************
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'COORD')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,COORD
      WRITE (IUIFO,*) COMMENT,COORD
      IF (COORD.NE.'BL'.AND.COORD.NE.'XY') THEN
          PRINT*, 'INVALID COORDINATE.'
          PAUSE
          STOP
      ENDIF
**********************READ INPUT AND OUTPUT DIRECTORY******************
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'IN_DIRE')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,IN_DIRE
      WRITE (IUIFO,*) COMMENT,IN_DIRE
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'OUT_DIRE')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,OUT_DIRE
      WRITE (IUIFO,*) COMMENT,OUT_DIRE
*************************READ THE TIME SETTINGS************************
	REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'DTITYPE')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,DTITYPE
      WRITE (IUIFO,*) COMMENT,DTITYPE
	IF (TRIM(DTITYPE).NE.'constant'.AND.TRIM(DTITYPE).NE.'variable') THEN
	    WRITE (*,*) 'WRONG DTI TYPE!!!!'
	    PAUSE
	    STOP
	ENDIF
	IF (DTITYPE.EQ.'constant') THEN
          REWIND(IUNML)
          DO WHILE (TRIM(COMMENT).NE.'DTI')
          READ (IUNML,*,END=161) COMMENT
          ENDDO
          BACKSPACE (IUNML)
          READ (IUNML,*) COMMENT,DTI
          WRITE (IUIFO,*) COMMENT,DTI
	ELSEIF (DTITYPE.EQ.'variable') THEN
          REWIND(IUNML)
          DO WHILE (TRIM(COMMENT).NE.'DTI')
          READ (IUNML,*,END=161) COMMENT
          ENDDO
          BACKSPACE (IUNML)
          READ (IUNML,*,END=161) COMMENT,DTIMAX1
          WRITE (IUIFO,*) COMMENT,DTIMAX1
	ENDIF
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'FACTOR')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,FACTOR
      WRITE (IUIFO,*) COMMENT,FACTOR
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'IRAMP')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,IRAMP
      WRITE (IUIFO,*) COMMENT,IRAMP
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'T_END')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,T_END
      WRITE (IUIFO,*) COMMENT,T_END
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'IYEAR')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,IYEAR
      WRITE (IUIFO,*) COMMENT,IYEAR
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'IMONTH')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,IMONTH
      WRITE (IUIFO,*) COMMENT,IMONTH
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'IDAY0')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,IDAY0
      WRITE (IUIFO,*) COMMENT,IDAY0
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'TIDE_LAG')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,TIDE_LAG
      WRITE (IUIFO,*) COMMENT,TIDE_LAG
*************************READ THE HOTSTART DATA************************
	REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'RESTAR')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,RESTAR
      WRITE (IUIFO,*) COMMENT,RESTAR
      IF (RESTAR.NE.'cold'.AND.RESTAR.NE.'hot') THEN
          PRINT*, 'INVALID RESTART PARAMETER !'
          PAUSE
          STOP
      ENDIF
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'N_RST')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,N_RST
      WRITE (IUIFO,*) COMMENT,N_RST
	REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'RST_SPEC')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,RST_SPEC
      WRITE (IUIFO,*) COMMENT,RST_SPEC
#ifdef RSTARC
	REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'RST_ARCH')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,RST_ARCH
      WRITE (IUIFO,*) COMMENT,RST_ARCH
#endif
*************READ ELEVATION INPUT DATA TYPE AT OPEN BOUNDARY***********
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'OPTEBC')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,OPTEBC !"harmconst" OR "data", ELEVATION OPEN BOUNDARY CONDITION INPUT TYPE.
      WRITE (IUIFO,*) COMMENT,OPTEBC
      IF (OPTEBC.NE.'harmconst'.AND.OPTEBC.NE.'data') THEN
          PRINT*, 'INVALID OPTEBC.'
          PAUSE
          STOP
      ENDIF
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'EL_RISEUP')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,EL_RISEUP
      WRITE (IUIFO,*) COMMENT,EL_RISEUP
      EL_RISEUP=EL_RISEUP/100.
#ifdef MODULE_SAL
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'S_BEG')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,S_BEG
      WRITE (IUIFO,*) COMMENT,S_BEG
#endif
#ifdef BPG
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'BPG_BEG')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,BPG_BEG
      WRITE (IUIFO,*) COMMENT,BPG_BEG
#endif
#ifdef MODULE_TMP
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'T_BEG')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,T_BEG
      WRITE (IUIFO,*) COMMENT,T_BEG
#endif
#ifdef MODULE_MATERIAL
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'M_BEG')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,M_BEG
      WRITE (IUIFO,*) COMMENT,M_BEG
#endif
#ifdef MODULE_SED
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'SED_BEG')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,SED_BEG
      WRITE (IUIFO,*) COMMENT,SED_BEG
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'F_SD50')
        READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,F_SD50
      WRITE (IUIFO,*) COMMENT,F_SD50
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'F_ALFA')
        READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,F_ALFA
      WRITE (IUIFO,*) COMMENT,F_ALFA	
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'F_MCSXS')
        READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,F_MCSXS
      WRITE (IUIFO,*) COMMENT,F_MCSXS
#  ifdef BED
          REWIND(IUNML)
          DO WHILE (TRIM(COMMENT).NE.'BED_BEG')
              READ (IUNML,*,END=161) COMMENT
          ENDDO
          BACKSPACE (IUNML)
          READ (IUNML,*) COMMENT,BED_BEG
          WRITE (IUIFO,*) COMMENT,BED_BEG
#  endif
#endif      
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'BFRIC')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,BFRIC
      WRITE (IUIFO,*) COMMENT,BFRIC
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'Z0B')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,Z0B
      WRITE (IUIFO,*) COMMENT,Z0B
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'HORCON')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,HORCON
      WRITE (IUIFO,*) COMMENT,HORCON
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'HPRNU')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,HPRNU
      WRITE (IUIFO,*) COMMENT,HPRNU
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'ADVECT')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,ADVECT
      WRITE (IUIFO,*) COMMENT,ADVECT
      IF (ADVECT.NE.'non-linear') THEN
          WRITE (*,*)' INVALID ADVECT!'
          PAUSE
          STOP
      ENDIF
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'UMOL')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,UMOL
      WRITE (IUIFO,*) COMMENT,UMOL
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'VPRNU')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,VPRNU
      WRITE (IUIFO,*) COMMENT,VPRNU
***********************READ STANDARD LEVEL DEPTH DATA******************
	REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'IKSL')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,KSL
      WRITE (IUIFO,*) COMMENT,KSL
      ALLOCATE (DPTHSL(KSL))
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'DPTHSL')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,(DPTHSL(K),K = 1,KSL)
**************************READ DRY/WET PARAMETERS**********************
	REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'DMIN')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,DMIN
      WRITE (IUIFO,*) COMMENT,DMIN
**********************READ TIMESERIES OUTPUT PARAMETER*****************
	REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'N_OPT')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,N_OPT
      WRITE (IUIFO,*) COMMENT,N_OPT
#ifdef FPT
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'N_FPT')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,N_FPT
      WRITE (IUIFO,*) COMMENT,N_FPT
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'FPTSTAR')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,FPTSTAR
      WRITE (IUIFO,*) COMMENT,FPTSTAR
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'FPTEND')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,FPTEND
      WRITE (IUIFO,*) COMMENT,FPTEND
#endif
******************READ RESIDUAL AND RESIDUAL FLUX PARAMETER************
      JHM=0
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'JHM')
          READ (IUNML,*,END=108) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,JHM
108   CONTINUE
      WRITE (IUIFO,*) COMMENT,JHM
      IF (JHM.NE.0) THEN
      ALLOCATE (HIST(JHM,2))
      ALLOCATE (AVGE(JHM))
      ALLOCATE (DEI(JHM))
      ALLOCATE (LOG_ARC(JHM))
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'HIST&AVGE')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE (IUNML)
      READ (IUNML,*) COMMENT,(HIST(I,2),AVGE(I),I=1,JHM)
      DO I=1,JHM
          LOG_ARC(I)=.TRUE.
          DEI(I)=1./AVGE(I)
      ENDDO
      DO I = 1,JHM
          HIST(I,1) = HIST(I,2) - AVGE(I)
          IF (HIST(I,1)/24..GE.T_END) THEN
            PRINT*, 'TIME FOR RESIDUAL COMPUTATION IS BEYOND THE 
     * MODEL TIME SETUP. PLEASE ENLARGE THE VALUE OF T_END'
              PAUSE
          ENDIF
      ENDDO
      ENDIF
      
#ifdef MOD_ONLINE_VIEW
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'CHC_VIEW')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE(IUNML)
      READ (IUNML,*) COMMENT,CHC_VIEW
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'LAYER')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE(IUNML)
      READ (IUNML,*) COMMENT,LAYER
#endif

c --- V2410:
#ifdef MODULE_LAG
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'NDT')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE(IUNML)
      READ (IUNML,*) COMMENT,NDT
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'I_LAG')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE(IUNML)
      READ (IUNML,*) COMMENT,I_LAG
#endif
      REWIND(IUNML)
      DO WHILE (TRIM(COMMENT).NE.'LOG_TH')
          READ (IUNML,*,END=161) COMMENT
      ENDDO
      BACKSPACE(IUNML)
      READ (IUNML,*) COMMENT,LOG_TH
      IF (LOG_TH) THEN
          READ (IUNML,*) COMMENT,YEAR_TH,MONTH_TH,DAY_TH
          READ (IUNML,*) COMMENT,HOUR_TH,MINUTE_TH,SECOND_TH
      ENDIF
c --- V2410.

      CLOSE (IUNML)
#ifdef OMP
!======================================================================
!							  SET THREADS
!======================================================================
      PRINT*, 'OMP_GET_NUM_THREADS: ', OMP_GET_NUM_THREADS()
      PRINT*, 'OMP_GET_MAX_THREADS: ', OMP_GET_MAX_THREADS()
      PRINT*, 'OMP_GET_NUM_PROCS: ', OMP_GET_NUM_PROCS()
      PRINT*, 'NTHREADS: ', NTHREADS
!======================================================================
!							END: SET THREADS
!======================================================================
#endif

! - del rem V2410

6002  FORMAT (A36)
6003  FORMAT (A37)
6004  FORMAT (A43)
6008  FORMAT (100A)
	RETURN 
161   WRITE(*,*) 'SOME PARAMETER CANNOT BE FOUND IN THE RUN FILE.'
      PAUSE
      STOP
	END SUBROUTINE READ_NML