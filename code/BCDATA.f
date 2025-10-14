#include "DEFS.h"
	SUBROUTINE BCDATA
      USE MOD_GLOBAL
      DIMENSION COM(80)
      DIMENSION TXX(N_CTRD), TYY(N_CTRD)
#ifdef MODULE_WAVE
      REAL WTEMP1(N_CTRD),WTEMP2(N_CTRD),WTEMP3(N_CTRD),WTEMP4(N_CTRD)
#endif
***********************************************************************
*                                                                     *
*                               EL_OBC                                *
*                                                                     *
***********************************************************************
#if defined TIDE_EL || defined TIDE_FLATHER  
      READ (IUBCS,5000) (COM(I),I = 1,80)
      WRITE (IUPRT,5100) (COM(I),I = 1,80)
      READ (IUBCS,*) NUM_HARMCONST
      READ (IUBCS,*) NUMEBC
      WRITE (IUPRT,5400) NUMEBC
      IF (NUMEBC.NE.0) THEN
        ALLOCATE(NETA(NUMEBC))
        ALLOCATE(NCON(NUMEBC)) 
        ALLOCATE(NAJA(NUMEBC)) !CBR: NETA - NCON - NAJA
        ALLOCATE (AMP(NUMEBC,NUM_HARMCONST))
        ALLOCATE (PHASE(NUMEBC,NUM_HARMCONST))
        ALLOCATE (EMEAN((NUMEBC)))
        ALLOCATE (TBDRYSL(NUMEBC,KB))
        ALLOCATE (SBDRYSL(NUMEBC,KB))
        ALLOCATE (TBDRY(NUMEBC,KBM1))
        ALLOCATE (DTBDRY(NUMEBC,KBM1,2))
        ALLOCATE (SBDRY(NUMEBC,KBM1))
        ALLOCATE (DSBDRY(NUMEBC,KBM1,2))
        ALLOCATE (EBDRY(NUMEBC))
        ALLOCATE (DEBDRY(NUMEBC,2))
        ALLOCATE (ELBC(NUMEBC))
        ALLOCATE (ELBC_TIDE(NUMEBC)) !CBR
        ALLOCATE (DELBC(NUMEBC,2))
        ALLOCATE (SEDBDRYSL(NUMEBC,KB))
        ALLOCATE (SEDBDRY(NUMEBC,KBM1))
        ALLOCATE (DSEDBDRY(NUMEBC,KBM1,2))
      ENDIF
      IF (NUMEBC.NE.0) THEN
      IF (TRIM(OPTEBC).NE.'data') THEN
      DO N = 1, NUMEBC
            READ (IUBCS,*) NETA(N), NCON(N), EMEAN(N)
            READ (IUBCS,*) (AMP(N,I),I = 1,NUM_HARMCONST)
            READ (IUBCS,*) (PHASE(N,I),I = 1,NUM_HARMCONST)
            WRITE(IUPRT,5601) NETA(N), NCON(N), EMEAN(N)
            WRITE (IUPRT,'(16F10.5)') (AMP(N,I),I = 1,NUM_HARMCONST)
            WRITE (IUPRT,'(16F10.5)') (PHASE(N,I),I = 1,NUM_HARMCONST)
      ENDDO
      ELSE
          READ (IUBCS,*) (NETA(N),NCON(N),N = 1,NUMEBC)
          WRITE (IUPRT,5801) (NETA(N),NCON(N),N = 1,NUMEBC)
          DO I = 1, 100000
              READ (IUBCS,6400,ERR=60,END=60) TIME1
              READ (IUBCS,6400) (EBDRY(N),N = 1,NUMEBC)
              WRITE (IUPRT,6400) TIME1
              WRITE (IUPRT,6400) (EBDRY(N),N = 1,NUMEBC)
              WRITE (IUT90,6500) TIME1
              WRITE (IUT90,6500) (EBDRY(N),N = 1,NUMEBC)
          ENDDO
60        CONTINUE
      ENDIF
      ENDIF
      CLOSE(IUBCS)
#ifdef TIDE_FLATHER
      READ (IUFOBC,*)
      ALLOCATE (VFOBC(NUMEBC,KB))
      ALLOCATE (DVFOBC(NUMEBC,KB,2))
      DO II=1,10000000
      READ (IUFOBC,*,END=133,ERR=133) TIME1
      WRITE (IUT96,6500) TIME1
      DO N=1,NUMEBC
      READ (IUFOBC,*) (VFOBC(N,K), K=1,KBM1)
      ENDDO
      DO N=1,NUMEBC
      WRITE (IUT96,6501) (VFOBC(N,K), K=1,KBM1)
      ENDDO
      ENDDO
133   CONTINUE
      CLOSE (IUFOBC)
#endif      
#elif defined TIDE_FLUX
      READ (IUBCS,5000) (COM(I),I = 1,80)
      WRITE (IUPRT,5100) (COM(I),I = 1,80)
      READ (IUBCS,*) NUM_HARMCONST
      READ (IUBCS,*) NUMEBC
      WRITE (IUPRT,5400) NUMEBC
      IF (NUMEBC.NE.0) THEN
        ALLOCATE(NETA(NUMEBC))
        ALLOCATE(NCON(NUMEBC)) 
        ALLOCATE(NAJA(NUMEBC)) !CBR: NETA - NCON - NAJA
        ALLOCATE (AMP(NUMEBC,NUM_HARMCONST))
        ALLOCATE (PHASE(NUMEBC,NUM_HARMCONST))
        ALLOCATE (EMEAN((NUMEBC)))
        ALLOCATE (TBDRYSL(NUMEBC,KB))
        ALLOCATE (SBDRYSL(NUMEBC,KB))
        ALLOCATE (TBDRY(NUMEBC,KBM1))
        ALLOCATE (DTBDRY(NUMEBC,KBM1,2))
        ALLOCATE (SBDRY(NUMEBC,KBM1))
        ALLOCATE (DSBDRY(NUMEBC,KBM1,2))
        ALLOCATE (EBDRY(NUMEBC))
        ALLOCATE (DEBDRY(NUMEBC,2))
        ALLOCATE (ELBC(NUMEBC))
        ALLOCATE (DELBC(NUMEBC,2))
        ALLOCATE (SEDBDRYSL(NUMEBC,KB))
        ALLOCATE (SEDBDRY(NUMEBC,KBM1))
        ALLOCATE (DSEDBDRY(NUMEBC,KBM1,2))
      ENDIF
      IF (NUMEBC.NE.0) THEN
      DO N = 1, NUMEBC
            READ (IUBCS,*) NCON(N),NETA(N) 
            READ (IUBCS,*) (AMP(N,I),I = 1,NUM_HARMCONST)
            READ (IUBCS,*) (PHASE(N,I),I = 1,NUM_HARMCONST)
            WRITE(IUPRT,5601) NETA(N), NCON(N), EMEAN(N)
            WRITE (IUPRT,'(16F10.5)') (AMP(N,I),I = 1,NUM_HARMCONST)
            WRITE (IUPRT,'(16F10.5)') (PHASE(N,I),I = 1,NUM_HARMCONST)
      ENDDO
      ENDIF
      CLOSE(IUBCS)
#endif
      ALLOCATE(HA(NUM_HARMCONST))
      ALLOCATE(GA(NUM_HARMCONST))
***********************************************************************
*                                                                     *
*                             REL_OBC                                 *
*                                                                     *
***********************************************************************
#ifdef ELB
	READ (IUREBC,*)
      READ (IUREBC,*)
      WRITE (IUPRT,5400) NUMEBC
      IF (NUMEBC.NE.0) THEN
          DO I = 1, 100000
            READ (IUREBC,*,ERR=41,END=41) TIME1
            WRITE (IUPRT,*) TIME1
!            READ (IUREBC,*) (ELBC(N),N = 1,NUMEBC)
            READ (IUREBC,"(5F10.5)") (ELBC(N),N = 1,NUMEBC) !CBR: Formatted input for standard fortran    
            WRITE (IUPRT,6400) (ELBC(N),N = 1,NUMEBC)
            WRITE (IUT95,6500) TIME1
            WRITE (IUT95,6500) (ELBC(N),N = 1,NUMEBC)
          ENDDO
41        CONTINUE
            ENDIF
	CLOSE(IUREBC)
#endif
***********************************************************************
*                                                                     *
*                               TS_OBC                                *
*                                                                     *
***********************************************************************
#if defined MODULE_SAL ||  defined MODULE_TMP
      READ (IUSBC,*)
      DO I = 1, 100000
!          READ (IUSBC,6402,ERR=110,END=110) TIME1
          READ (IUSBC,*,ERR=110,END=110) TIME1
          WRITE (IUPRT,6402) TIME1
          WRITE (IUT94,6500) TIME1
          DO N = 1,NUMEBC
!              READ (IUSBC,*) (TBDRYSL(N,K),K = 1,KBM1)
!              READ (IUSBC,*) (SBDRYSL(N,K),K = 1,KBM1) 
              READ (IUSBC,"(10F5.1)") (TBDRYSL(N,K),K = 1,KBM1)
              READ (IUSBC,"(10F5.1)") (SBDRYSL(N,K),K = 1,KBM1) 
              WRITE (IUPRT,5900) (TBDRYSL(N,K),K = 1,KBM1)
              WRITE (IUPRT,5900) (SBDRYSL(N,K),K = 1,KBM1)
              WRITE (IUT94,6500) (TBDRYSL(N,K),K = 1,KBM1)   !WUHUI
              WRITE (IUT94,6500) (SBDRYSL(N,K),K = 1,KBM1)
          ENDDO
      ENDDO
110   CONTINUE
      
      CLOSE(IUSBC)
#endif
***********************************************************************
*                                                                     *
*                               FLUX                                  *
*                                                                     *
***********************************************************************
!======================================================================
!       RIVER/DAM AND ONSHORE INTAKE/OUTFALL DISCHARGE BOUNDARY
!======================================================================
      READ (IUFBC,5000) (COM(I),I = 1,80)
      WRITE (IUPRT,5100) (COM(I),I = 1,80)      
      READ (IUFBC,*) NUMQBC
      WRITE (IUPRT,5500) NUMQBC
      IF (NUMQBC.NE.0) THEN
          ALLOCATE (NQD(NUMQBC))
          ALLOCATE (NQC(NUMQBC))
          ALLOCATE (VQDIST(NUMQBC,KB))
          ALLOCATE (QDIS(NUMQBC))
          ALLOCATE (TDIS(NUMQBC))
          ALLOCATE (SDIS(NUMQBC))
          ALLOCATE (DQDIS(NUMQBC,2))
          ALLOCATE (DTDIS(NUMQBC,2))
          ALLOCATE (DSDIS(NUMQBC,2))
      ENDIF      
      IF (NUMQBC.NE.0) THEN
      DO N=1,NUMQBC          
          READ(IUFBC,*) NQD(N),NQC(N),(VQDIST(N,K),K = 1,KBM1)
          WRITE(IUPRT,6001) NQD(N),NQC(N),(VQDIST(N,K),K = 1,KBM1)
      ENDDO
      DO I = 1, 100000
          READ (IUFBC,*,ERR=150,END=150) TIME1
          READ (IUFBC,*) (QDIS(N),N = 1,NUMQBC)
          READ (IUFBC,*) (TDIS(N),N = 1,NUMQBC)
          READ (IUFBC,*) (SDIS(N),N = 1,NUMQBC)
          WRITE (IUPRT,6402) TIME1
          WRITE (IUPRT,6400) (QDIS(N),N = 1,NUMQBC)
          WRITE (IUPRT,6400) (TDIS(N),N = 1,NUMQBC)
          WRITE (IUPRT,6400) (SDIS(N),N = 1,NUMQBC)
          WRITE (IUT91,6500) TIME1
          WRITE (IUT91,6500) (QDIS(N),N = 1,NUMQBC)
          WRITE (IUT91,6500) (TDIS(N),N = 1,NUMQBC)
          WRITE (IUT91,6500) (SDIS(N),N = 1,NUMQBC)
      ENDDO
150   CONTINUE
      ENDIF
      CLOSE(IUFBC)
***********************************************************************
*                                                                     *
*                            METEOROLOGICAL                           *
*                                                                     *
***********************************************************************
!======================================================================
!    METEOROLOGICAL BOUNDARY CONDITION
!    PRECIPITATION (M/YEAR), EVAPORATION (M/YEAR), WIND FROM DIRECTION
!======================================================================
#ifdef WDTYP_UNIFORM
      READ (IUWDS,*)  
      DO I = 1, 100000
          READ (IUWDS,6402,ERR=190,END=190) TIME1
          WRITE (IUPRT,6402) TIME1
          READ (IUWDS,6400) QPREC, QEVAP, WDS, WDD, HFLUX
          WRITE (IUPRT,6400) QPREC, QEVAP, WDS, WDD, HFLUX
          QPREC = QPREC / (86400.*365.)
          QEVAP = QEVAP / (86400.*365.)
          WDD = 180. + WDD
          WDD = AMOD(WDD,360.)
          VWIND = WDS * COS(6.28319*WDD/360.)
          UWIND = WDS * SIN(6.28319*WDD/360.)
          CD = 1.2E-3
          IF (WDS.GE.11.) CD = (0.49+0.065*WDS) * 1.E-3
          IF (WDS.GE.25.) CD = (0.49+0.065*25.) * 1.E-3
          TX = 1.2 * CD * UWIND * WDS
          TY = 1.2 * CD * VWIND * WDS
          WRITE (IUT93) TIME1
          WRITE (IUT93) QPREC, QEVAP, TX, TY, HFLUX
      ENDDO
190   CONTINUE
#elif defined WDTYP_FIELD    
      DO I=1,1000000
          READ (IUWDS,ERR=191,END=191) TIME1 !BIN
          READ (IUWDS) WINDU,WINDV
          DO II = 1,N_CTRD
              CD = 1.2E-3   
	        WDS = SQRT(WINDU(II)**2+WINDV(II)**2)
              IF (WDS.GE.11.) CD = (0.49+0.065*WDS) * 1.E-3   
              IF (WDS.GE.25.) CD = (0.49+0.065*25.) * 1.E-3
              TXX(II) = 1.2 * CD * WINDU(II) * WDS
              TYY(II) = 1.2 * CD * WINDV(II) * WDS
	    ENDDO
          WRITE (IUT93) TIME1
          WRITE (IUT93) TXX,TYY
      ENDDO
191   CONTINUE
#else
      PRINT*, 'Wind data type undefined!'
	PAUSE
	STOP
#endif
#ifdef HEATFLUX_BULK
      DO I = 1, 100000
          READ (IUATP,END=192) TIME1
          READ (IUATP) ATP
          WRITE (IUT86) TIME1
          WRITE (IUT86) ATP
      ENDDO
192   CONTINUE
      DO I = 1, 100000
          READ (IURHM,END=193) TIME1
          READ (IURHM) RHM
          WRITE (IUT87) TIME1
          WRITE (IUT87) RHM
      ENDDO
193   CONTINUE
      DO I = 1, 100000
          READ (IUCLD,END=194) TIME1
          READ (IUCLD) CLOUD
          WRITE (IUT88) TIME1
          WRITE (IUT88) CLOUD
      ENDDO
194   CONTINUE
#endif
#if defined HEATFLUX_BULK  ||  defined AIRPRESSURE 
      DO I=1,1000000
          READ (IUAPR,ERR=195,END=195) TIME1
          READ (IUAPR) TXX 
          WRITE (IUT89) TIME1
          WRITE (IUT89) TXX
      ENDDO
195   CONTINUE
#endif
***********************************************************************
*                                                                     *
*                                SED                                  *
*                                                                     *
***********************************************************************
#ifdef MODULE_SED
      READ (IUSDFLX,*)
      READ (IUSDFLX,*) NUMQBC
      ALLOCATE (VSDDIST(NUMQBC,KB))
      ALLOCATE (SEDDIS(NUMQBC))
      ALLOCATE (DSEDDIS(NUMQBC,2))           
      IF (NUMQBC.NE.0) THEN
      DO 125 N = 1, NUMQBC
          READ (IUSDFLX,*) NQD(N),NQC(N),(VSDDIST(N,K),K = 1,KBM1)
          WRITE (IUPRT,6000) NQD(N),NQC(N),(VSDDIST(N,K),K = 1,KBM1)
125   CONTINUE
      DO 135 I = 1, 100000
          READ (IUSDFLX,*,ERR=145,END=145) TIME
          WRITE (IUPRT,6402) TIME
          READ (IUSDFLX,*) (SEDDIS(N),N = 1,NUMQBC)
          WRITE (IUPRT,6400) (SEDDIS(N),N = 1,NUMQBC)
          WRITE (IUT20,6500) TIME
          WRITE (IUT20,6500) (SEDDIS(N),N = 1,NUMQBC)
135   CONTINUE
145   BACKSPACE IUSDFLX
      ENDIF
      CLOSE(IUSDFLX)
      READ (IUSDBC,*)
      READ (IUSDBC,*)
      IF (NUMEBC.NE.0) THEN
      DO 109 I = 1, 100000
          READ (IUSDBC,*,ERR=119,END=119) TIME1
          WRITE (IUPRT,6402) TIME1
          WRITE (IUT21,6500) TIME1
          DO 99 N = 1, NUMEBC
              READ (IUSDBC,*) (SEDBDRYSL(N,K),K = 1,KBM1)
              WRITE (IUPRT,5900) (SEDBDRYSL(N,K),K = 1,KBM1)
              WRITE (IUT21,6500) (SEDBDRYSL(N,K),K = 1,KBM1)   !WUHUI
99        CONTINUE
109   CONTINUE
119   BACKSPACE IUSDBC
      CLOSE(IUSDBC)
      ENDIF
#endif
***********************************************************************
*                                                                     *
*                               WAVE                                  *
*                                                                     *
***********************************************************************
#ifdef MODULE_WAVE
      DO 8010 M=1,1000000
          READ (IUWAVE,Err=8100,END=8100) TIME1
          READ (IUWAVE) WTEMP1   ! WAVE HEIGHT
          READ (IUWAVE) WTEMP2   ! WAVE PERIOD
          READ (IUWAVE) WTEMP3   ! WAVE DIRECTION
          WRITE (IUT92) TIME1
          WRITE (IUT92) (WTEMP1(I),I=1,N_CTRD)
          WRITE (IUT92) (WTEMP2(I),I=1,N_CTRD)
          WRITE (IUT92) (WTEMP3(I),I=1,N_CTRD)
 8010   CONTINUE           
 8100   CLOSE (IUWAVE)
#endif
***********************************************************************
*                                                                     *
*                               MATRIAL                               *
*                                                                     *
***********************************************************************
#ifdef MODULE_MATERIAL
c --- V2410:
      READ (IUMAT,5000) (COM(I),I = 1,80)
      READ (IUMAT,*) NUMMAT
      IF (NUMMAT.NE.0) THEN
          ALLOCATE (I_MAT(NUMMAT))
          ALLOCATE (XR_MAT(NUMMAT))
          ALLOCATE (YR_MAT(NUMMAT))
          ALLOCATE (VQDIST_MAT(NUMMAT,KB))
          ALLOCATE (MDIS(NUMMAT))
          ALLOCATE (DMDIS(NUMMAT,2))
      ENDIF
      IF (NUMMAT.NE.0) THEN 
      DO N=1,NUMMAT          
          READ(IUMAT,*) XR_MAT(N),YR_MAT(N),(VQDIST_MAT(N,K),K=1,KBM1)
      ENDDO
c --- V2410.
      DO 8020 M=1,1000000
          READ (IUMAT,*,Err=8200,END=8200) TIME1
          WRITE (IUT98,6500) TIME1
          READ (IUMAT,*) (MDIS(I),I=1,NUMQBC)
          WRITE (IUT98,6500) (MDIS(I),I=1,NUMQBC)
 8020   CONTINUE           
 8200   CLOSE (IUMAT)
      ENDIF !V2410
#endif        
      RETURN      
5000   FORMAT (80A1)
5100   FORMAT (/1X,80A1/)
5200   FORMAT (' KSL = ',I5,/)
5300   FORMAT (//' NUMBER OF STANDARD LEVELS IN RUN_DATA',I5,' (IKSL)'/
     *    '           DO NOT EQUAL'/
     *    ' NUMBER OF STANDARD LEVELS IN NOCTOM_INC ',I5,' (KSL)'/
     *    ' PLEASE CORRECT THIS PROBLEM AND TRY AGAIN'//)
5400   FORMAT (I5)
5500   FORMAT (2I5,4F10.5)
5600   FORMAT (4I5,1F10.5)
5601   FORMAT (2I7,1F10.5)
5700   FORMAT (8F10.5)
5800   FORMAT (16I5)
5801   FORMAT (16I7)
5900   FORMAT (16F5.1)
6000   FORMAT (4I5,20F5.2)
6001   FORMAT (2I7,20F5.2)
6100   FORMAT (4I5,/,20F5.0)
6200   FORMAT (2I5,20F5.2)
6300   FORMAT (2I5,/,20F5.0)
6400   FORMAT (8F10.5)
6402   FORMAT (F11.5)
6401   FORMAT (8F9.2)
6500   FORMAT (8E14.7)
6501   Format (50F14.7)
6502   Format (50F14.7)
6405   FORMAT (1000F10.5)  
       END SUBROUTINE BCDATA
