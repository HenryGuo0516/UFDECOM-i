#include "DEFS.h"
	
	SUBROUTINE FIRST

      USE MOD_GLOBAL
      
      REAL DIS,DISMIN !V2410      
      DIMENSION TXX(N_CTRD), TYY(N_CTRD)
      
      IF (NUMEBC.NE.0) THEN
      IF (TRIM(OPTEBC).EQ.'data') THEN
      REWIND IUT90
      DO 20 I = 1, 100000
          READ (IUT90,5100,ERR=170) T2E
          READ (IUT90,5100) (DEBDRY(N,2),N = 1,NUMEBC)
          IF (THOUR.LT.T2E) GO TO 30
          T1E = T2E
          DO 10 N = 1, NUMEBC
              DEBDRY(N,1) = DEBDRY(N,2)
10        CONTINUE
20    CONTINUE
30    CONTINUE
      ENDIF
	ENDIF
      
#ifdef ELB
      IF (NUMEBC.NE.0) THEN
      REWIND IUT95
      DO 21 I = 1, 100000
          READ (IUT95,5100,ERR=170) T2EL
          READ (IUT95,5100) (DELBC(N,2),N = 1,NUMEBC)
          IF (THOUR.LT.T2EL) GO TO 31
          T1EL = T2EL
          DO 11 N = 1, NUMEBC
              DELBC(N,1) = DELBC(N,2)
11        CONTINUE
21    CONTINUE
31    CONTINUE
      ENDIF
#endif


#if defined MODULE_SAL ||  defined MODULE_TMP
      IF (NUMEBC.NE.0) THEN
      REWIND IUT94
      DO 70 I = 1, 100000
          READ (IUT94,5100,ERR=170) T2TS
          DO 40 N = 1, NUMEBC
              READ (IUT94,5100) (DTBDRY(N,K,2),K = 1,KBM1)
              READ (IUT94,5100) (DSBDRY(N,K,2),K = 1,KBM1)
40        CONTINUE
          IF (THOUR.LT.T2TS) GO TO 80
          T1TS = T2TS
          DO 60 N = 1, NUMEBC
          DO 50 K = 1, KBM1
              DTBDRY(N,K,1) = DTBDRY(N,K,2)
              DSBDRY(N,K,1) = DSBDRY(N,K,2)
50        CONTINUE
60        CONTINUE
70    CONTINUE
80    CONTINUE
      ENDIF
#endif

      IF (NUMQBC.NE.0) THEN
      REWIND IUT91
      DO 100 I = 1, 100000
          READ (IUT91,5100,ERR=170) T2Q
          READ (IUT91,5100) (DQDIS(N,2),N = 1,NUMQBC)
          READ (IUT91,5100) (DTDIS(N,2),N = 1,NUMQBC)
          READ (IUT91,5100) (DSDIS(N,2),N = 1,NUMQBC)
          IF (THOUR.LT.T2Q) GO TO 110
          T1Q = T2Q
          DO 90 N = 1, NUMQBC
              DQDIS(N,1) = DQDIS(N,2)
              DTDIS(N,1) = DTDIS(N,2)
              DSDIS(N,1) = DSDIS(N,2)
90        CONTINUE
100   CONTINUE
      ENDIF
110   CONTINUE
      REWIND IUT93
#ifdef WDTYP_UNIFORM
      DO 150 II = 1, 100000
          READ (IUT93,ERR=170) T2M
          READ (IUT93) DQPREC(2), DQEVAP(2), DTX(2), DTY(2),DHFLUX(2)
          IF (THOUR.LT.T2M) GO TO 160
          T1M = T2M
          DQPREC(1) = DQPREC(2)
          DQEVAP(1) = DQEVAP(2)
          DTX(1) = DTX(2)
          DTY(1) = DTY(2)
        DHFLUX(1) = DHFLUX(2)
150   CONTINUE
#elif defined WDTYP_FIELD
      DO 151 II = 1, 100000
          READ (IUT93,ERR=170) T2M
          READ (IUT93) TXX,TYY
          DO I = 1,N_CTRD
              DTXX(I,2) = TXX(I)
              DTYY(I,2) = TYY(I)
          ENDDO
          IF (THOUR.LT.T2M) GO TO 160
          T1M=T2M
          DO I = 1,N_CTRD
              DTXX(I,1) = DTXX(I,2)
              DTYY(I,1) = DTYY(I,2)
          ENDDO
151   CONTINUE
#else
	PRINT*, 'Wind data type undefined!'
	PAUSE
	STOP
#endif
160   CONTINUE
#ifdef HEATFLUX_BULK
      REWIND (IUT86)
      DO II = 1, 100000
          READ (IUT86,END=170,ERR=170) T2A
          READ (IUT86) ATP
          DO I=1,N_CTRD
              DATP(I,2)=ATP(I)
          ENDDO
          IF (THOUR.LT.T2A) GO TO 161
          T1A=T2A
          DO I=1,N_CTRD
              DATP(I,1)=DATP(I,2)
          ENDDO
      ENDDO
161   CONTINUE
      REWIND (IUT87)
      DO II = 1, 100000
          READ (IUT87,END=170,ERR=170) T2R
          READ (IUT87) RHM
          DO I=1,N_CTRD
              DRHM(I,2)=RHM(I)
          ENDDO
          IF (THOUR.LT.T2R) GO TO 162
          T1R=T2R
          DO I=1,N_CTRD
              DRHM(I,1)=DRHM(I,2)
          ENDDO
      ENDDO
162   CONTINUE
      REWIND (IUT88)
      DO II = 1, 100000
          READ (IUT88,END=170,ERR=170) T2C
          READ (IUT88) CLOUD
          DO I=1,N_CTRD
              DCLOUD(I,2)=CLOUD(I)
          ENDDO
          IF (THOUR.LT.T2C) GO TO 163
          T1C=T2C
          DO I=1,N_CTRD
              DCLOUD(I,1)=DCLOUD(I,2)
          ENDDO
      ENDDO
163   CONTINUE
#endif
#if defined HEATFLUX_BULK  ||  defined AIRPRESSURE 
      REWIND (IUT89)
      DO II = 1, 100000
          READ (IUT89,END=170,ERR=170) T2APR
          READ (IUT89) APR
          DO I=1,N_CTRD
              DAPR(I,2)=APR(I)
          ENDDO
          IF (THOUR.LT.T2APR) GO TO 164
          T1APR=T2APR
          DO I=1,N_CTRD
              DAPR(I,1)=DAPR(I,2)
          ENDDO
      ENDDO
164   CONTINUE
#endif
   
   
#ifdef MODULE_SED
      IF (NUMQBC.NE.0) THEN
      REWIND IUT20
      DO 503 I = 1, 100000
          READ (IUT20,*,ERR=170,END=170) T2SEDQ
          READ (IUT20,*) (DSEDDIS(N,2),N = 1,NUMQBC)
          IF (THOUR.LT.T2SEDQ) GO TO 513
          T1SEDQ = T2SEDQ
          DO 93 N = 1, NUMQBC
              DSEDDIS(N,1) = DSEDDIS(N,2)
93        CONTINUE
503   CONTINUE
      ENDIF
513	CONTINUE
      IF (NUMEBC.NE.0) THEN
      REWIND IUT21
      DO 571 I = 1, 100000
          READ (IUT21,*,ERR=170,END=170) T2SED
          DO 540 N = 1, NUMEBC
              READ (IUT21,*) (DSEDBDRY(N,K,2),K = 1,KBM1)
540       CONTINUE
          IF (THOUR.LT.T2SED) GO TO 580
          T1SED = T2SED
          DO 563 N = 1, NUMEBC
          DO 550 K = 1, KBM1
              DSEDBDRY(N,K,1) = DSEDBDRY(N,K,2)
550       CONTINUE
563       CONTINUE
571   CONTINUE
580   CONTINUE
	ENDIF 
#endif

#ifdef MODULE_WAVE
      REWIND IUT92
      DO 681 K = 1, 100000  
          READ (IUT92,ERR=170,END=170) T2WAVE
          READ (IUT92) (DWHT (I,2),I=1,N_CTRD)
          READ (IUT92) (DWPER(I,2),I=1,N_CTRD)
          READ (IUT92) (DWDIR(I,2),I=1,N_CTRD)   
          IF (THOUR.Lt.T2WAVE) GO TO 690
          T1WAVE = T2WAVE
681   CONTINUE
690   CONTINUE 
#endif
      ISTART = NSTEP + 1
      IEND = 100000000
      THOUR = TDAY * 24.
	SECOND = THOUR*3600.
	NN_OPT = INT(SECOND/N_OPT)
      NN_OPT2 = INT(SECOND/N_OPT) !V2410
#ifdef RSTARC
	NN_OPTRST = INT(TDAY/RST_ARCH)
#endif
#ifdef TIDE_FLATHER
      REWIND IUT96
      DO 108 I = 1, 100000
          READ (IUT96,*,END=170,ERR=170) T2F
          DO N=1,NUMEBC
              READ (IUT96,*) (DVFOBC(N,K,2), K=1,KBM1)
          ENDDO
          IF (THOUR.LT.T2F) GO TO 117
          T1F = T2F
          DO N = 1, NUMEBC
          DO K=1,KBM1
              DVFOBC(N,K,1) =DVFOBC(N,K,2)
          ENDDO
          ENDDO
108   CONTINUE
117   CONTINUE      
#endif      
#ifdef FPT
	NN_FPT = INT(SECOND/N_FPT)
#endif   
      WRITE (IUGDH) XR,YR,LON,LAT,H,DJ
	CLOSE (IUGDH)
#ifdef MODULE_MATERIAL      
      REWIND IUT98
      DO N=1,NUMMAT !V2410
          DISMIN=50000.
          DO I=1,N_CTRD_AG
              IF (FSM(I).GT.0) THEN
                  DIS=SQRT((XR(I)-XR_MAT(N))**2+(YR(I)-YR_MAT(N))**2)
                  IF (DIS.LT.DISMIN) THEN
                      DISMIN=DIS
                      I_MAT(N)=I
                  ENDIF
              ENDIF
          ENDDO
      ENDDO
      
      DO 111 I = 1, 100000
          READ (IUT98,5100,ERR=170) T2MAT
          READ (IUT98,5100) (DMDIS(N,2),N = 1,NUMQBC)
          IF (THOUR.LT.T2MAT) GO TO 112
          T1MAT = T2MAT
          DO 91 N = 1, NUMQBC
              DMDIS(N,1) = DMDIS(N,2)
91        CONTINUE
111   CONTINUE
112   CONTINUE
#endif      
      RETURN      
170   WRITE (6,5000)
	PAUSE
      STOP
5000  FORMAT (//' THERE IS INSUFFICIENT TEMPORAL DATA FOR THIS RUN'/,
     *    '        REVISE INPUT DECK AND RESUBMIT '//)
5100  FORMAT (8E14.7)
      END SUBROUTINE FIRST
