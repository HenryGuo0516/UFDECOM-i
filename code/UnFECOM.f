#include "DEFS.h"	
	PROGRAM UnFECOM
      
      USE ieee_arithmetic      
      USE DATETIME1
      USE MOD_GLOBAL
	USE MOD_OUTPUT
#ifdef MODULE_SED
      USE MOD_SED
#endif
#ifdef WEIR
      USE MOD_WEIR
#endif
#ifdef MODULE_WAVE
      USE MOD_WAVE
#endif
#ifdef MODULE_LAG
      USE MOD_LAG !V2410
#endif
      USE ADVT
      USE TGSpeed !CBR

      IMPLICIT NONE
      
      INTEGER :: DATE_TIME(8)
      INTEGER I,J,K,KK

      REAL SECOND,UV_ABS
c      REAL TG1,TG2,TG3,TG4,TG5,TG6,TG7,TG8,TG9,TG10,TGA
c      REAL TGA1,TGA2,TGA3,TGA4,TGA5
c      INTEGER NGA1,NGA2,NGA3,NGA4,NGA5
      
c      CHARACTER (LEN = 12) REAL_CLOCK (3)
c      INTEGER(KIND=8) ::COUNT1,COUNT2,COUNT_RATE,COUNT_MAX
      TYPE(TIME) TIME00,TIME01 !V2410
      
      LOG_RSTSPEC = .FALSE.
      
	CALL READ_NML      
	CALL OPENFILE(1)
	CALL ALLOC_VARS
      CALL BCDATA
      CALL SETDOM 
      CALL OPENFILE(2)
      CALL FIND_INTERFACE_AG
      CALL DISTINGUISH
	CALL PRE_IE
	CALL PRE_TSR
      CALL VARBD(1)
      CALL VARBD(11) !CBR: EBC FSMADD, HU, HV moved from VARBD(1) to VARBD(11)
      CALL RESTART
      CALL GEN_N
#ifdef Z0BTYP_CONSTANT
      CALL BROUGH
#elif defined Z0BTYP_USERDEF
      CALL BROUGH_USERDEFINED
#else
      PRINT*, 'Bottom drag coefficient type undefined!'
      PAUSE
      STOP
#endif
      CALL FIRST
c   --- V2410:
      TIME00=NEWTIME(IYEAR,IMONTH,IDAY0,0,0,0)
      TIME01=NEWTIME(YEAR_TH,MONTH_TH,DAY_TH,HOUR_TH,
     *MINUTE_TH,SECOND_TH)
      TG_TH=DBLE(TIME01-TIME00)
      
c      Print*, NNIVGFEALL,NNIVGFWALL,NNIVGFNALL,NNIVGFSALL
c      pause
      
c --- CBR: Time tests variables:      
      NGA1=0; NGA2=0; NGA3=0; NGA4=0; NGA5=0
      TG1=0.0; TG2=0.0; TG3=0.0; TG4=0.0; TG5=0.0; TG6=0.0
      TG7=0.0; TG8=0.0; TG9=0.0; TG10=0.0; TGA=0.0
      TGA1=0.0; TGA2=0.0; TGA3=0.0; TGA4=0.0; TGA5=0.0
      TG31=0.0; TG32=0.0; TG33=0.0; TG34=0.0; TG35=0.0; TG36=0.0 
*******************************************************************************
*                            ENTER DATA TO DEVICE                             *
*******************************************************************************
*  =============== Constant Arrays:
!$ACC ENTER DATA COPYIN(DZR,Z,ZZ,DZ,DZZ,H,HU,HV,H1,H2,H3)     
!$ACC ENTER DATA COPYIN(DX1,DX2,DX3,DX4,DY1,DY2,DY3,DY4)     
!$ACC ENTER DATA COPYIN(DJ,ARU,ARV,COR)     
!$ACC ENTER DATA COPYIN(H,H1,H2,H3)     
!$ACC ENTER DATA COPYIN(NIM1,NIM2,NIM3,NIM4)     
!$ACC ENTER DATA COPYIN(NIP1,NIP2,NIP3,NIP4)     
!$ACC ENTER DATA COPYIN(NJM1,NJM2,NJM3,NJM4)     
!$ACC ENTER DATA COPYIN(NJP1,NJP2,NJP3,NJP4)     
!$ACC ENTER DATA COPYIN(NIM1JP1,NIM1JM1,NIP1JP1,NIP1JM1)     
!$ACC ENTER DATA COPYIN(NIM1JM2,NIM2JM1)     

!$ACC ENTER DATA COPYIN(NETA,NCON,NAJA,NQC,NQD,VQDIST)     
!$ACC ENTER DATA COPYIN(NIVGFWALL,NIVGFEALL)
!$ACC ENTER DATA COPYIN(NIVGFSALL,NIVGFNALL)     
!$ACC ENTER DATA COPYIN(XNODE,YNODE,XR,YR)     
!$ACC ENTER DATA COPYIN(WTU_IVG,WTV_IVG,WTC_IVG)     
!$ACC ENTER DATA COPYIN(WTU_EVG,WTV_EVG,WTC_EVG)     
!$ACC ENTER DATA COPYIN(NAGU_IVG,NAGV_IVG,NAGC_IVG)     
!$ACC ENTER DATA COPYIN(NAGU_EVG,NAGV_EVG,NAGC_EVG)     
!$ACC ENTER DATA COPYIN(FSMADD,KM_COF,KH_COF)     
!$ACC ENTER DATA COPYIN(NAG4VGC,IAPOI) 
!$ACC ENTER DATA COPYIN(CBC_COF_A1,CBC_COF_B1)
!$ACC ENTER DATA COPYIN(CBC_COF_A2,CBC_COF_B2,CBCADJ)
!$ACC ENTER DATA COPYIN(XETA,XXI,YETA,YXI) 
c!$ACC ENTER DATA COPYIN() 
      
*  =============== Variable Arrays:
!$ACC ENTER DATA COPYIN(MASK4VAR)     
!$ACC ENTER DATA COPYIN(XPOI,RHSPOI,AIJ,APOI) 
c!$ACC ENTER DATA COPYIN(JAPOI,DPOI) !CBR: DPOI scalarized
!$ACC ENTER DATA COPYIN(KAPOI)
      
!$ACC ENTER DATA COPYIN(CBC_U,CBC_V,UC,VC) 
      
*  ----- Frequently Required:
!$ACC ENTER DATA COPYIN(FSM,FSM11,DUM,DVM,U,V)
!$ACC ENTER DATA COPYIN(UF,VF,EL,ELF,CURV4)
!$ACC ENTER DATA COPYIN(UR,VR,UNN,VNN,W,WNN,WR)     
!$ACC ENTER DATA COPYIN(D,DU,DV)
!$ACC ENTER DATA COPYIN(Q2,Q2L)     
!$ACC ENTER DATA COPYIN(RHO,DRHOX,DRHOY)     
!$ACC ENTER DATA COPYIN(KAX,KAY,KAZ,S,T,SWRAD,RAD)     
!$ACC ENTER DATA COPYIN(BOYGR,KN,GH,SH,SM)     
!$ACC ENTER DATA COPYIN(KM,KH,L,KQ)     
!$ACC ENTER DATA COPYIN(WTSURF)     
     
      
*  ----- Updated in SMAG.f:
!$ACC ENTER DATA COPYIN(AAM)
*  ----- Updated in eltest_sor.f:
!$ACC ENTER DATA COPYIN(AAUU,AAVV)
!$ACC ENTER DATA COPYIN(XMFLUX,YMFLUX,XMFLUX0,YMFLUX0)
!$ACC ENTER DATA COPYIN(XFLUX,YFLUX,DTEF,PROD) 
!$ACC ENTER DATA COPYIN(A,C,VH,VHP,TPS) 

*  ----- Updated to device in BCOND(9):
!$ACC ENTER DATA COPYIN(DEBDRY,DELBC,DTBDRY,DSBDRY)
!$ACC ENTER DATA COPYIN(EBDRY,ELBC,TBDRY,SBDRY)
!$ACC ENTER DATA COPYIN(DQDIS,DTDIS,DSDIS,DMDIS) !V2410
!$ACC ENTER DATA COPYIN(QDIS,TDIS,SDIS,MDIS) !V2410
!$ACC ENTER DATA COPYIN(I_MAT,VQDIST_MAT) !V2410
!$ACC ENTER DATA COPYIN(CONCENTRATION) !V2410
!$ACC ENTER DATA COPYIN(DTXX,DTYY)
!$ACC ENTER DATA COPYIN(WUSURF,WVSURF,WINDU,WINDV)
!$ACC ENTER DATA COPYIN(WUBOT,WVBOT)
!$ACC ENTER DATA COPYIN(DATP,DRHM,DCLOUD,DAPR,DVFOBC)
!$ACC ENTER DATA COPYIN(ATP,RHM,CLOUD,APR,VFOBC)

*  ----- for INFO_EXCH_GPU.f:
!$ACC ENTER DATA COPYIN(MASK4VAR,XIU,XIV,XIC)
*  ----- for bcond(3):
!$ACC ENTER DATA COPYIN(ELBC_TIDE)
*  ----- for eltest_sor.f:
!$ACC ENTER DATA COPYIN(NIAGCW,NIAGCE,NIAGCN,NIAGCS)
!$ACC ENTER DATA COPYIN(NIVGFW,NIVGFE,NIVGFN,NIVGFS)
!$ACC ENTER DATA COPYIN(NIVGFWS,NIVGFWN,NIVGFES,NIVGFEN)
!$ACC ENTER DATA COPYIN(NVGF4AG)
      
      PRINT*, '------ BEGIN NUMERICAL INTEGRATION ------'
*******************************************************************************
*                                                                             *
*                         BEGIN NUMERICAL INTEGRATION                         *
*                                                                             *
******************************************************************************* 
      DO 160 NSTEP=ISTART,IEND
	IF (TDAY .GT. T_END) THEN
	    GOTO 170
      ELSE
******************************************************************************* 
c -------------- Record Speed:
      CALL SYSTEM_CLOCK(COUNT1,COUNT_RATE) !ver 02G moved here
***************************** ADJUST THE TIME STEP ****************************
      CALL TIMESTEP_GPU
c      CALL TIMESTEP_MC
c      CALL TIMESTEP
************************ REMEMBER INFORMATION OF STEP_N ***********************
      TDAY0 = TDAY
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 1, KB
c      DO I = 1, N_CTRDP1
      do concurrent (I = 1:N_CTRDP1,K = 1:KB)
          XMFLUX0(I,K) = XMFLUX(I,K)
          YMFLUX0(I,K) = YMFLUX(I,K)
      ENDDO
c      ENDDO
c!$OMP END TARGET      
      
      COUNT_LOOP=0
********************************* ONLINE VIEW *********************************  
#ifdef MOD_ONLINE_VIEW
      IF (NSTEP.EQ.ISTART) CALL INITVIEW
      IF(MOD(NSTEP,50).EQ.0.OR.NSTEP.EQ.ISTART) THEN
          CALL XVIEW(CHC_VIEW)
      ENDIF	
#endif
******************************* CALCULATE AGAIN *******************************
1011  CONTINUE
      IFCOUNT=0
c      XMFLUX = XMFLUX0 !CBR: moved to before goto 1011
c      YMFLUX = YMFLUX0 !CBR: moved to before goto 1011
      TDAY = TDAY0+DTI*DAYI
      THOUR = TDAY*24.
      SECOND=THOUR*3600.
      RAMP = TANH(FLOAT(NSTEP)/FLOAT(IRAMP+1))
      RAMP_TIDE = TANH((TDAY-TIDE_LAG)/0.5)  
******************************* MOVING BOUNDARY *******************************
c -------------- Record Speed:
c      CALL SYSTEM_CLOCK(COUNT1,COUNT_RATE) !ver 02G moved to beginning
      CALL BCOND(9)
c --- Variables below are updated from HOST to DEVICE in BCOND(9):      
c!$ACC UPDATE DEVICE (EBDRY,ELBC,TBDRY,SBDRY,QDIS,TDIS,SDIS,
c     *WUSURF,WVSURF,WINDU,WINDV,ATP,RHM,CLOUD,APR,VFOBC)
c --- Target both EBDRY and DEBDRY might be better, yet to test.
      
c -------------- Record Speed:
      CALL SYSTEM_CLOCK(COUNT2,COUNT_RATE)
      TG1=TG1+(COUNT2-COUNT1)/REAL(COUNT_RATE)
#ifdef MODULE_SED
      IF (THOUR .GE. SED_BEG) THEN
          CALL BOTSTRESS
      ENDIF
#endif 
***************************** CALCULATE BAROCLINE *****************************

#ifdef BPG
        IF(THOUR .GT. BPG_BEG) THEN
            CALL BAROPG5B
c!$ACC UPDATE DEVICE (DRHOX,DRHOY)
        ENDIF
#endif
c -------------- Record Speed:
      CALL SYSTEM_CLOCK(COUNT1,COUNT_RATE)
      TG9=TG9+(COUNT1-COUNT2)/REAL(COUNT_RATE)
*********************** CALCULATE DIFFUSION COEFFICIENT  **********************
#ifdef HORZMIX_CLOSURE
	CALL SMAG
#endif
c -------------- Record Speed:
      CALL SYSTEM_CLOCK(COUNT2,COUNT_RATE)
      TG10=TG10+(COUNT2-COUNT1)/REAL(COUNT_RATE)
#ifdef WEIR
      CALL ADDWEIR(1)
#endif
************************* CALCULATE MOMENTUM EQUATION *************************
************************************ STEP_H ***********************************
      CALL ADVU_TVD_3RD
      CALL ADVV_TVD_3RD
c -------------- Record Speed:
      CALL SYSTEM_CLOCK(COUNT1,COUNT_RATE)
      TG2=TG2+(COUNT1-COUNT2)/REAL(COUNT_RATE)
************************************ STEP_V ***********************************
      CALL PROFU
      CALL PROFV
*******************************************************************************
      VTP = 2
      CALL INFO_EXCH_GPU(VTP,1,WUBOT,DUM
     * ,NAGU_EVG,WTU_EVG,NAGU_IVG,WTU_IVG)
      VTP = 3
      CALL INFO_EXCH_GPU(VTP,1,WVBOT,DVM
     * ,NAGV_EVG,WTV_EVG,NAGV_IVG,WTV_IVG)
c -------------- Record Speed:
      CALL SYSTEM_CLOCK(COUNT2,COUNT_RATE)
      TG3=TG3+(COUNT2-COUNT1)/REAL(COUNT_RATE)
********************************** ELEVATION **********************************
#ifdef EXP_EL
      CALL CAL_EL
      VTP = 1
      CALL INFO_EXCH(VTP,1,ELF,FSM
     * ,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
#elif defined IMP_EL
      CALL ELTEST
      VTP = 1
      CALL INFO_EXCH_GPU(VTP,1,ELF,FSM
     * ,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
      CALL BCOND(1)
#endif
c -------------- Record Speed:
      CALL SYSTEM_CLOCK(COUNT1,COUNT_RATE)
      TG4=TG4+(COUNT1-COUNT2)/REAL(COUNT_RATE)

************************************ CHECK ************************************
c      DO I = 1,N_CTRD_AG
c      UV_ABS = SQRT(UF(I,1)**2+VF(I,1)**2)
c      IF (UV_ABS .GE.15..OR.ISNAN(UV_ABS).OR.ABS(ELF(I)).GT.10) THEN
c          WRITE (IUUVK,*) NSTEP-1,DTI
c          WRITE (IUUVK,*) I,UV_ABS,H(I),ELF(I)
c          WRITE (IUUVK,*) (UF(I,K),K=1,KBM1)
c          WRITE (IUUVK,*) (VF(I,K),K=1,KBM1)
c          WRITE (IUUVK,*) (KM(I,K),K=1,KBM1)
c          IFCOUNT=1
c      ENDIF
c      ENDDO
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO reduction(.OR.:IFCOUNT)
c      DO I = 1,N_CTRD_AG
      DO concurrent (I = 1:N_CTRD_AG) REDUCE(.OR.:IFCOUNT) local(UV_ABS)
c      DO concurrent (I = 1:N_CTRD_AG) local(UV_ABS)
      UV_ABS = UF(I,1)**2+VF(I,1)**2
!      IF (UV_ABS .GE.225..OR.ISNAN(UV_ABS).OR.ABS(ELF(I)).GT.10) THEN
      IF (UV_ABS>=225..OR.IEEE_IS_NAN(UV_ABS).OR.ABS(ELF(I))>20.) THEN !V2410:10->20
!          WRITE (IUUVK,*) NSTEP-1,DTI
!          WRITE (IUUVK,*) I,UV_ABS,H(I),ELF(I)
!          WRITE (IUUVK,*) (UF(I,K),K=1,KBM1)
!          WRITE (IUUVK,*) (VF(I,K),K=1,KBM1)
!          WRITE (IUUVK,*) (KM(I,K),K=1,KBM1)
          IFCOUNT=1
      ENDIF
      ENDDO
c!$OMP END TARGET 
      
      IF (IFCOUNT.EQ.1) THEN
!$ACC UPDATE HOST (UF,VF,ELF,KM) !GPU to CPU
      DO I = 1,N_CTRD_AG
      UV_ABS = SQRT(UF(I,1)**2+VF(I,1)**2)
!      IF (UV_ABS .GE.15..OR.ISNAN(UV_ABS).OR.ABS(ELF(I)).GT.10) THEN
      IF (UV_ABS>=15..OR.IEEE_IS_NAN(UV_ABS).OR.ABS(ELF(I))>20.) THEN !V2410:10->20
          WRITE (IUUVK,*) NSTEP-1,DTI
          WRITE (IUUVK,*) I,UV_ABS,H(I),ELF(I)
          WRITE (IUUVK,*) (UF(I,K),K=1,KBM1)
          WRITE (IUUVK,*) (VF(I,K),K=1,KBM1)
          WRITE (IUUVK,*) (KM(I,K),K=1,KBM1)
!          IFCOUNT=1
      ENDIF
      ENDDO
      DTI = DTI*0.5
      COUNT_LOOP = COUNT_LOOP+1
      IF (COUNT_LOOP.GE.2) THEN
          PRINT*, 'SIMULATION BECOMES UNSTABLE, PLEASE TRY
     *VARIABLE TIMESTEP OR REDUCE THE TIMESTEP'
          PAUSE
          STOP
      ENDIF
c!$OMP TARGET DEFAULTMAP(present: allocatable)     
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 1, KB
c      DO I = 1, N_CTRDP1
      do concurrent (I = 1:N_CTRDP1,K = 1:KB)
          XMFLUX(I,K) = XMFLUX0(I,K) !CBR: moved to before goto 1011
          YMFLUX(I,K) = YMFLUX0(I,K) !CBR: moved to before goto 1011
      ENDDO
c      ENDDO
c!$OMP END TARGET      
      GOTO 1011
      ENDIF


*******************************************************************************
      CALL REUV
      CALL VARBD_GPU(1)
      
#ifdef Z0BTYP_CONSTANT
      CALL BROUGH
#elif defined Z0BTYP_USERDEF
      CALL BROUGH_USERDEFINED_GPU
#else
      PRINT*, 'Bottom drag coefficient type undefined!'
      PAUSE
      STOP
#endif
      
************************** END: COMPUTE UF AND VF *****************************
      CALL VERTVL
      CALL BCOND(2)          
c      CALL WREAL  !CBR: moved to before OUTPUT_FIELD
      VTP = 0
!      DO K = 1,KB
!          CALL INFO_EXCH(VTP,K,W(:,K),FSM
!     *,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
!          CALL INFO_EXCH(VTP,K,WR(:,K),FSM
!     *,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
!      ENDDO
          CALL INFO_EXCH3D_GPU(W,FSM
     *,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
c          CALL INFO_EXCH3D_GPU(WR,FSM
c     *,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)  !CBR: moved to before OUTPUT_FIELD
c -------------- Record Speed:
c      CALL SYSTEM_CLOCK(COUNT1,COUNT_RATE)
c      TG4=TG4+(COUNT1-COUNT2)/REAL(COUNT_RATE)
c -------------- Record Speed:
      CALL SYSTEM_CLOCK(COUNT2,COUNT_RATE)
      TG5=TG5+(COUNT2-COUNT1)/REAL(COUNT_RATE)
      
********* COMPUTE Q2 AND Q2L USING UF AND VF AS TEMPORARY VARIABLES ***********
#ifdef VERTMIX_CLOSURE
      CALL ADVQ(Q2,UF)
      CALL ADVQ(Q2L,VF)
c -------------- Record Speed:
c      CALL SYSTEM_CLOCK(COUNT2,COUNT_RATE)
c      TG5=TG5+(COUNT2-COUNT1)/REAL(COUNT_RATE)
      CALL PROFQ
      CALL BCOND(5)
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 1,KB
c      DO I = 1,N_CTRD_AG
      do concurrent (I = 1:N_CTRD_AG,K = 1:KB)
          Q2(I,K) = UF(I,K)
          Q2L(I,K) = VF(I,K)
      ENDDO
c      ENDDO
c!$OMP END TARGET
      
      VTP = 0
!      DO K = 1,KB
!          CALL INFO_EXCH(VTP,K,Q2(:,K),FSMADD
!     *,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
!          CALL INFO_EXCH(VTP,K,Q2L(:,K),FSMADD
!     *,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
!          CALL INFO_EXCH(VTP,K,KM(:,K),FSMADD
!     *,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
!          CALL INFO_EXCH(VTP,K,KH(:,K),FSMADD
!     *,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
!      ENDDO 
          CALL INFO_EXCH3D_GPU(Q2,FSMADD
     *,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
          CALL INFO_EXCH3D_GPU(Q2L,FSMADD
     *,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
          CALL INFO_EXCH3D_GPU(KM,FSMADD
     *,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
          CALL INFO_EXCH3D_GPU(KH,FSMADD
     *,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 1,KBM1
c      DO I = 1,N_CTRD
      do concurrent (I = 1:N_CTRD,K = 1:KBM1) local_init(UMOL1)
!          IF (KM(I,K).LT.0..OR.ISNAN(KM(I,K))) KM(I,K) = UMOL1
!          IF (KH(I,K).LT.0..OR.ISNAN(KH(I,K))) KH(I,K) = UMOL1
          IF (KM(I,K).LT.0..OR.IEEE_IS_NAN(KM(I,K))) KM(I,K) = UMOL1
          IF (KH(I,K).LT.0..OR.IEEE_IS_NAN(KH(I,K))) KH(I,K) = UMOL1
      ENDDO
c      ENDDO
c!$OMP END TARGET
      
#endif
c -------------- Record Speed:
      CALL SYSTEM_CLOCK(COUNT1,COUNT_RATE)
      TG6=TG6+(COUNT1-COUNT2)/REAL(COUNT_RATE)
      
********** COMPUTE T AND S USING UF AND VF AS TEMPORARY VARIABLES *************
#ifdef MODULE_SAL
      IF (THOUR .GE. S_BEG) THEN
      CALL CAL_SAL
      
      VTP = 0
!      DO K = 1,KB
!      CALL INFO_EXCH(VTP,K,S(:,K),FSM
!     * ,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
!      ENDDO
      CALL INFO_EXCH3D_GPU(S,FSM
     * ,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
      ENDIF
c -------------- Record Speed:
      CALL SYSTEM_CLOCK(COUNT2,COUNT_RATE)
      TG7=TG7+(COUNT2-COUNT1)/REAL(COUNT_RATE)
#endif

#ifdef MODULE_TMP
      IF (THOUR .GE. T_BEG) THEN
      CALL CAL_TMP
      VTP = 0
!      DO K = 1,KB
!      CALL INFO_EXCH(VTP,K,T(:,K),FSM
!     * ,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
!      ENDDO
      CALL INFO_EXCH3D_GPU(T,FSM
     * ,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
      ENDIF
#endif
c -------------- Record Speed:
      CALL SYSTEM_CLOCK(COUNT1,COUNT_RATE)
      TG8=TG8+(COUNT1-COUNT2)/REAL(COUNT_RATE)
      
#ifdef MODULE_MATERIAL
      IF (THOUR .GE. M_BEG) THEN
      CALL CAL_MATERIAL
      VTP = 0
!      DO K = 1,KB
!      CALL INFO_EXCH(VTP,K,CONCENTRATION(:,K),FSM
!     * ,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
!      ENDDO
      CALL INFO_EXCH3D_GPU(CONCENTRATION,FSM
     * ,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG) !V2410
      ENDIF
#endif
*************************** SOLVE SEDIMENT PROCESS ****************************
#ifdef MODULE_SED
      IF (THOUR .GE. SED_BEG) THEN
      CALL CAL_SED
      VTP = 0
!      DO K = 1,KB
!          CALL INFO_EXCH(VTP,K,SED(:,K),FSM
!     * ,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
!      ENDDO
          CALL INFO_EXCH(SED,FSM
     * ,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
      ENDIF
#endif
*******************************************************************************
#if defined MODULE_SAL  ||  defined MODULE_TMP  ||  defined MODULE_SED  ||  defined MODULE_MATERIAL
      CALL VARBD_GPU(2)     
      CALL DENS_GPU
#endif

c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I = 1,N_CTRD !CBR: moved before CALL RESIDUAL
      do concurrent (I = 1:N_CTRD)
      IF(H(I) .GE. -10.) THEN
c           EL1(I) = EL(I)
           EL(I) = ELF(I)
      ENDIF
      ENDDO
c!$OMP END TARGET

      
      
*******************************************************************************
*                            OUTPUT SUBROUTINES                               *
*******************************************************************************
c -------------- Record Speed:
c      CALL SYSTEM_CLOCK(COUNT2,COUNT_RATE)
c      TGA1=TGA1+(COUNT2-COUNT1)/REAL(COUNT_RATE)


c -------------- Record Speed:
      CALL SYSTEM_CLOCK(COUNT2,COUNT_RATE)
      TGA1=TGA1+(COUNT2-COUNT1)/REAL(COUNT_RATE)

      CALL RESIDUAL
      CALL SEC
c -------------- Record Speed:
        CALL SYSTEM_CLOCK(COUNT1,COUNT_RATE)
        TGA2=TGA2+(COUNT1-COUNT2)/REAL(COUNT_RATE)
*******************************************************************************
	SECOND = THOUR*3600.
	IF ((INT(SECOND/N_OPT) .GT. NN_OPT) .OR. (NSTEP.EQ.ISTART)) THEN
      IF (N_TSR .NE. 0) THEN
          IF (RESTAR .EQ. 'cold') THEN
              CALL OUTPUT_POINT
              NGA1=NGA1+1 !CBR: Test output number
          ELSEIF ((RESTAR .EQ. 'hot') .AND. (THOUR .GT. THOUR_HOT)) THEN
              CALL OUTPUT_POINT
              NGA1=NGA1+1 !CBR: Test output number
          ENDIF
          IF (N_SEC_N.NE.0) THEN !      N_TSR .NE. 0
		    IF (RESTAR.EQ.'cold') THEN
		        CALL OUTPUT_SECFLUX
                  NGA2=NGA2+1 !CBR: Test output number
		    ELSE IF (RESTAR.EQ.'hot'.AND.THOUR.GT.THOUR_HOT) THEN
		        CALL OUTPUT_SECFLUX
                  NGA2=NGA2+1 !CBR: Test output number
		    ENDIF
		ENDIF
      ENDIF
      NN_OPT = NN_OPT+1
      ENDIF
#ifdef MODULE_LAG
      CALL PARTICLE_TRACK !V2410
#endif      
      
      
c -------------- Record Speed:
        CALL SYSTEM_CLOCK(COUNT2,COUNT_RATE)
        TGA3=TGA3+(COUNT2-COUNT1)/REAL(COUNT_RATE)
#ifdef FPT   
      IF (INT(SECOND/N_FPT) .GT. NN_FPT) THEN
      NN_FPT = NN_FPT+1
      IF ((TDAY .GE. FPTSTAR) .AND. (TDAY .LE. FPTEND)) THEN
          CALL WREAL  !CBR: moved here for WR
          CALL INFO_EXCH3D_GPU(WR,FSM
     *,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)  !CBR: moved here for WR
          CALL OUTPUT_FIELD
          NGA3=NGA3+1 !CBR: Test output number
      ENDIF
      ENDIF
#endif
c -------------- Record Speed:
      CALL SYSTEM_CLOCK(COUNT1,COUNT_RATE)
      TGA4=TGA4+(COUNT1-COUNT2)/REAL(COUNT_RATE)

      IF (MOD(NSTEP,N_RST).EQ.0) THEN 
          FN = 'RESTART'
          WRITE (IUPRT,8300) TDAY
#include "WRITE_RST.h"
          NGA4=NGA4+1 !CBR: Test output number
	ENDIF

#ifdef RSTARC
      IF (INT(TDAY/RST_ARCH).GT.NN_OPTRST) THEN
          NN_OPTRST = NN_OPTRST + 1
          WRITE(FN,553) NN_OPTRST
553   FORMAT('./RSTFILES/RESTART_',I5.5)
#include "WRITE_RST.h"
          NGA4=NGA4+1 !CBR: Test output number
	ENDIF
#endif
	IF (.NOT.LOG_RSTSPEC) THEN
      IF (TDAY.GT.RST_SPEC) THEN
          WRITE(FN,554) INT(TDAY)
554		FORMAT('RESTART_',I5.5) 
#include "WRITE_RST.h"
          NGA4=NGA4+1 !CBR: Test output number
		LOG_RSTSPEC = .TRUE.
      ENDIF
      ENDIF

      IF (LOG_TH) THEN 
          CALL CAL_TH !V2410
      ENDIF
      
c      WRITE (IUTSP,*) tday*24,LOG_SLUICE_ON !V2410
      
      
c -------------- Record Speed:
      CALL SYSTEM_CLOCK(COUNT2,COUNT_RATE)
      TGA5=TGA5+(COUNT2-COUNT1)/REAL(COUNT_RATE)
      IF (MOD(NSTEP,50).EQ.0) THEN
          TGA=TGA1+TGA2+TGA3+TGA4+TGA5
c --- CBR add print:          
       PRINT*, "NStep=",NSTEP
       WRITE (*,'(11F8.2," | ",5F8.2," | ",6F8.2)')
     *TG1,TG2,TG3,TG4,TG5,TG6,TG7,TG8,TG9,TG10,TGA,
     *TGA1,TGA2,TGA3,TGA4,TGA5,TG31,TG32,TG33,TG34,TG35,TG36
c       WRITE (*,*) TG1,TG2,TG3,TG4,TG5,TG6,TG7,
c     *TG8,TG9,TG10,TGA,TGA1,TGA2,TGA3,TGA4,TGA5
c       WRITE (*,'(4I8," | ",4(F8.2,"%"))') NGA1,NGA2,NGA3,NGA4,
c     * float(NGA1*100)/NSTEP,float(NGA2*100)/NSTEP,
c     * float(NGA3*100)/NSTEP,float(NGA4*100)/NSTEP
	 WRITE (IUTSP,'(I10,11F8.2," | ",5F8.2," | ",6F8.2)')
     *NSTEP,TG1,TG2,TG3,TG4,TG5,TG6,TG7,TG8,TG9,TG10,TGA,
     *TGA1,TGA2,TGA3,TGA4,TGA5,TG31,TG32,TG33,TG34,TG35,TG36
      ENDIF
      
      ENDIF
*******************************************************************************
*                                                                             *
*                          END NUMERICAL INTEGRATION                          *
*                                                                             *
******************************************************************************* 	
160   CONTINUE
170   CONTINUE  
      CLOSE (IUT90)
      CLOSE (IUT91)
      CLOSE (IUT93)
      CLOSE (IUT94)
      WRITE(*,*)'CONGRATULATION! SIMULATION IS FULFILLED SUCCESSFULLY!'
8300  FORMAT (/2X,'  ATTENZIONE :JOB SUCCESSFULLY COMPLETED; TIME = ',1
     *    P1E10.2,' DAYS',//)
      
      END PROGRAM UnFECOM
