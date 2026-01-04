#include "DEFS.h"
      
      MODULE MOD_GLOBAL
      IMPLICIT NONE
      SAVE        

!======================================================================
! PARAMS
	REAL, PARAMETER :: DAYI = 1. / 86400.
      REAL, PARAMETER :: GRAV = 9.806
      REAL, PARAMETER :: RHEAT=0.60
      REAL, PARAMETER :: ZETA1=0.35
      REAL, PARAMETER :: ZETA2=23.0
      REAL, PARAMETER :: ER = 6371000.
      
      INTEGER, PARAMETER :: IUNML = 5
      INTEGER, PARAMETER :: IUPRT = 7
      INTEGER, PARAMETER :: IUGRD = 8
      INTEGER, PARAMETER :: IUPLT = 10
      INTEGER, PARAMETER :: IUWRS = 13
      INTEGER, PARAMETER :: IURRS = 14
	INTEGER, PARAMETER :: IUIFO = 85
      INTEGER, PARAMETER :: IUMAT = 99669
      INTEGER, PARAMETER :: IUT90 = 90
      INTEGER, PARAMETER :: IUT91 = 91
      INTEGER, PARAMETER :: IUT92 = 92
      INTEGER, PARAMETER :: IUT93 = 93
      INTEGER, PARAMETER :: IUT94 = 94
      INTEGER, PARAMETER :: IUT95 = 95 
      INTEGER, PARAMETER :: IUT96 = 96
      INTEGER, PARAMETER :: IUT97 = 97
      INTEGER, PARAMETER :: IUT98 = 98
      INTEGER, PARAMETER :: IUT86 = 86
      INTEGER, PARAMETER :: IUT87 = 87
      INTEGER, PARAMETER :: IUT88 = 88
      INTEGER, PARAMETER :: IUT89 = 89
      
	INTEGER, PARAMETER :: IUT20 = 120   !SEDIMENT FLUX BOND
      INTEGER, PARAMETER :: IUT21 = 121   !SEDIMENT OBC  
	INTEGER, PARAMETER :: IUBCS = 220
	INTEGER, PARAMETER :: IUDEP = 221
	INTEGER, PARAMETER :: IUEBC = 222
	INTEGER, PARAMETER :: IUREBC = 223
	INTEGER, PARAMETER :: IUSBC = 224
	INTEGER, PARAMETER :: IUFBC = 225
	INTEGER, PARAMETER :: IUWDS = 226
	INTEGER, PARAMETER :: IUTSR = 227
	INTEGER, PARAMETER :: IUITS = 228
	INTEGER, PARAMETER :: IUGDH = 230
	INTEGER, PARAMETER :: IUSEC = 231
	INTEGER, PARAMETER :: IUUVK = 232
	INTEGER, PARAMETER :: IULAG = 233
	INTEGER, PARAMETER :: IUWEIR = 234
	INTEGER, PARAMETER :: IUDCHG = 235
      INTEGER, PARAMETER :: IULBA = 236
      INTEGER, PARAMETER :: IUTSP = 239 !DZL FOR TIMESTEP
	INTEGER, PARAMETER :: IUVER = 249  !XK
      !SED
	INTEGER, PARAMETER :: IUISED =254
	INTEGER, PARAMETER :: IUSDBC =255
      INTEGER, PARAMETER :: IUSDFLX=256
      INTEGER, PARAMETER :: IUD50 = 257
	INTEGER, PARAMETER :: IUTAUE = 258
	INTEGER, PARAMETER :: IUWAVE = 259
      !INFO EXCH
	INTEGER, PARAMETER :: IUNAG = 262
	INTEGER, PARAMETER :: IUNEG = 263
	INTEGER, PARAMETER :: IUNIG = 264
      INTEGER, PARAMETER :: IUTHR = 261   !OMP THREADS
	!CBCADJ
	INTEGER, PARAMETER :: IUCBCADJ = 91849      
	INTEGER, PARAMETER :: IUCBCADJS = 91894   
	INTEGER, PARAMETER :: IUPOI = 266
      INTEGER, PARAMETER :: IUFOBC = 267
      INTEGER, PARAMETER :: IUFTID = 268
      INTEGER, PARAMETER :: IUFSBC = 269
      
      INTEGER, PARAMETER :: IUAPR = 271
      INTEGER, PARAMETER :: IUATP = 272
      INTEGER, PARAMETER :: IURHM = 273
      INTEGER, PARAMETER :: IUCLD = 274
	
! END: PARAMS
!======================================================================
	
	
!======================================================================
! CORE VARS
	
	INTEGER KB
      INTEGER KBM1
      INTEGER KBM2
      INTEGER KSL
      INTEGER ISTART
      INTEGER IEND
	INTEGER N_DTI
	INTEGER IYEAR
      INTEGER IMONTH
      INTEGER IDAY0
	INTEGER N_RST
	INTEGER N_OPT
	INTEGER N_TSR
      INTEGER N_FPT
	INTEGER NSTEPS
      INTEGER NSTEP
	INTEGER NN_OPTRST
	INTEGER NN_OPT
      INTEGER NN_OPT2 !V2410
      INTEGER NN_FPT
	INTEGER NN_VIEW
      INTEGER IRAMP
	INTEGER N_SEC_XY
      INTEGER NTHREADS
      INTEGER IFCOUNT
      INTEGER COUNT_LOOP
      INTEGER SPECIES
      INTEGER NUMCBCADJ
      INTEGER ANUM
      INTEGER NPOI
      INTEGER NIT !V2410
      INTEGER, ALLOCATABLE :: NN_TSR(:)
	INTEGER, ALLOCATABLE :: NN_SEC(:)
	INTEGER, ALLOCATABLE :: NXY_SEC(:,:)
	
      REAL TG
	REAL DTI
	REAL HPRNU
      REAL VPRNU
      REAL UMOL
      REAL UMOL1
      REAL RAMP_TIDE
      REAL RAMP
      REAL HORCON
      REAL TLAG
      REAL TDAY
      REAL TDAY0 !FOR UV CHECK !DZL
      REAL THOUR
      REAL BFRIC
      REAL Z0B
      REAL Z00
      REAL CBCAA
      REAL CBCBB
      REAL CBCCC
      REAL CBCJG
      REAL AAVVEE
      REAL SSUUMM
      REAL CBCJG1
      REAL CBCJG2
      REAL DMIN
      REAL TIDE_LAG
      REAL FACTOR
	REAL S_BEG
	REAL BPG_BEG
	REAL FPTSTAR
      REAL FPTEND
	REAL DTIMAX
      REAL DTIMAX1
	REAL THOUR_HOT
	REAL T_END
	REAL ACCUM_SALT
	REAL RST_SPEC
	REAL RST_ARCH
      REAL T1F
      REAL T2F
      REAL T1FS
      REAL T2FS
c      REAL T1MAT !c V2410
c      REAL T2MAT !c V2410
      REAL EL_RISEUP
      
            

	REAL, ALLOCATABLE :: DZR(:)
      REAL, ALLOCATABLE :: Z(:)
      REAL, ALLOCATABLE :: ZZ(:)
      REAL, ALLOCATABLE :: DZ(:)
      REAL, ALLOCATABLE :: DZZ(:)
      REAL, ALLOCATABLE :: DPTHSL(:)
	REAL, ALLOCATABLE :: H(:)
      REAL, ALLOCATABLE :: H1(:)
      REAL, ALLOCATABLE :: H2(:)
      REAL, ALLOCATABLE :: H3(:)
      REAL, ALLOCATABLE :: DX1(:)
      REAL, ALLOCATABLE :: DX2(:)
      REAL, ALLOCATABLE :: DX3(:)
      REAL, ALLOCATABLE :: DX4(:)
      REAL, ALLOCATABLE :: DY1(:)
      REAL, ALLOCATABLE :: DY2(:)
      REAL, ALLOCATABLE :: DY3(:)
      REAL, ALLOCATABLE :: DY4(:)  
      REAL, ALLOCATABLE :: D(:)
      REAL, ALLOCATABLE :: ANG(:)
      REAL, ALLOCATABLE :: ARU(:)
      REAL, ALLOCATABLE :: ARV(:)
      REAL, ALLOCATABLE :: CBC_U(:)
      REAL, ALLOCATABLE :: CBC_V(:)
      REAL, ALLOCATABLE :: CBC_UV(:) 
      REAL, ALLOCATABLE :: CBC_COF(:)
      REAL, ALLOCATABLE :: CBC_COF_A1(:)
      REAL, ALLOCATABLE :: CBC_COF_B1(:)
      REAL, ALLOCATABLE :: CBC_COF_A2(:)
      REAL, ALLOCATABLE :: CBC_COF_B2(:)
      REAL, ALLOCATABLE :: DUM(:)
      REAL, ALLOCATABLE :: DVM(:)
      REAL, ALLOCATABLE :: FSM(:)
      REAL, ALLOCATABLE :: FSM11(:)
      REAL, ALLOCATABLE :: FSM12(:)
      REAL, ALLOCATABLE :: FSMADD(:)
      REAL, ALLOCATABLE :: FSMTEMP(:)
      REAL, ALLOCATABLE :: COR(:)
      REAL, ALLOCATABLE :: WTSURF(:)
      REAL, ALLOCATABLE :: WSSURF(:)
      REAL, ALLOCATABLE :: WUSURF(:)
      REAL, ALLOCATABLE :: WVSURF(:)
      REAL, ALLOCATABLE :: WUBOT(:)
      REAL, ALLOCATABLE :: WVBOT(:)
      REAL, ALLOCATABLE :: ELF(:)
      REAL, ALLOCATABLE :: EL(:)
      REAL, ALLOCATABLE :: EL1(:)
      REAL, ALLOCATABLE :: RR(:)
      REAL, ALLOCATABLE :: RRT(:)
      REAL, ALLOCATABLE :: AAUU(:)
      REAL, ALLOCATABLE :: AAVV(:)
      REAL, ALLOCATABLE :: AAZZ(:)
      REAL, ALLOCATABLE :: BBBB(:)
      REAL, ALLOCATABLE :: HU(:)
      REAL, ALLOCATABLE :: HV(:)
      REAL, ALLOCATABLE :: DU(:)
      REAL, ALLOCATABLE :: DV(:)
      REAL, ALLOCATABLE :: DJ(:)
      REAL, ALLOCATABLE :: XXI(:)
      REAL, ALLOCATABLE :: YXI(:)
      REAL, ALLOCATABLE :: XETA(:)
      REAL, ALLOCATABLE :: YETA(:)
      REAL, ALLOCATABLE :: SWRAD(:)
      REAL, ALLOCATABLE :: HTEMP(:)
      REAL, ALLOCATABLE :: XR(:)
      REAL, ALLOCATABLE :: YR(:)
      REAL, ALLOCATABLE :: XU(:)
      REAL, ALLOCATABLE :: YU(:)
      REAL, ALLOCATABLE :: XV(:)
      REAL, ALLOCATABLE :: YV(:)
      REAL, ALLOCATABLE :: LON(:)
      REAL, ALLOCATABLE :: LAT(:)
      REAL, ALLOCATABLE :: LONU(:)
      REAL, ALLOCATABLE :: LATU(:)
      REAL, ALLOCATABLE :: LONV(:)
      REAL, ALLOCATABLE :: LATV(:)
      REAL, ALLOCATABLE :: TPS(:)
      REAL, ALLOCATABLE :: P(:)
      REAL, ALLOCATABLE :: AP(:)
	REAL, ALLOCATABLE :: X_TSR(:)
      REAL, ALLOCATABLE :: Y_TSR(:)
	REAL, ALLOCATABLE :: X_SECS(:)
      REAL, ALLOCATABLE :: Y_SECS(:)
      REAL, ALLOCATABLE :: X_SECE(:)
      REAL, ALLOCATABLE :: Y_SECE(:)
      REAL, ALLOCATABLE :: DIS_SEC(:)
      REAL, ALLOCATABLE :: ANG_SEC(:)
	REAL, ALLOCATABLE :: FLUXSUM(:) 
      REAL, ALLOCATABLE :: SFLUXSUM(:)
      REAL, ALLOCATABLE :: SEDFLUXSUM(:)
	REAL, ALLOCATABLE :: KM_COF(:) !XK
      REAL, ALLOCATABLE :: KH_COF(:) !XK
	REAL, ALLOCATABLE :: WINDU(:)  
      REAL, ALLOCATABLE :: WINDV(:)
      REAL, ALLOCATABLE :: XD_SEC(:,:)
      REAL, ALLOCATABLE :: YD_SEC(:,:)
      REAL, ALLOCATABLE :: AAMEBC(:,:)
	REAL, ALLOCATABLE :: DTXX(:,:)
      REAL, ALLOCATABLE :: DTYY(:,:)
	REAL, ALLOCATABLE :: A(:,:)
      REAL, ALLOCATABLE :: C(:,:)
      REAL, ALLOCATABLE :: VH(:,:)
      REAL, ALLOCATABLE :: VHP(:,:)
      REAL, ALLOCATABLE :: PROD(:,:)
      REAL, ALLOCATABLE :: DTEF(:,:)
      REAL, ALLOCATABLE :: T(:,:)
      REAL, ALLOCATABLE :: S(:,:)
      REAL, ALLOCATABLE :: KM(:,:)
      REAL, ALLOCATABLE :: KH(:,:)
      REAL, ALLOCATABLE :: KQ(:,:)
      REAL, ALLOCATABLE :: L(:,:)
      REAL, ALLOCATABLE :: Q2(:,:)
      REAL, ALLOCATABLE :: Q2L(:,:)
      REAL, ALLOCATABLE :: AAM(:,:)
      REAL, ALLOCATABLE :: XMFLUX(:,:)
      REAL, ALLOCATABLE :: YMFLUX(:,:)
      REAL, ALLOCATABLE :: XMFLUX0(:,:) !FOR UV CHECK
      REAL, ALLOCATABLE :: YMFLUX0(:,:)
      REAL, ALLOCATABLE :: XFLUX(:,:)
      REAL, ALLOCATABLE :: YFLUX(:,:)
      REAL, ALLOCATABLE :: UF(:,:)
      REAL, ALLOCATABLE :: VF(:,:)
      REAL, ALLOCATABLE :: U(:,:)
      REAL, ALLOCATABLE :: V(:,:)
      REAL, ALLOCATABLE :: W(:,:)
      REAL, ALLOCATABLE :: UNN(:,:)
      REAL, ALLOCATABLE :: VNN(:,:)
      REAL, ALLOCATABLE :: WNN(:,:)
      REAL, ALLOCATABLE :: UR(:,:)
      REAL, ALLOCATABLE :: VR(:,:)
      REAL, ALLOCATABLE :: UC(:,:)
      REAL, ALLOCATABLE :: VC(:,:)
      REAL, ALLOCATABLE :: WR(:,:)
      REAL, ALLOCATABLE :: RHO(:,:)
      REAL, ALLOCATABLE :: DRHOX(:,:)
      REAL, ALLOCATABLE :: DRHOY(:,:)
      REAL, ALLOCATABLE :: CURV4(:,:)
      REAL, ALLOCATABLE :: T0ZZ(:,:)
      REAL, ALLOCATABLE :: S0ZZ(:,:)
	REAL, ALLOCATABLE :: DTT22(:,:)
      REAL, ALLOCATABLE :: DEI22(:,:)
      REAL, ALLOCATABLE :: GA(:)
      REAL, ALLOCATABLE :: HA(:)
      REAL, ALLOCATABLE :: CBCADJ(:)
      REAL, ALLOCATABLE :: XPOI(:)
c      REAL, ALLOCATABLE :: DPOI(:)  !CBR: scalarized
	REAL(8), ALLOCATABLE ::FLUX_ACM(:)
      REAL(8), ALLOCATABLE ::SFLUX_ACM(:)
      REAL(8), ALLOCATABLE ::SEDFLUX_ACM(:) !LXY
	
	
	CHARACTER*20 CASENAME
      CHARACTER*10 TOR
      CHARACTER*10 ADVECT
      CHARACTER*10 HORZMIX
      CHARACTER*10 VERTMIX
      CHARACTER*20 OPTEBC 
      CHARACTER*20 RESTAR
      CHARACTER*80 IN_DIRE
      CHARACTER*80 OUT_DIRE
      CHARACTER*20 OSTYPE
      CHARACTER*10 XG
      CHARACTER*20 WDSTYPE
	CHARACTER*20 Z0BTYP
	CHARACTER*20 DTITYPE
      CHARACTER*20 COORD
      CHARACTER*120 FN
	CHARACTER*20, ALLOCATABLE :: FN_TSR(:)
	
	LOGICAL LOG_MMTPT           
      LOGICAL LOG_CBCADJ
	LOGICAL LOG_BPG
      LOGICAL LOG_FSMPT
	LOGICAL LOG_SPT
	LOGICAL LOG_VPT
	LOGICAL LOG_EPT
	LOGICAL LOG_FPT
	LOGICAL LOG_DCHG
	LOGICAL LOG_BEDPT
	!LOGICAL LOG_RSTARC
	LOGICAL LOG_RSTSPEC
      LOGICAL HOT_TRY
	
! END: CORE VARS
!======================================================================
	
	
!======================================================================
! NEST
	INTEGER N_CTRD,N_CTRDP1
      INTEGER N_CTRD_AG,N_CTRD_VG,N_CTRD_EVG,N_CTRD_IVG
      INTEGER, ALLOCATABLE :: NIM1(:)
      INTEGER, ALLOCATABLE :: NIM2(:)
      INTEGER, ALLOCATABLE :: NIM3(:)
      INTEGER, ALLOCATABLE :: NIM4(:)
      INTEGER, ALLOCATABLE :: NIP1(:)
      INTEGER, ALLOCATABLE :: NIP2(:)
      INTEGER, ALLOCATABLE :: NIP3(:)
      INTEGER, ALLOCATABLE :: NIP4(:)
      INTEGER, ALLOCATABLE :: NJM1(:)
      INTEGER, ALLOCATABLE :: NJM2(:)
      INTEGER, ALLOCATABLE :: NJM3(:)
      INTEGER, ALLOCATABLE :: NJM4(:)
      INTEGER, ALLOCATABLE :: NJP1(:)
      INTEGER, ALLOCATABLE :: NJP2(:)
      INTEGER, ALLOCATABLE :: NJP3(:)
      INTEGER, ALLOCATABLE :: NJP4(:)
      INTEGER, ALLOCATABLE :: NIM1JP1(:)
      INTEGER, ALLOCATABLE :: NIM1JM1(:)
      INTEGER, ALLOCATABLE :: NIP1JP1(:)
      INTEGER, ALLOCATABLE :: NIP1JM1(:)
      INTEGER, ALLOCATABLE :: NIM1JM2(:)
      INTEGER, ALLOCATABLE :: NIM2JM1(:)
	
	! INFO EXCHANGE
	INTEGER, PARAMETER :: NMAX = 30
	INTEGER VTP
	!CHARACTER*10 UDSCH,INTSCH
      INTEGER, ALLOCATABLE :: NAGC_EVG(:,:)
      INTEGER, ALLOCATABLE :: NAGU_EVG(:,:)
      INTEGER, ALLOCATABLE :: NAGV_EVG(:,:)
	INTEGER, ALLOCATABLE :: NAGC_IVG(:,:)
      INTEGER, ALLOCATABLE :: NAGU_IVG(:,:)
      INTEGER, ALLOCATABLE :: NAGV_IVG(:,:)
	REAL, ALLOCATABLE :: ELM(:)
      REAL, ALLOCATABLE :: WTC_EVG(:,:)
      REAL, ALLOCATABLE :: WTU_EVG(:,:)
      REAL, ALLOCATABLE :: WTV_EVG(:,:)
	REAL, ALLOCATABLE :: WTC_IVG(:,:)
      REAL, ALLOCATABLE :: WTU_IVG(:,:)
      REAL, ALLOCATABLE :: WTV_IVG(:,:)
      
	! IDW
	INTEGER NPOW
      REAL SEARCHRADIUS
	! DR\UNI
      INTEGER NB
	REAL, ALLOCATABLE :: XNODE(:,:),YNODE(:,:)
      ! HPI
#ifdef INT_HPI
      REAL, ALLOCATABLE :: XIC(:)
      REAL, ALLOCATABLE :: XIU(:)
      REAL, ALLOCATABLE :: XIV(:)
#endif

#ifdef INT_HAEI
      REAL, ALLOCATABLE :: EPSC(:)
      REAL, ALLOCATABLE :: EPSU(:)
      REAL, ALLOCATABLE :: EPSV(:)
#endif
	
	!TSR
	INTEGER, ALLOCATABLE :: NAGC_TSR(:,:)
      INTEGER, ALLOCATABLE :: NAGU_TSR(:,:)
      INTEGER, ALLOCATABLE :: NAGV_TSR(:,:)
	REAL, ALLOCATABLE :: WTC_TSR(:,:)
      REAL, ALLOCATABLE :: WTU_TSR(:,:)
      REAL, ALLOCATABLE :: WTV_TSR(:,:)
      
      !Find interface AG
      integer NniAGcW, NniAGcE, NniAGcN, NniAGcS
      integer, allocatable :: niAGcW(:)
      integer, allocatable :: niAGcE(:)
      integer, allocatable :: niAGcN(:)
      integer, allocatable :: niAGcS(:)
      integer NniAGfW, NniAGfE, NniAGfN, NniAGfS
      integer, allocatable :: niAGfW(:)
      integer, allocatable :: niAGfE(:)
      integer, allocatable :: niAGfN(:)
      integer, allocatable :: niAGfS(:)
      
#ifdef NISpg
      integer N_NIS
      integer, allocatable :: nispg(:)
      integer, allocatable :: NAGC_NIS(:,:)
      integer, allocatable :: NAGU_NIS(:,:)
      integer, allocatable :: NAGV_NIS(:,:)
      real, allocatable :: WTC_NIS(:,:)
      real, allocatable :: WTU_NIS(:,:)
      real, allocatable :: WTV_NIS(:,:)
      real, allocatable :: mu(:)
#endif
	
! END: NEST
!======================================================================
	
	
!======================================================================
!SED

      REAL :: F_SD50
      REAL :: F_ALFA
      REAL :: F_MCSXS
      REAL :: SED_BEG
      REAL :: BED_BEG
      REAL :: T1SED
      REAL :: T2SED      
      REAL :: T1SEDQ
      REAL :: T2SEDQ         
	REAL, ALLOCATABLE :: SUMH(:) !LXY
      REAL, ALLOCATABLE :: TAUE0(:) !LXY
      REAL, ALLOCATABLE :: TAUD0(:) !LXY
	REAL, ALLOCATABLE :: D50(:)
      REAL, ALLOCATABLE :: TAU(:)
      REAL, ALLOCATABLE :: TAU_WAVE(:)
      REAL, ALLOCATABLE :: TAU_TIDE(:)
      REAL, ALLOCATABLE :: TAU_U(:)
      REAL, ALLOCATABLE :: TAU_V(:)
      REAL, ALLOCATABLE :: TAUE(:)
      REAL, ALLOCATABLE :: TAUD(:)
      REAL, ALLOCATABLE :: QERO(:)
      REAL, ALLOCATABLE :: QDEP(:)
      REAL, ALLOCATABLE :: QSI(:)
      REAL(8), ALLOCATABLE :: ZBED(:)
      REAL, ALLOCATABLE :: ZBEDD(:)
      REAL, ALLOCATABLE :: VSDDIST(:,:)
      REAL, ALLOCATABLE :: SEDDIS(:)
      REAL, ALLOCATABLE :: DSEDDIS(:,:)
      REAL, ALLOCATABLE :: SEDBDRYSL(:,:)
      REAL, ALLOCATABLE :: SEDBDRY(:,:)
      REAL, ALLOCATABLE :: DSEDBDRY(:,:,:)
	REAL, ALLOCATABLE :: SED(:,:)
      REAL, ALLOCATABLE :: WFF(:,:)
      REAL, ALLOCATABLE :: DWS(:,:)
      REAL, ALLOCATABLE :: VNIAN(:,:)
      REAL, ALLOCATABLE :: USIN(:)
	REAL, ALLOCATABLE :: FU1SED(:,:,:)
      REAL, ALLOCATABLE :: FU2SED(:,:,:)
      REAL, ALLOCATABLE :: FU3SED(:,:,:)
      REAL, ALLOCATABLE :: FU4SED(:,:,:)
      REAL, ALLOCATABLE :: FU5SED(:,:,:)
      REAL, ALLOCATABLE :: FV1SED(:,:,:)
      REAL, ALLOCATABLE :: FV2SED(:,:,:)
      REAL, ALLOCATABLE :: FV3SED(:,:,:)
      REAL, ALLOCATABLE :: FV4SED(:,:,:)
      REAL, ALLOCATABLE :: FV5SED(:,:,:)
	
	CHARACTER*20 TAUTYP
	LOGICAL :: LOG_ISED
	LOGICAL :: LOG_SEDPT
    
	
! END: SED
!======================================================================
	
	
!======================================================================
! FLUX ANALYSIS
	
	INTEGER JHM
	
      REAL, ALLOCATABLE :: AVGE(:)
      REAL, ALLOCATABLE :: DEI(:)
	REAL, ALLOCATABLE :: ARC_UA(:,:)
      REAL, ALLOCATABLE :: ARC_CA(:,:)
      REAL, ALLOCATABLE :: ARC_VA(:,:)
      REAL, ALLOCATABLE :: ARC_UHH(:,:)
      REAL, ALLOCATABLE :: ARC_CHH(:,:)
      REAL, ALLOCATABLE :: ARC_UCHH(:,:)
      REAL, ALLOCATABLE :: ARC_VHH(:,:)
      REAL, ALLOCATABLE :: ARC_VCHH(:,:)
      REAL, ALLOCATABLE :: ARC_UCDH(:,:)
      REAL, ALLOCATABLE :: ARC_VCDH(:,:)
      REAL, ALLOCATABLE :: ARC_H(:,:)
      REAL, ALLOCATABLE :: FU1_WATER(:,:)
      REAL, ALLOCATABLE :: FU2_WATER(:,:)
      REAL, ALLOCATABLE :: FU1(:,:)
      REAL, ALLOCATABLE :: FU2(:,:)
      REAL, ALLOCATABLE :: FU3(:,:)
      REAL, ALLOCATABLE :: FU4(:,:)
      REAL, ALLOCATABLE :: FU5(:,:)
      !REAL, ALLOCATABLE :: FU6(:,:)
      !REAL, ALLOCATABLE :: FU7(:,:)
      REAL, ALLOCATABLE :: FV1_WATER(:,:)
      REAL, ALLOCATABLE :: FV2_WATER(:,:)
      REAL, ALLOCATABLE :: FV1(:,:)
      REAL, ALLOCATABLE :: FV2(:,:)
      REAL, ALLOCATABLE :: FV3(:,:)
      REAL, ALLOCATABLE :: FV4(:,:)
      REAL, ALLOCATABLE :: FV5(:,:)
      !REAL, ALLOCATABLE :: FV6(:,:)
      !REAL, ALLOCATABLE :: FV7(:,:)
	REAL, ALLOCATABLE :: ARCET(:,:)
      REAL, ALLOCATABLE :: HIST(:,:)
	REAL, ALLOCATABLE :: ARCUX(:,:,:)
      REAL, ALLOCATABLE :: ARCVX(:,:,:)
      REAL, ALLOCATABLE :: ARCS(:,:,:)
      REAL, ALLOCATABLE :: ARCUR(:,:,:)
      REAL, ALLOCATABLE :: ARCVR(:,:,:)
      REAL, ALLOCATABLE :: ARCSUX(:,:,:)
      REAL, ALLOCATABLE :: ARCSVX(:,:,:)
	
	LOGICAL, ALLOCATABLE :: LOG_ARC(:)
	
! END: FLUX ANALYSIS
!======================================================================
	
	
!======================================================================
! BCOND
	
      INTEGER NUM_HARMCONST
      INTEGER NUMEBC
	INTEGER NUMQBC
      INTEGER NUMFOBC1
      INTEGER, ALLOCATABLE :: NETA(:)
      INTEGER, ALLOCATABLE :: NCON(:)
      INTEGER, ALLOCATABLE :: NAJA(:) !CBR: NETA - NCON - NAJA
      INTEGER, ALLOCATABLE :: NETA0(:)
      INTEGER, ALLOCATABLE :: NCON0(:)
	INTEGER, ALLOCATABLE :: NQD(:)
      INTEGER, ALLOCATABLE :: NQC(:)
      
	REAL T1E
      REAL T2E
	REAL T1Q
      REAL T2Q
	REAL T1M
      REAL T2M
	REAL T1TS
      REAL T2TS
	REAL T1EL
      REAL T2EL
      REAL QPREC
      REAL DQPREC(2)
      REAL QEVAP
      REAL DQEVAP(2)
      REAL TX
      REAL DTX(2)
      REAL TY
      REAL DTY(2)
      REAL HFLUX
      REAL DHFLUX(2)
	REAL, ALLOCATABLE :: EMEAN(:)
      REAL, ALLOCATABLE :: AMP(:,:)
      REAL, ALLOCATABLE :: PHASE(:,:)
      REAL, ALLOCATABLE :: TBDRYSL(:,:)
      REAL, ALLOCATABLE :: SBDRYSL(:,:)
	REAL, ALLOCATABLE :: VQDIST(:,:)
	REAL, ALLOCATABLE :: EBDRY(:)
      REAL, ALLOCATABLE :: DEBDRY(:,:)
	REAL, ALLOCATABLE :: QDIS(:)
      REAL, ALLOCATABLE :: DQDIS(:,:)
      REAL, ALLOCATABLE :: TDIS(:)
      REAL, ALLOCATABLE :: DTDIS(:,:)
      REAL, ALLOCATABLE :: SDIS(:)
      REAL, ALLOCATABLE :: DSDIS(:,:)
	REAL, ALLOCATABLE :: TBDRY(:,:)
      REAL, ALLOCATABLE :: DTBDRY(:,:,:)
      REAL, ALLOCATABLE :: SBDRY(:,:)
      REAL, ALLOCATABLE :: DSBDRY(:,:,:)
	REAL, ALLOCATABLE :: ELBC(:)
	REAL, ALLOCATABLE :: ELBC_TIDE(:) !CBR
      REAL, ALLOCATABLE :: DELBC(:,:)
      REAL, ALLOCATABLE :: VFOBC(:,:)
      REAL, ALLOCATABLE :: SFOBC(:,:)
      REAL, ALLOCATABLE :: TFOBC(:,:)
      REAL, ALLOCATABLE :: SDFOBC(:,:)
      REAL, ALLOCATABLE :: FTBDRY(:,:)
      REAL, ALLOCATABLE :: FSBDRY(:,:)
      REAL, ALLOCATABLE :: FSDBDRY(:,:)
      REAL, ALLOCATABLE :: DVFOBC(:,:,:)
      REAL, ALLOCATABLE :: FDTBDRY(:,:,:)
      REAL, ALLOCATABLE :: FDSBDRY(:,:,:)
	REAL, ALLOCATABLE :: FDSDBDRY(:,:,:)
! END: BCOND
!======================================================================
	
	
!======================================================================
! VIEW
      INTEGER LAYER
      REAL XMIN
      REAL YMIN
      REAL XMAX
      REAL YMAX
	CHARACTER*20 CHC_VIEW
!======================================================================
	
	
!======================================================================
! TMP
	
      REAL T_BEG
	REAL, ALLOCATABLE :: QS_RAD(:)
      !REAL, ALLOCATABLE :: QSMAX1(:,:)
      !REAL, ALLOCATABLE :: QSMAX2(:,:)
      REAL, ALLOCATABLE :: QB_RAD(:)
      REAL, ALLOCATABLE :: QE_RAD(:)
      REAL, ALLOCATABLE :: QH_RAD(:)
	
! END: TMP
!======================================================================
      REAL M_BEG
	
	
!======================================================================
! WEIR
	
	!LOGICAL LOG_WEIR
      INTEGER ND
      INTEGER, ALLOCATABLE :: WEIR_N(:)
      INTEGER, ALLOCATABLE :: WEIR_I(:,:)
      REAL, ALLOCATABLE :: WEIR_H(:,:)
      REAL, ALLOCATABLE :: DELTA0(:,:)
      REAL, ALLOCATABLE :: DELTA1(:,:)
      REAL, ALLOCATABLE :: OCSL(:,:)  !!
      INTEGER, ALLOCATABLE :: KB_WEIR_0(:,:)  !!
      INTEGER, ALLOCATABLE :: KB_WEIR_1(:,:)  !!
      INTEGER, ALLOCATABLE :: KB_WEIR_5(:,:)  !!
      CHARACTER*40, ALLOCATABLE :: WEIR_TYP(:) 
	
! END: WEIR
!======================================================================
	
	
!======================================================================
! WAVE
	REAL T1WAVE
      REAL T2WAVE
      
	REAL, ALLOCATABLE :: HSIG(:)
      REAL, ALLOCATABLE :: TSIG(:)
      REAL, ALLOCATABLE :: WAVEDIR(:)
      REAL, ALLOCATABLE :: UBM(:)
      REAL, ALLOCATABLE :: UBM_WB(:)

	REAL, ALLOCATABLE :: ACB_VARIAT(:)
      REAL, ALLOCATABLE :: VARIAT1(:)
      REAL, ALLOCATABLE :: UBAR(:,:)
      REAL, ALLOCATABLE :: VBAR(:,:)
      REAL, ALLOCATABLE :: DWHT(:,:)
      REAL, ALLOCATABLE :: DWPER(:,:)
      REAL, ALLOCATABLE :: DWDIR(:,:)
      REAL, ALLOCATABLE :: DWLEN(:,:)
      REAL, ALLOCATABLE :: WAVEFFP(:)
      
	REAL, ALLOCATABLE :: LEN_WAVE(:,:)
	LOGICAL :: LOG_WAVE
	
! END: WAVE
!======================================================================
      
!======================================================================
! EL
      INTEGER NNIVGCW
      INTEGER NNIVGCE
      INTEGER NNIVGCS
      INTEGER NNIVGCN
      INTEGER NNIVGFW
      INTEGER NNIVGFE
      INTEGER NNIVGFS
      INTEGER NNIVGFN
      INTEGER NNIVGFWALL
      INTEGER NNIVGFEALL
      INTEGER NNIVGFSALL
      INTEGER NNIVGFNALL
      INTEGER NNIVGFWS
      INTEGER NNIVGFES
      INTEGER NNIVGFWN
      INTEGER NNIVGFEN
      INTEGER NNVGF4AG
      INTEGER NNIVGW
      INTEGER NNIVGS
      
      INTEGER, ALLOCATABLE :: NIVGFW(:)
      INTEGER, ALLOCATABLE :: NIVGFE(:)
      INTEGER, ALLOCATABLE :: NIVGFN(:)
      INTEGER, ALLOCATABLE :: NIVGFS(:)
      INTEGER, ALLOCATABLE :: NIVGFWS(:)
      INTEGER, ALLOCATABLE :: NIVGFES(:)
      INTEGER, ALLOCATABLE :: NIVGFWN(:)
      INTEGER, ALLOCATABLE :: NIVGFEN (:)
      INTEGER, ALLOCATABLE :: NIVGFWALL(:,:)
      INTEGER, ALLOCATABLE :: NIVGFEALL(:,:)
      INTEGER, ALLOCATABLE :: NIVGFNALL(:,:)
      INTEGER, ALLOCATABLE :: NIVGFSALL(:,:)
      INTEGER, ALLOCATABLE :: NVGF4AG (:,:)
      INTEGER, ALLOCATABLE :: NAG4VGC (:,:)
      INTEGER, ALLOCATABLE :: NIVGW (:,:)
      INTEGER, ALLOCATABLE :: NIVGS (:,:)
  
      REAL, ALLOCATABLE :: RHSPOI(:)
      REAL, ALLOCATABLE :: APOI(:)
      INTEGER, ALLOCATABLE :: IAPOI(:)
      INTEGER, ALLOCATABLE :: JAPOI(:)
      INTEGER, ALLOCATABLE :: JAPOI2(:)
      INTEGER, ALLOCATABLE :: KAPOI(:) !CBR add
      INTEGER, ALLOCATABLE :: AIJ (:,:)
      
! END: EL
!======================================================================      

!====================================================================== !c V2410
c! MAT    
c      REAL, ALLOCATABLE :: MDIS (:)
c      REAL, ALLOCATABLE :: DMDIS(:,:)
c      REAL, ALLOCATABLE :: CONCENTRATION(:,:)
c
c! END: MAT
!====================================================================== !c V2410
 
!======================================================================
! METEOROLOGICAL    
      REAL T1A
      REAL T2A
      REAL T1R
      REAL T2R
      REAL T1C
      REAL T2C
      REAL T1APR
      REAL T2APR
      REAL, ALLOCATABLE :: ATP(:)
      REAL, ALLOCATABLE :: RHM(:)
      REAL, ALLOCATABLE :: CLOUD(:)
      REAL, ALLOCATABLE :: APR (:)
      REAL, ALLOCATABLE :: DATP (:,:)
      REAL, ALLOCATABLE :: DRHM (:,:)
      REAL, ALLOCATABLE :: DCLOUD (:,:)
      REAL, ALLOCATABLE :: DAPR (:,:)
      
! END: METEOROLOGICAL
!======================================================================      
      
      
!======================================================================
! SEC  
      INTEGER N_SEC_N
      INTEGER , ALLOCATABLE ::TOLNUM_SEC(:)
      INTEGER , ALLOCATABLE ::NUM_SEC(:,:)
      
! END: SEC
!======================================================================


!====================================================================== !V2410
! LAGRANGE TRACKING
      REAL EARLIEST_PARTICLE
      
      INTEGER NDT
      INTEGER I_LAG
      INTEGER N_LAG
      
      INTEGER, ALLOCATABLE :: PI_LAG(:)
      INTEGER, ALLOCATABLE :: PF_LAG(:)
      
      REAL, ALLOCATABLE :: PIX_LAG(:)
      REAL, ALLOCATABLE :: PIY_LAG(:)
      REAL, ALLOCATABLE :: PZ_LAG(:)
      
      REAL, ALLOCATABLE :: PB_LAG(:)
      REAL, ALLOCATABLE :: PS_LAG(:)
      REAL, ALLOCATABLE :: PT_LAG(:)
      REAL, ALLOCATABLE :: PENDX_LAG(:)
      REAL, ALLOCATABLE :: PENDY_LAG(:)
      
      
      
! END: LAGRANGE TRACKING
!======================================================================  !V2410     
!====================================================================== !V2410
! Tonghua
      LOGICAL LOG_TH
      LOGICAL OUT_T0
      LOGICAL OUT_T01
      LOGICAL OUT_T07
      INTEGER YEAR_TH
      INTEGER MONTH_TH
      INTEGER DAY_TH
      INTEGER HOUR_TH
      INTEGER MINUTE_TH
      INTEGER SECOND_TH
      
      REAL DELTA_TGTH_S
      REAL DELTA_TGTH_S1
      REAL DELTA_TGTH_S7

      REAL(8) TG_TH 
  
! END: Tonghua
!====================================================================== !V2410
      
!====================================================================== !V2410
! SLUICE
      LOGICAL, ALLOCATABLE :: LOG_SLUICE_ON(:)
      INTEGER, PARAMETER :: IUSLUICE = 275
      INTEGER N_SLUICE
      INTEGER, ALLOCATABLE :: SLUICE_NUM(:)
      INTEGER, ALLOCATABLE :: SLUICE_GRD(:,:)
      INTEGER, ALLOCATABLE :: SLUICE_FSM(:,:)
      CHARACTER*40, ALLOCATABLE :: SLUICE_TYP(:) 
! END: SLUICE
!====================================================================== !V2410
      
!====================================================================== !V2410
! MATERIAL      
      INTEGER NUMMAT
      
      REAL T1MAT
      REAL T2MAT
      
      INTEGER, ALLOCATABLE :: I_MAT(:)
      
      REAL, ALLOCATABLE :: XR_MAT(:)
      REAL, ALLOCATABLE :: YR_MAT(:)
      REAL, ALLOCATABLE :: VQDIST_MAT(:,:)
      REAL, ALLOCATABLE :: MDIS(:)
      REAL, ALLOCATABLE :: DMDIS(:,:)
      REAL, ALLOCATABLE :: CONCENTRATION(:,:)

      
! END: MATERIAL
!====================================================================== !V2410
      
      
! CBR :      
!======================================================================
! INFO_EXCH 
      REAL, ALLOCATABLE :: MASK4VAR(:,:)
! END: INFO_EXCH
!======================================================================
! PROFQ
      REAL, ALLOCATABLE :: GH(:,:),SM(:,:),SH(:,:)
      REAL, ALLOCATABLE :: KN(:,:),BOYGR(:,:)
! END: PROFQ
!======================================================================
! ADVT_HSIMT
      REAL, ALLOCATABLE :: KAX(:,:),KAY(:,:),KAZ(:,:)
! END: ADVT_HSIMT
!======================================================================
!======================================================================
! PROFT
      REAL, ALLOCATABLE :: RAD(:,:)
! END: PROFT
!======================================================================
      
      
      
      
      END MODULE MOD_GLOBAL