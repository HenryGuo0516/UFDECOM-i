#INCLUDE "DEFS.h"
	
	SUBROUTINE BCOND(IDX)
      
      USE MOD_GLOBAL
      
      INTEGER, INTENT(IN) :: IDX
      INTEGER FINDCTRD
      REAL FACT,FACT1,AMEAN,FUCCK,FORCE
      DIMENSION TXX(N_CTRD), TYY(N_CTRD)
      DATA PI2 /6.283185307/
      
#if defined LBCEL_RAD || defined LBCUV_RAD
      INTEGER NTEMP
      REAL(KIND=8), DIMENSION(N_CTRD_AG) :: GRAD
      REAL(KIND=8), PARAMETER :: EPS = 1.0E-20
      REAL(KIND=8) :: CE, CX, CFF, DUDE, DUDT, DUDX
#endif

      SELECT CASE (IDX)    
          

      CASE(1)
************************************************************************
*                              CASE (1)                                *
*                              U,V  OBC                                *          
************************************************************************
C   HERE ARE THE OPEN BOUNDARY CONDITIONS FOR U AND V. NOTE THAT
C   ALL THE BOUNDARIES ARE SPECIFIED AT THE LINES OF SEA ELEVATION.
C   IN THE C-GRID, THE ELEVATION AND V ARE ON THE SAME LINE ON THE
C   WEST AND EAST BOUNDARIES, WHILE THE ELEVATION AND U ARE ON THE
C   SAME LINE ON THE NORTH AND SOUTH BOUDARIES.
#ifdef TIDE_EL

#ifdef LBCUV_CHA

c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO N = 1, NUMEBC
      DO concurrent (N=1:NUMEBC) local(IE,IC,A0,A1,K,CP,CP1,FBF)
      IE = NETA(N)
      IC = NCON(N)
c      DO K = 1, KBM1
      IF (IE==NIP1(IC)) THEN  !EAST BOUNDARY
          A0 = 0.0
          A1 = 0.0
          DO K=1,KBM1
              A0 = A0+V(IE,K)
              A1 = A1+VF(IC,K)
          ENDDO
          A0 = A0/FLOAT(KBM1) !V(IE) AVERAGED ALONG DEPTH
          A1 = A1/FLOAT(KBM1) !VF(IC) AVERAGED ALONG DEPTH
          DO K = 1,KBM1
              VF(IC,K) = VF(IC,K)-A1
          ENDDO
          CP = SQRT(GRAV*H(IC))*DTI/H1(IC)
          CP1 = 1./(1.+CP)
          FBF = CP1*(A0+CP*A1)
          DO K = 1, KBM1
              VF(IE,K) = VF(NIM1(IE),K)+FBF
              VF(NIM1(IE),K) = VF(IC,K)+A1
              UF(NIP1(IE),K) = UF(IE,K)   
          ENDDO
          WVBOT(IE) = WVBOT(NIM1(IE))
      ELSEIF (IE==NIM1(IC)) THEN  !WEST BOUNDARY
          A0 = 0.0
          A1 = 0.0
          DO K = 1,KBM1
             A0 = A0+V(IE,K)
             A1 = A1+VF(IC,K)
          ENDDO
          A0 = A0/FLOAT(KBM1)
          A1 = A1/FLOAT(KBM1)
          DO K=1,KBM1
             VF(IC,K) = VF(IC,K)-A1
          ENDDO
c          IF(H(IC).GT.0) THEN
             CP = SQRT(GRAV*H(IC))*DTI/H1(IC)
             CP1 = 1./(1.+CP)
             FBF = CP1*(A0+CP*A1)
              DO K = 1, KBM1
                 VF(IE,K) = VF(NIP1(IE),K)+FBF
                 VF(IC,K) = VF(IC,K)+A1
                 UF(IE,K) = UF(NIP1(IE),K)   
              ENDDO
              WVBOT(IE) = WVBOT(NIP1(IE))
c         ENDIF
      ELSEIF (IE==NJP1(IC)) THEN  !NORTH BOUNDARY
          A0 = 0.0
          A1 = 0.0
          DO K = 1,KBM1
              A0 = A0+U(IE,K)
              A1 = A1+UF(IC,K)
          ENDDO
          A0 = A0/FLOAT(KBM1)
          A1 = A1/FLOAT(KBM1)
          DO K = 1,KBM1
             UF(IC,K) = UF(IC,K)-A1
          ENDDO
c          IF(H(IC).GT.0) THEN
              CP = SQRT(GRAV*H(IC))*DTI/H2(IC)
              CP1 = 1./(1.+CP)
              FBF = CP1*(A0+CP*A1)
              DO K = 1, KBM1
                  UF(IE,K) = UF(NJM1(IE),K)+FBF
                  UF(NJM1(IE),K) = UF(IC,K)+A1
                  VF(NJP1(IE),K) = VF(IE,K)   
              ENDDO
              WUBOT(IE) = WUBOT(NJM1(IE))
c          ENDIF
      ELSEIF (IE==NJM1(IC)) THEN  !SOUTH BOUNDARY
          A0 = 0.0
          A1 = 0.0
          DO K=1,KBM1
              A0 = A0+U(IE,K)
              A1 = A1+UF(IC,K)
          ENDDO
          A0 = A0/FLOAT(KBM1)
          A1 = A1/FLOAT(KBM1)
          DO K = 1,KBM1
              UF(IC,K) = UF(IC,K)-A1
          ENDDO
c          IF(H(IC).GT.0) THEN
              CP = SQRT(GRAV*H(IC))*DTI/H2(IC)
              CP1 = 1./(1.+CP)
              FBF = CP1*(A0+CP*A1)
              DO K = 1, KBM1
                  UF(IE,K) = UF(NJP1(IE),K)+FBF
                  UF(NJP1(IE),K) = UF(IC,K)+A1
                  VF(IE,K) = VF(NJP1(IE),K)   
              ENDDO
              WUBOT(IE) = WUBOT(NJP1(IE))
c          ENDIF
      ENDIF
c      ENDDO
      ENDDO 
c!$OMP END TARGET

#elif defined LBCUV_GRA

      DO N = 1, NUMEBC
c      DO concurrent (N=1:NUMEBC) local(IE,IC,K)
      IE = NETA(N)
      IC = NCON(N)
      DO K = 1, KBM1
      IF (IE==NIP1(IC)) THEN  !EAST BOUNDARY
          VF(IE,K) = VF(NIM1(IE),K)
      ELSEIF (IE==NIM1(IC)) THEN  !WEST BOUNDARY
          VF(IE,K) = VF(NIP1(IE),K)
      ELSEIF (IE==NJP1(IC)) THEN  !NORTH BOUNDARY
          UF(IE,K) = UF(NJM1(IE),K)
      ELSEIF (IE==NJM1(IC)) THEN  !SOUTH BOUNDARY
          UF(IE,K) = UF(NJP1(IE),K)
      ENDIF
      ENDDO
      ENDDO 

#elif defined LBCUV_RAD

      DO N = 1, NUMEBC
c      DO concurrent (N=1:NUMEBC) local(IE,IC,K,J,I,NTEMP,
c     * DDDT,DDDX,DDDE,CFF,CX,CE)
      IE = NETA(N)
      IC = NCON(N)
      DO K = 1, KBM1
      IF (IE==NIP1(IC)) THEN  !EAST BOUNDARY
          IF (DVM(IE)*DVM(NJP1(IE))*DVM(NJM1(IE)).EQ.0.) THEN
              VF(IE,K) = VF(NIM1(IE),K)
          ELSE
          DO J = 0,1
          DO I = -1,0
              NTEMP = FINDCTRD(IE,I,J)
              GRAD(NTEMP) = V(NTEMP,K)-V(NJM1(NTEMP),K)
          ENDDO
          ENDDO
          DDDT = V(NIM1(IE),K)-VF(NIM1(IE),K)
          DDDX = VF(NIM1(IE),K)-VF(NIM2(IE),K)
          IF ((DDDT*DDDX).LT.0.0) DDDT = 0.0
          IF ((DDDT*(GRAD(NIM1(IE))+GRAD(FINDCTRD(IE,-1,1)))).GT.0.0) 
     *    THEN
              DDDE = GRAD(NIM1(IE))
          ELSE 
              DDDE = GRAD(FINDCTRD(IE,-1,1))
          ENDIF
          CFF = MAX(DDDX*DDDX+DDDE*DDDE,EPS)
          CX = DDDT*DDDX
          CE = MIN(CFF,MAX(DDDT*DDDE,-CFF))
          VF(IE,K) = (CFF*V(IE,K)+CX*VF(NIM1(IE),K) 
     *    -MAX(CE,0.0)*GRAD(IE) -MIN(CE,0.0)*GRAD(NJP1(IE)))/(CFF+CX)
          ENDIF
      ELSEIF (IE==NIM1(IC)) THEN  !WEST BOUNDARY
          IF (DVM(IE)*DVM(NJP1(IE))*DVM(NJM1(IE)).EQ.0.) THEN
              VF(IE,K) = VF(NIP1(IE),K)
          ELSE
          DO J = 0,1
          DO I = 0,1
              NTEMP = FINDCTRD(IE,I,J)
              GRAD(NTEMP) = V(NTEMP,K)-V(NJM1(NTEMP),K)
          ENDDO
          ENDDO
          DDDT = V(NIP1(IE),K)-VF(NIP1(IE),K)
          DDDX = VF(NIP1(IE),K)-VF(NIP2(IE),K)
          IF ((DDDT*DDDX).LT.0.0) DDDT = 0.0
          IF ((DDDT*(GRAD(NIP1(IE))+GRAD(FINDCTRD(IE,1,1)))).GT.0.0) 
     * THEN
              DDDE = GRAD(NIP1(IE))
          ELSE 
              DDDE = GRAD(FINDCTRD(IE,1,1))
          ENDIF
          CFF = MAX(DDDX*DDDX+DDDE*DDDE,EPS)
          CX = DDDT*DDDX
          CE = MIN(CFF,MAX(DDDT*DDDE,-CFF))
          VF(IE,K) = (CFF*V(IE,K)+CX*VF(NIP1(IE),K) 
     *    -MAX(CE,0.0)*GRAD(IE)-MIN(CE,0.0)*GRAD(NJP1(IE)))/(CFF+CX)
          ENDIF 
      ELSEIF (IE==NJP1(IC)) THEN  !NORTH BOUNDARY
          IF (DUM(NIM1(IE))*DUM(IE)*DUM(NIP1(IE)).EQ.0.) THEN
              UF(IE,K) = UF(NJM1(IE),K)
          ELSE
          DO J = -1,0
          DO I = -1,0
              NTEMP = FINDCTRD(IE,I,J)
              GRAD(NTEMP) = U(NIP1(NTEMP),K)-U(NTEMP,K)
          ENDDO
          ENDDO
          DDDT = U(NJM1(IE),K)-UF(NJM1(IE),K)
          DDDE = UF(NJM1(IE),K)-UF(NJM2(IE),K)
          IF ((DDDT*DDDE).LT.0.0) DDDT = 0.0
          IF ((DDDT*(GRAD(FINDCTRD(IE,-1,-1))+GRAD(NJM1(IE)))).GT.0.0)
     * THEN
              DDDX = GRAD(FINDCTRD(IE,-1,-1))
          ELSE
              DDDX = GRAD(NJM1(IE))
          ENDIF
          CFF = MAX(DDDX*DDDX+DDDE*DDDE,EPS)
          CX = MIN(CFF,MAX(DDDT*DDDX,-CFF))
          CE = DDDT*DDDE
          UF(IE,K) = (CFF*U(IE,K)+CE*UF(NJM1(IE),K) 
     *    -MAX(CX,0.0)*GRAD(NIM1(IE))-MIN(CX,0.0)*GRAD(IE))/(CFF+CE)
          ENDIF
      ELSEIF (IE==NJM1(IC)) THEN  !SOUTH BOUNDARY
          IF (DUM(NIM1(IE))*DUM(IE)*DUM(NIP1(IE)).EQ.0.) THEN
              UF(IE,K) = UF(NJP1(IE),K)
          ELSE
          DO J = 0,1
          DO I = -1,0
              NTEMP = FINDCTRD(IE,I,J)
              GRAD(NTEMP) = U(NIP1(NTEMP),K)-U(NTEMP,K)
          ENDDO
          ENDDO
          DDDT = U(NJP1(IE),K)-UF(NJP1(IE),K)
          DDDE = UF(NJP1(IE),K)-UF(NJP2(IE),K)
          IF ((DDDT*DDDE).LT.0.0) DDDT = 0.0
          IF ((DDDT*(GRAD(FINDCTRD(IE,-1,1))+GRAD(NJP1(IE)))).GT.0.0) 
     * THEN
          DDDX = GRAD(FINDCTRD(IE,-1,1))
          ELSE 
              DDDX = GRAD(NJP1(IE))
          ENDIF
          CFF = MAX(DDDX*DDDX+DDDE*DDDE,EPS)
          CX = MIN(CFF,MAX(DDDT*DDDX,-CFF))
          CE = DDDT*DDDE
          UF(IE,K) = (CFF*U(IE,K)+CE*UF(NJP1(IE),K) 
     *    -MAX(CX,0.0)*GRAD(NIM1(IE))-MIN(CX,0.0)*GRAD(IE))/(CFF+CE)
        ENDIF             
      ENDIF
      ENDDO
      ENDDO 
#endif

#endif

!======================================================================
!                         END: U,V OBC
!======================================================================  
      
     
      CASE(2)
************************************************************************
*                              CASE (2)                                *
*                               W OBC                                  *
************************************************************************
#ifdef TIDE_EL
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO N = 1,NUMEBC
c      DO K = 1,KBM1
      DO concurrent (N=1:NUMEBC,K=1:KBM1) local(IE,IC)
      IE = NETA(N)
      IC = NCON(N)
          IF (FSM(NIP1(IE)).EQ.0.0.AND.
     *    (IE==NIM1(IC).OR.IE==NIP1(IC))) THEN
              W(IE,K) = W(NIM1(IE),K)
          ELSEIF (FSM(NIM1(IE)).EQ.0.0.AND.
     *    (IE==NIM1(IC).OR.IE==NIP1(IC))) THEN
              W(IE,K) = W(NIP1(IE),K)
          ELSEIF (FSM(NJP1(IE)).EQ.0.0.AND.
     *    (IE==NJP1(IC).OR.IE==NJM1(IC))) THEN
              W(IE,K) = W(NJM1(IE),K)
          ELSEIF (FSM(NJM1(IE)).EQ.0.0.AND.
     *    (IE==NJP1(IC).OR.IE==NJM1(IC))) THEN
              W(IE,K) = W(NJP1(IE),K)
          ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET
#endif
!======================================================================   
!                             END: W OBC 
!======================================================================
      
      
      CASE(3)
************************************************************************
*                              CASE (3)                                *
*                            ELEVATION OBC                             *
************************************************************************
      IF (TDAY.GE.TIDE_LAG) THEN 
          IDAY = IDAY0+INT(TDAY)
	    TIME0 = TDAY*24.*3600.
	    IF(TIME0.GE.86400.) TIME0=(TDAY-INT(TDAY))*86400.0 
c	    DO I = 1,NUM_HARMCONST
c	      HA(I) = AMP(N,I)
c	      GA(I) = PHASE(N,I)
c          ENDDO
	    AMEAN = 0.0
	    IF (NUM_HARMCONST.EQ.16) THEN
	        CALL GETEL16EBC(IYEAR,IMONTH,IDAY,TIME0,AMEAN,
     *NUMEBC,AMP,PHASE,ELBC_TIDE)
	    ELSEIF(NUM_HARMCONST.EQ.17) THEN
	        CALL GETEL17EBC(IYEAR,IMONTH,IDAY,TIME0,AMEAN,
     *NUMEBC,AMP,PHASE,ELBC_TIDE)
          ENDIF	
!$ACC UPDATE DEVICE (ELBC_TIDE)
      ENDIF
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO N = 1, NUMEBC
      DO concurrent (N=1:NUMEBC) local(II,II1,K,FUCCK,AMEAN,FRESH)
      II = NETA(N)
      II1= NCON(N)
          
****************************************************************       
#ifdef TIDE_EL   
****************************************************************  
#ifdef ELB
      ELBC_TIDE(N) = ELBC_TIDE(N)*0.01+ELBC(N)*RAMP
#else
      ELBC_TIDE(N) = ELBC_TIDE(N)*0.01
#endif
      ELF(II) = ELBC_TIDE(N)+EL_RISEUP
      EL(II) = ELBC_TIDE(N)+EL_RISEUP
****************************************************************       
#elif defined TIDE_FLATHER   
****************************************************************
#ifdef ELB
      AMEAN=H(II)+ELBC(N)*RAMP+EL_RISEUP
#else      
      AMEAN=H(II)+EL_RISEUP
#endif      
      DO K=1,KBM1
          FUCCK=VFOBC(N,K)+ELBC_TIDE(N)*0.01
          IF (II.EQ.NJP1(II1)) THEN
              VF(NJP1(II),K)=FUCCK
              YMFLUX(NJP1(II),K)=FUCCK*AMEAN*H1(II)
          ELSEIF(II.EQ.NJM1(II1)) THEN
              VF(II,K)=FUCCK
              YMFLUX(II,K)=FUCCK*AMEAN*H1(II)
          ELSEIF(II.EQ.NIM1(II1)) THEN
              UF(II,K)=FUCCK
              XMFLUX(II,K)=FUCCK*AMEAN*H2(II)
          ELSEIF(II.EQ.NIP1(II1)) THEN
              UF(NIP1(II),K)=FUCCK
              XMFLUX(NIP1(II),K)=FUCCK*AMEAN*H2(II)
          ENDIF
      ENDDO
****************************************************************  
#elif defined TIDE_FLUX
****************************************************************
      DO K=1,KBM1
          FRESH=0.
#ifdef shelf_circulation
      IF (TDAY.GE.TIDE_LAG) THEN
          FRESH=VFOBC(N,K)+ELBC_TIDE(N)*0.01 !V2412b
      ENDIF
#else
      IF (TDAY.GE.TIDE_LAG) THEN
          FRESH=ELBC_TIDE(N)*0.01 !V2412b
      ENDIF
#endif
      IF (II.EQ.NJP1(II1)) THEN
          YMFLUX(II,K)=FRESH*H(II1)*H1(II1)
      ELSEIF(II.EQ.NJM1(II1)) THEN
          YMFLUX(II1,K)=FRESH*H(II1)*H1(II1)
      ELSEIF(II.EQ.NIM1(II1)) THEN
          XMFLUX(II1,K)=FRESH*H(II1)*H2(II1)
      ELSEIF(II.EQ.NIP1(II1)) THEN
          XMFLUX(II,K)=FRESH*H(II1)*H2(II1)
      ENDIF
      ENDDO
****************************************************************
#endif
      ENDDO
c!$OMP END TARGET      
!======================================================================   
!                        END: ELEVATION OBC
!======================================================================
      
************************************************************************
*                              CASE (4)                                *
*                              FLUX BC                                 *
************************************************************************
      CASE(4)
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO N = 1, NUMQBC
c      DO K = 1, KBM1
      DO concurrent (N=1:NUMQBC,K=1:KBM1) local(ID,IC,FRESH)
      ID = NQD(N)
      IC = NQC(N) 
      FRESH = QDIS(N)*RAMP*VQDIST(N,K)
      IF (ID==NIM1(IC) .OR. ID==NIP1(IC)) THEN !U DIRECTION
          IF (ID==NIM1(IC)) THEN  !NEGATIVE
              XMFLUX(IC,K) = FRESH / DZ(K)
          ELSE                      !POSITIVE
              XMFLUX(ID,K) = FRESH / DZ(K)
              !UF(ID,K)=XMFLUX(ID,K)/D(ID)/H2(ID) !MR
          ENDIF
      ELSE  !V DIRECTION 
!      ELSEIF (ID==NJP1(IC) .OR. ID==NJM1(IC)) THEN !V DIRECTION 
          IF (ID==NJM1(IC)) THEN  !NEGATIVE 
              YMFLUX(IC,K) = FRESH / DZ(K)
          ELSE                      !POSITIVE
              YMFLUX(ID,K) = FRESH / DZ(K)
          ENDIF
!      ELSE
!          WRITE(*,*) 'ERROR IN FLUXBOND!'
!          PAUSE
!          STOP
      ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET  
!======================================================================   
!                          END: FLUX BC
!======================================================================

      
      CASE(5)
!======================================================================   
!                         KM,KH,Q2,Q2L,L OBC
!======================================================================
  !180 CONTINUE
C------------------------- Q2 AND Q2L B.C.'S --------------------------
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO N = 1,NUMQBC
c        DO K = 1,KBM1
        DO concurrent (N=1:NUMQBC,K=1:KBM1) local(ID,IC)
         ID = NQD(N)
         IC = NQC(N)
          UF(IC,K) = UF(ID,K)
          VF(IC,K) = VF(ID,K)
          L(IC,K) = L(ID,K)
          KM(IC,K) = KM(ID,K)
          KH(IC,K) = KH(ID,K)
          KQ(IC,K) = KQ(ID,K)
        ENDDO
c      ENDDO
c!$OMP END TARGET 

      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO N = 1,NUMEBC
c        DO K = 1,KBM1
        DO concurrent (N=1:NUMEBC,K=1:KBM1) local(IE,IC)
        IE = NETA(N)
        IC = NCON(N)
#if defined TIDE_FLUX  || defined TIDE_FLATHER
          UF(IE,K)=UF(IC,K)
          VF(IE,K)=VF(IC,K)
          KM(IE,K)=KM(IC,K)
          KH(IE,K)=KH(IC,K)
          KQ(IE,K)=KQ(IC,K)
          L (IE,K)=L (IC,K)
#elif defined TIDE_EL
          UF(IE,K)=UF(IC,K)
          VF(IE,K)=VF(IC,K)
          KM(IE,K)=KM(IC,K)
          KH(IE,K)=KH(IC,K)
          KQ(IE,K)=KQ(IC,K)
          L (IE,K)=L (IC,K)
     !!     IF(FSM(NIP1(IE)).EQ.0.0.AND.
     !!*    (IE==NIM1(IC).OR.IE==NIP1(IC))) THEN
     !!       UF(IE,K) = UF(NIM1(IE),K)
     !!       VF(IE,K) = VF(NIM1(IE),K)
     !!       L(IE,K) = L(NIM1(IE),K)
     !!       KM(IE,K) = KM(NIM1(IE),K)
     !!       KH(IE,K) = KH(NIM1(IE),K)
     !!       KQ(IE,K) = KQ(NIM1(IE),K)
     !!     ELSEIF(FSM(NIM1(IE)).EQ.0.0.AND.
     !!*          (IE==NIM1(IC).OR.IE==NIP1(IC))) THEN
     !!       UF(IE,K) = UF(NIP1(IE),K)
     !!       VF(IE,K) = VF(NIP1(IE),K)
     !!       L(IE,K) = L(NIP1(IE),K)
     !!       KM(IE,K) = KM(NIP1(IE),K)
     !!       KH(IE,K) = KH(NIP1(IE),K)
     !!       KQ(IE,K) = KQ(NIP1(IE),K)
     !!     ELSEIF(FSM(NJP1(IE)).EQ.0.0.AND.
     !!*          (IE==NJP1(IC).OR.IE==NJM1(IC))) THEN
     !!       UF(IE,K) = UF(NJM1(IE),K)
     !!       VF(IE,K) = VF(NJM1(IE),K)
     !!       L(IE,K) = L(NJM1(IE),K)
     !!       KM(IE,K) = KM(NJM1(IE),K)
     !!       KH(IE,K) = KH(NJM1(IE),K)
     !!       KQ(IE,K) = KQ(NJM1(IE),K)
     !!     ELSEIF(FSM(NJM1(IE)).EQ.0.0.AND.
     !!*          (IE==NJP1(IC).OR.IE==NJM1(IC))) THEN
     !!       UF(IE,K) = UF(NJP1(IE),K)
     !!       VF(IE,K) = VF(NJP1(IE),K)
     !!       L(IE,K) = L(NJP1(IE),K)
     !!       KM(IE,K) = KM(NJP1(IE),K)
     !!       KH(IE,K) = KH(NJP1(IE),K)
     !!       KQ(IE,K) = KQ(NJP1(IE),K)
     !!     ENDIF
#endif
      ENDDO
c      ENDDO
c!$OMP END TARGET 
!======================================================================   
!                     END: KM,KH,Q2,Q2L,L OBC
!======================================================================
      
      
      CASE(6)
!======================================================================   
!                           TEMP & SAL OBC 
!======================================================================
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO N = 1,NUMQBC
c        DO K = 1,KBM1
        DO concurrent (N=1:NUMQBC,K=1:KBM1) local(ID,IC)
        ID = NQD(N)
        IC = NQC(N)
          IF (VQDIST(N,K).NE.0.0.AND.QDIS(N).GE.0.0) THEN  
            IF (FSM(NIP1(ID)).EQ.0.0.AND.
     *         (ID==NIM1(IC).OR.ID==NIP1(IC))) THEN  ! EAST
              UF(IC,K) = UF(ID,K) !TEMPERATURE
              VF(IC,K) = VF(ID,K) !SALINITY
            ELSEIF (FSM(NIM1(ID)).EQ.0.0.AND.
     *              (ID==NIM1(IC).OR.ID==NIP1(IC))) THEN  !WEST
              UF(IC,K) = TDIS(N)
              VF(IC,K) = SDIS(N)
            ELSEIF (FSM(NJP1(ID)).EQ.0.0.AND.
     *              (ID==NJP1(IC).OR.ID==NJM1(IC))) THEN  !NORTH
              UF(IC,K) = UF(ID,K)
              VF(IC,K) = VF(ID,K)
            ELSEIF (FSM(NJM1(ID)).EQ.0.0.AND.
     *              (ID==NJP1(IC).OR.ID==NJM1(IC))) THEN  !SOUTH
            UF(IC,K) = TDIS(N)
            VF(IC,K) = SDIS(N)
            ENDIF
          ELSEIF (VQDIST(N,K).NE.0.0.AND.QDIS(N).LT.0.0) THEN
            IF (FSM(NIP1(ID)).EQ.0.0.AND.
     *         (ID==NIM1(IC).OR.ID==NIP1(IC))) THEN  ! EAST
              UF(IC,K) = TDIS(N)
              VF(IC,K) = SDIS(N)
            ELSEIF (FSM(NIM1(ID)).EQ.0.0.AND.
     *              (ID==NIM1(IC).OR.ID==NIP1(IC))) THEN  !WEST
              UF(IC,K) = UF(ID,K)
              VF(IC,K) = VF(ID,K)            
            ELSEIF (FSM(NJP1(ID)).EQ.0.0.AND.
     *              (ID==NJP1(IC).OR.ID==NJM1(IC))) THEN  !NORTH
              UF(IC,K) = TDIS(N)
              VF(IC,K) = SDIS(N)            
            ELSEIF (FSM(NJM1(ID)).EQ.0.0.AND.
     *              (ID==NJP1(IC).OR.ID==NJM1(IC))) THEN  !SOUTH
              UF(IC,K) = UF(ID,K)
              VF(IC,K) = VF(ID,K)             
            ENDIF
	    ENDIF
        ENDDO
c      ENDDO
c!$OMP END TARGET 
C----------------------------------------------------------------------
#if defined TIDE_FLUX  || defined TIDE_FLATHER
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO N = 1, NUMEBC   
c      DO K = 1, KBM1 
      DO concurrent (N=1:NUMEBC,K=1:KBM1) local(IC)
          IC = NETA(N)
c          ID = NCON(N)
          UF(IC,K) = TBDRY(N,K)
          VF(IC,K) = SBDRY(N,K)
      ENDDO
c      ENDDO
c!$OMP END TARGET 
#elif defined TIDE_EL

c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)

c      DO N = 1, NUMEBC   
c      DO  K = 1, KBM1 
c      IC = NETA(N)
c      ID = NCON(N)
c          UF(IC,K) = TBDRY(N,K)
c          VF(IC,K) = SBDRY(N,K)
c      ENDDO
c      ENDDO

c      DO N = 1,NUMEBC
c      DO  K = 1, KBM1 
      DO concurrent (N=1:NUMEBC,K=1:KBM1) local(IE,IC,TBDY,SBDY,
     * VEL,CPH,CPH1)
        IE = NETA(N)
        IC = NCON(N)
          TBDY = TBDRY(N,K)
          SBDY = SBDRY(N,K)
     
          IF (FSM(NIP1(IE)).EQ.0.0.AND.
     *       (IE==NIM1(IC).OR.IE==NIP1(IC))) THEN  !EAST
            VEL = U(IE,K)
            IF (VEL.LE.0.0) THEN
              UF(IE,K) = TBDY  
              VF(IE,K) = SBDY 
            ELSE
              CPH = ABS(VEL)*DTI*2.0/(H1(IE)+H1(NIM1(IE)))
              CPH1=1./(1.+CPH)
              UF(IE,K) = CPH1*(T(IE,K)+CPH*UF(IC,K))
              VF(IE,K) = CPH1*(S(IE,K)+CPH*VF(IC,K))
            ENDIF
     
          ELSEIF (FSM(NIM1(IE)).EQ.0.0.AND.
     *            (IE==NIM1(IC).OR.IE==NIP1(IC))) THEN !WEST
            VEL = U(NIP1(IE),K)
            IF (VEL.GE.0.0) THEN
              UF(IE,K) = TBDY  
              VF(IE,K) = SBDY
            ELSE
              CPH = ABS(VEL)*DTI*2.0/(H1(IE)+H1(NIP1(IE)))
              CPH1=1./(1.+CPH)
              UF(IE,K) = CPH1*(T(IE,K)+CPH*UF(IC,K))
              VF(IE,K) = CPH1*(S(IE,K)+CPH*VF(IC,K)) 
            ENDIF
     
          ELSEIF (FSM(NJP1(IE)).EQ.0.0.AND.
     *            (IE==NJP1(IC).OR.IE==NJM1(IC))) THEN  !NORTH
            VEL = V(IE,K)
            IF (VEL.LE.0.0) THEN
              UF(IE,K) = TBDY
              VF(IE,K) = SBDY 
            ELSE
              CPH = ABS(VEL)*DTI*2.0/(H2(IE)+H2(NJM1(IE)))
              CPH1=1./(1.+CPH)
             UF(IE,K)=CPH1*(T(IE,K)+CPH*UF(NJM1(IE),K))
             VF(IE,K)=CPH1*(S(IE,K)+CPH*VF(NJM1(IE),K))
            ENDIF
     
          ELSEIF (FSM(NJM1(IE)).EQ.0.0.AND.
     *            (IE==NJP1(IC).OR.IE==NJM1(IC))) THEN  !SOUTH
            VEL = V(NJP1(IE),K)
            IF (VEL.GE.0.0) THEN
              UF(IE,K) = TBDY 
              VF(IE,K) = SBDY 
            ELSE
              CPH = ABS(VEL)*DTI*2.0/(H2(IE)+H2(NJP1(IE)))
              CPH1=1./(1.+CPH)
              UF(IE,K)=CPH1*(T(IE,K)+CPH*UF(NJP1(IE),K))
              VF(IE,K)=CPH1*(S(IE,K)+CPH*VF(NJP1(IE),K))
            ENDIF
     
          ENDIF
        ENDDO
c      ENDDO 
c!$OMP END TARGET 
#endif
      RETURN
!======================================================================   
!                          END: TEMP & SAL OBC 
!======================================================================
      
      
      CASE(7)
!======================================================================   
!                            SEDIMENT OBC
!======================================================================

	DO N = 1,NUMQBC    
        ID = NQD(N)
        IC = NQC(N)
        DO K = 1,KBM1
          IF (VQDIST(N,K).NE.0.0 .AND. QDIS(N).GE.0.0) THEN  
            IF (FSM(NIP1(ID)).EQ.0.0 .AND. 
     * (ID==NIM1(IC) .OR. ID==NIP1(IC))) THEN  ! EAST
              WFF(IC,K) = WFF(ID,K)
            ELSEIF (FSM(NIM1(ID)).EQ.0.0 .AND. 
     * (ID==NIM1(IC) .OR. ID==NIP1(IC))) THEN  !WEST
              WFF(IC,K) = SEDDIS(N)
            ELSEIF (FSM(NJP1(ID)).EQ.0.0 .AND. 
     * (ID==NJP1(IC) .OR. ID==NJM1(IC))) THEN  !NORTH
              WFF(IC,K) = WFF(ID,K)
            ELSEIF (FSM(NJM1(ID)).EQ.0.0 .AND. 
     * (ID==NJP1(IC) .OR. ID==NJM1(IC))) THEN  !SOUTH
              WFF(IC,K) = SEDDIS(N)
            ENDIF
          ELSEIF (VQDIST(N,K).NE.0.0 .AND. QDIS(N).LT.0.0) THEN
            IF (FSM(NIP1(ID)).EQ.0.0 .AND. 
     * (ID==NIM1(IC) .OR. ID==NIP1(IC))) THEN  ! EAST
              WFF(IC,K) = SEDDIS(N)
            ELSEIF (FSM(NIM1(ID)).EQ.0.0 .AND. 
     * (ID==NIM1(IC) .OR. ID==NIP1(IC))) THEN  !WEST
              WFF(IC,K) = WFF(ID,K)
            ELSEIF (FSM(NJP1(ID)).EQ.0.0 .AND. 
     * (ID==NJP1(IC) .OR. ID==NJM1(IC))) THEN  !NORTH
              WFF(IC,K) = SEDDIS(N)
            ELSEIF (FSM(NJM1(ID)).EQ.0.0 .AND. 
     * (ID==NJP1(IC) .OR. ID==NJM1(IC))) THEN  !SOUTH
              WFF(IC,K) = WFF(ID,K)
            ENDIF
	    ENDIF
        ENDDO
      ENDDO

	DO 743 N = 1,NUMEBC   
        IE = NETA(N)
        IC = NCON(N)
        DO 733 K = 1,KBM1
          SEDBDY = SEDBDRY(N,K)

		IF (FSM(NIP1(IE)).EQ.0.0.AND.
     *    (IE==NIM1(IC).OR.IE==NIP1(IC))) THEN  !EAST
            VEL = U(IE,K)
            IF (VEL.LE.0.0) THEN
              WFF(IE,K) = SEDBDY  
            ELSE
              CPH = ABS(VEL)*DTI*2.0/(H1(IE)+H1(NIM1(IE)))
              CPH1 = 1./(1+CPH)
			WFF(IE,K) = CPH1*(SED(IE,K)+CPH*WFF(IC,K))
            ENDIF
            GO TO 733
		  
	    ELSEIF (FSM(NIM1(IE)).EQ.0.0 .AND.
     *    (IE==NIM1(IC).OR.IE==NIP1(IC))) THEN  !WEST
            VEL = U(NIP1(IE),K)
            IF (VEL.GE.0.0) THEN
              WFF(IE,K) = SEDBDY
            ELSE
              CPH = ABS(VEL)*DTI*2.0/(H1(IE)+H1(NIP1(IE)))
              CPH1 = 1./(1+CPH)
              WFF(IE,K) = CPH1*(SED(IE,K)+CPH*WFF(IC,K)) 
            ENDIF
            GO TO 733

		ELSEIF (FSM(NJP1(IE)).EQ.0.0.AND.
     *    (IE==NJP1(IC).OR.IE==NJM1(IC))) THEN  !NORTH
            VEL = V(IE,K)
            IF (VEL.LE.0.0) THEN
              WFF(IE,K) = SEDBDY 
		  ELSE
              CPH = ABS(VEL)*DTI*2.0/(H2(IE)+H2(NJM1(IE)))
              CPH1 = 1./(1+CPH)
              WFF(IE,K) = CPH1*(SED(IE,K)+CPH*WFF(NJM1(IE),K))
            ENDIF
		  GO TO 733
		 
		ELSEIF (FSM(NJM1(IE)).EQ.0.0.AND.
     *    (IE==NJP1(IC).OR.IE==NJM1(IC))) THEN !SOUTH  
            VEL = V(NJP1(IE),K)
            IF (VEL.GE.0.0) THEN
              WFF(IE,K) = SEDBDY 
            ELSE
              CPH = ABS(VEL)*DTI*2.0/(H2(IE)+H2(NJP1(IE)))
              CPH1 = 1./(1+CPH)
              WFF(IE,K) = CPH1*(SED(IE,K)+CPH*WFF(NJP1(IE),K))
            ENDIF

          ENDIF

733     CONTINUE
743	CONTINUE
!======================================================================   
!                           END: SEDIMENT OBC
!======================================================================
      
  
      CASE(8)
          
      CASE(9)
************************************************************************
*                              CASE (9)                                *
*                         TEMPORAL CYCLING                             *
************************************************************************  
      IF (NUMEBC.NE.0) THEN
      IF (TRIM(OPTEBC).EQ.'DATA') THEN
      IF (THOUR.GE.T2E) THEN
          T1E = T2E
          DO 240 N = 1, NUMEBC
              DEBDRY(N,1) = DEBDRY(N,2)
240       ENDDO
          READ (IUT90,5000,END=370) T2E
          READ (IUT90,5000) (DEBDRY(N,2),N = 1,NUMEBC)
!$ACC UPDATE DEVICE (DEBDRY)
      ENDIF
      FACT = (THOUR-T1E) / (T2E-T1E)
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO N = 1, NUMEBC
      DO concurrent (N=1:NUMEBC)
           EBDRY(N) = DEBDRY(N,1) + FACT * (DEBDRY(N,2)-DEBDRY(N,1))
      ENDDO
c!$OMP END TARGET
c!$ACC UPDATE DEVICE (EBDRY)
      ENDIF
      ENDIF
      
#ifdef ELB
      IF (NUMEBC.NE.0) THEN
      IF (THOUR.GE.T2EL) THEN
          T1EL = T2EL
          DO 241 N = 1, NUMEBC
              DELBC(N,1) = DELBC(N,2)
241       ENDDO
          READ (IUT95,5000,END=371) T2EL
          READ (IUT95,5000) (DELBC(N,2),N = 1,NUMEBC)
!$ACC UPDATE DEVICE (DELBC)
      ENDIF
      FACT = (THOUR-T1EL) / (T2EL-T1EL)
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO N = 1, NUMEBC
      DO concurrent (N=1:NUMEBC)
            ELBC(N) = DELBC(N,1) + FACT * (DELBC(N,2)-DELBC(N,1))
      ENDDO
c!$OMP END TARGET
c!$ACC UPDATE DEVICE (ELBC)
      ENDIF
#endif

#if defined MODULE_SAL ||  defined MODULE_TMP
      IF (NUMEBC.NE.0) THEN
      IF (THOUR.GE.T2TS) THEN
          T1TS = T2TS
          DO 270 N = 1, NUMEBC
          DO 260 K = 1, KBM1
              DTBDRY(N,K,1) = DTBDRY(N,K,2)
              DSBDRY(N,K,1) = DSBDRY(N,K,2)
260       ENDDO
270       ENDDO
          READ (IUT94,5000,END=380) T2TS
          DO 280 N = 1, NUMEBC
              READ (IUT94,5000) (DTBDRY(N,K,2),K = 1,KBM1)
              READ (IUT94,5000) (DSBDRY(N,K,2),K = 1,KBM1)
280       ENDDO
!$ACC UPDATE DEVICE (DTBDRY,DSBDRY)
      ENDIF
      FACT = (THOUR-T1TS) / (T2TS-T1TS)
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO N = 1, NUMEBC
c      DO K = 1, KBM1
      DO concurrent (N=1:NUMEBC,K=1:KBM1)
          TBDRY(N,K) = DTBDRY(N,K,1) + FACT * (DTBDRY(N,K,2)-
     *    DTBDRY(N,K,1))
          SBDRY(N,K) = DSBDRY(N,K,1) + FACT * (DSBDRY(N,K,2)-
     *    DSBDRY(N,K,1))
      ENDDO
c      ENDDO
c!$OMP END TARGET
c!$ACC UPDATE DEVICE (TBDRY,SBDRY)
      ENDIF
#endif

#ifdef MODULE_SED      
      IF (NUMQBC.NE.0) THEN
      IF (THOUR.GE.T2SEDQ) THEN
          T1SEDQ = T2SEDQ
          DO 315 N = 1, NUMQBC
              DSEDDIS(N,1) = DSEDDIS(N,2)
315       CONTINUE
          READ (IUT20,5000,END=390) T2SEDQ
          READ (IUT20,5000) (DSEDDIS(N,2),N = 1,NUMQBC)
      ENDIF
      FACT = (THOUR-T1SEDQ) / (T2SEDQ-T1SEDQ)
      DO 325 N = 1, NUMQBC
          SEDDIS(N) = DSEDDIS(N,1) + FACT * (DSEDDIS(N,2)-DSEDDIS(N,1))
325   CONTINUE
	ENDIF
      IF (NUMEBC.NE.0) THEN
      IF (THOUR.GE.T2SED) THEN
          T1SED = T2SED
          DO 271 N = 1, NUMEBC
          DO 261 K = 1, KBM1
              DSEDBDRY(N,K,1) = DSEDBDRY(N,K,2)
261       CONTINUE
271       CONTINUE
          READ (IUT21,5000,END=385) T2SED
          DO 281 N = 1, NUMEBC
	        READ (IUT21,5000) (DSEDBDRY(N,K,2),K = 1,KBM1)
281       CONTINUE
      ENDIF
      FACT = (THOUR-T1SED) / (T2SED-T1SED)
      DO 301 N = 1, NUMEBC
      DO 291 K = 1, KBM1
          SEDBDRY(N,K) = DSEDBDRY(N,K,1) + FACT * (DSEDBDRY(N,K,2)-
     *    DSEDBDRY(N,K,1))
291   CONTINUE
301   CONTINUE
      ENDIF
#endif	 
      IF (NUMQBC.NE.0) THEN
      IF (THOUR.GE.T2Q) THEN
          T1Q = T2Q
          DO 310 N = 1, NUMQBC
              DQDIS(N,1) = DQDIS(N,2)
              DTDIS(N,1) = DTDIS(N,2)
              DSDIS(N,1) = DSDIS(N,2)
310       ENDDO
          READ (IUT91,5000,END=390) T2Q
          READ (IUT91,5000) (DQDIS(N,2),N = 1,NUMQBC)
          READ (IUT91,5000) (DTDIS(N,2),N = 1,NUMQBC)
          READ (IUT91,5000) (DSDIS(N,2),N = 1,NUMQBC)
!$ACC UPDATE DEVICE (DQDIS,DTDIS,DSDIS)
      ENDIF
      FACT = (THOUR-T1Q) / (T2Q-T1Q)
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO N = 1, NUMQBC
      DO concurrent (N=1:NUMQBC)
          QDIS(N) = DQDIS(N,1) + FACT * (DQDIS(N,2)-DQDIS(N,1))
          TDIS(N) = DTDIS(N,1) + FACT * (DTDIS(N,2)-DTDIS(N,1))
          SDIS(N) = DSDIS(N,1) + FACT * (DSDIS(N,2)-DSDIS(N,1))
      ENDDO
c!$OMP END TARGET
c!$ACC UPDATE DEVICE (QDIS,TDIS,SDIS)
      ENDIF
#ifdef MODULE_MATERIAL
      IF (THOUR.GE.T2MAT) THEN
          T1MAT = T2MAT
          DO N = 1, NUMMAT !V2410
              DMDIS(N,1) = DMDIS(N,2)
          ENDDO
          READ (IUT98,5000,END=390) T2MAT
          READ (IUT98,5000) (DMDIS(N,2),N = 1,NUMMAT) !V2410
!$ACC UPDATE DEVICE (DMDIS)
      ENDIF 
      FACT = (THOUR-T1MAT) / (T2MAT-T1MAT) !V2410
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO N = 1, NUMMAT !V2410
      DO concurrent (N=1:NUMMAT) !V2410
          MDIS(N)=DMDIS(N,1)+FACT*(DMDIS(N,2)-DMDIS(N,1))
      ENDDO
c!$OMP END TARGET
#endif
      
#ifdef WDTYP_UNIFORM
      IF (THOUR.GE.T2M) THEN
          T1M = T2M
          DQPREC(1) = DQPREC(2)
          DQEVAP(1) = DQEVAP(2)
          DTX(1) = DTX(2)
          DTY(1) = DTY(2)
          DHFLUX(1) = DHFLUX(2)
          READ (IUT93,END=410) T2M
          READ (IUT93) DQPREC(2), DQEVAP(2), DTX(2), DTY(2),DHFLUX(2)
      ENDIF
      FACT = (THOUR-T1M) / (T2M-T1M)
      QPREC = DQPREC(1) + FACT * (DQPREC(2)-DQPREC(1))
      QEVAP = DQEVAP(1) + FACT * (DQEVAP(2)-DQEVAP(1))
      TX = DTX(1) + FACT * (DTX(2)-DTX(1))
      TY = DTY(1) + FACT * (DTY(2)-DTY(1))
      WSTR=SQRT(TX**2+TY**2)
      CD = 1.2E-3
      WDS=SQRT(WSTR/(1.2*CD)) !WDS<11
      IF (WSTR.GE.0.1750) THEN !WDS>=11
      DO K=1,10
          CD = (0.49+0.065*WDS) * 1.E-3
	    WDS=SQRT(WSTR/(1.2*CD))
      ENDDO
      ENDIF
      IF (WSTR.GE.1.5862) THEN  !WDS>=25
          CD = (0.49+0.065*25.) * 1.E-3
          WDS=SQRT(WSTR/(1.2*CD))
      ENDIF
      WINDU = TX/(1.2*CD)/WDS
      WINDV = TY/(1.2*CD)/WDS
      HFLUX = DHFLUX(1) + FACT * (DHFLUX(2)-DHFLUX(1))
      SPCP = 4.186E3
      ROSEA = 1.0E3
      SPRO = SPCP*ROSEA
      PI = 3.1415926 
      TTIME = FLOAT(NSTEP)*DTI/3600.
      TSDAY = 0.0*DTI/3600.
      DO I = 1, N_CTRD
          WUSURF(I) = -1.E-3*H2(I)*(TX*XXI(I)+TY*YXI(I))/DJ(I)
          WVSURF(I) = -1.E-3*H1(I)*(TX*XETA(I)+TY*YETA(I))/DJ(I)
          IF(TTIME.LT.TSDAY) THEN
              WTSURF(I) = 0.0
              SWRAD (I) = 0.0
          ELSE
              WTSURF(I) = 0.
              SWRAD (I) = 0.
          ENDIF
      ENDDO
#elif defined WDTYP_FIELD
      IF (THOUR.GE.T2M) THEN
          T1M = T2M
          DO I = 1,N_CTRD
              DTXX(I,1) = DTXX(I,2)
              DTYY(I,1) = DTYY(I,2)
          ENDDO
          READ (IUT93,END=410) T2M
          READ (IUT93) TXX,TYY
          DO I = 1,N_CTRD
              DTXX(I,2) = TXX(I)
              DTYY(I,2) = TYY(I)
          ENDDO
!$ACC UPDATE DEVICE (DTXX,DTYY)
      ENDIF
      FACT = (THOUR-T1M) / (T2M-T1M)    
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I = 1, N_CTRD
      DO concurrent (I=1:N_CTRD) local(TX,TY,WSTR,CD,WDS)
          TX = DTXX(I,1) + FACT * (DTXX(I,2)-DTXX(I,1))
          TY = DTYY(I,1) + FACT * (DTYY(I,2)-DTYY(I,1))
          WUSURF(I) = -1.E-3*H2(I)*(TX*XXI(I)+TY*YXI(I))/DJ(I)
          WVSURF(I) = -1.E-3*H1(I)*(TX*XETA(I)+TY*YETA(I))/DJ(I)
          WSTR = SQRT(TX**2+TY**2)
          
c          CD = 1.2E-3
c          WDS = SQRT(WSTR/(1.2*CD)) !WDS<11
c          IF (WSTR.GE.0.1750) THEN !WDS>=11
c              CD = (0.49+0.065*WDS) * 1.E-3
c              WDS = SQRT(WSTR/(1.2*CD))
c          ENDIF
c          IF (WSTR.GE.1.5862) THEN  !WDS>=25
c              CD = (0.49+0.065*25.) * 1.E-3
c              WDS = SQRT(WSTR/(1.2*CD))
c          ENDIF 
          IF (WSTR<0.1750) THEN  !NORMAL
              CD = 1.2E-3
          ELSEIF (WSTR<1.5862) THEN  !WDS>=11,<25    
              WDS = SQRT(WSTR/(1.2*1.2E-3))
              CD = (0.49+0.065*WDS) * 1.E-3
          ELSE !WDS>=25
              CD = (0.49+0.065*25.) * 1.E-3
          ENDIF
c          WDS = SQRT(WSTR/(1.2*CD))
c          WDS = 1.2*CD*WDS
          WDS = SQRT(WSTR*1.2*CD)
          IF (WDS.EQ.0.) THEN
			WINDU(I) = 0.0
			WINDV(I) = 0.0
          ELSE
              WINDU(I) = TX/WDS
              WINDV(I) = TY/WDS
          ENDIF
      ENDDO
c!$OMP END TARGET
c!$ACC UPDATE DEVICE (WUSURF,WVSURF,WINDU,WINDV)
#ELSE
      PRINT*, 'WIND DATA TYPE UNdefined!'
	PAUSE
	STOP
#endif
#ifdef HEATFLUX_BULK
      IF (THOUR.GE.T2A) THEN
          T1A = T2A
          DO I=1,N_CTRD
              DATP(I,1)=DATP(I,2)
          ENDDO
          READ (IUT86,END=410) T2A
          READ (IUT86) TXX
          DO I=1,N_CTRD
              DATP(I,2)=TXX(I)
          ENDDO
!$ACC UPDATE DEVICE (DATP)
      ENDIF
      FACT=(THOUR-T1A)/(T2A-T1A)
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I=1,N_CTRD
      DO concurrent (I=1:N_CTRD)
          ATP(I)=DATP(I,1)+FACT*(DATP(I,2)-DATP(I,1))
      ENDDO
c!$OMP END TARGET
      IF (THOUR.GE.T2R) THEN
          T1R = T2R
          DO I=1,N_CTRD
              DRHM(I,1)=DRHM(I,2)
          ENDDO
          READ (IUT87,END=410) T2R
          READ (IUT87) TXX
          DO I=1,N_CTRD
              DRHM(I,2)=TXX(I)
          ENDDO
!$ACC UPDATE DEVICE (DRHM)
      ENDIF
      FACT=(THOUR-T1R)/(T2R-T1R)
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I=1,N_CTRD
      DO concurrent (I=1:N_CTRD)
          RHM(I)=DRHM(I,1)+FACT*(DRHM(I,2)-DRHM(I,1))
      ENDDO
c!$OMP END TARGET
      IF (THOUR.GE.T2C) THEN
          T1C = T2C
          DO I=1,N_CTRD
              DCLOUD(I,1)=DCLOUD(I,2)
          ENDDO
          READ (IUT88,END=410) T2C
          READ (IUT88) TXX
          DO I=1,N_CTRD
              DCLOUD(I,2)=TXX(I)
          ENDDO
!$ACC UPDATE DEVICE (DCLOUD)
      ENDIF
      FACT=(THOUR-T1C)/(T2C-T1C)
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I=1,N_CTRD
      DO concurrent (I=1:N_CTRD)
          CLOUD(I)=DCLOUD(I,1)+FACT*(DCLOUD(I,2)-DCLOUD(I,1))
          CLOUD(I)=CLOUD(I)/10
      ENDDO
c!$OMP END TARGET
c!$ACC UPDATE DEVICE (ATP,RHM,CLOUD)
#endif
#if defined HEATFLUX_BULK  ||  defined AIRPRESSURE 
      IF (THOUR.GE.T2APR) THEN
          T1APR=T2APR
          DO I=1,N_CTRD
              DAPR(I,1)=DAPR(I,2)
          ENDDO
          READ (IUT89,END=410) T2APR
          READ (IUT89) TXX
          DO I=1,N_CTRD
              DAPR(I,2)=TXX(I)
          ENDDO
!$ACC UPDATE DEVICE (DAPR)
      ENDIF
      FACT=(THOUR-T1APR) / (T2APR-T1APR)
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I=1,N_CTRD
      DO concurrent (I=1:N_CTRD)
          APR(I)=DAPR(I,1)+FACT*(DAPR(I,2)-DAPR(I,1))
      ENDDO
c!$OMP END TARGET
c!$ACC UPDATE DEVICE (APR)
#endif

#ifdef MODULE_WAVE
      IF (THOUR.LT.T2WAVE) GOTO 8010
          T1WAVE=T2WAVE
          DO 8020 I=1,N_CTRD
              DWHT (I,1)=DWHT (I,2) 
              DWPER(I,1)=DWPER(I,2) 
              DWDIR(I,1)=DWDIR(I,2) 
 8020     CONTINUE
          READ (IUT92,END=5016) T2WAVE
          READ (IUT92) (DWHT (I,2),I=1,N_CTRD)
          READ (IUT92) (DWPER(I,2),I=1,N_CTRD)
          READ (IUT92) (DWDIR(I,2),I=1,N_CTRD)
 8010     CONTINUE
          FACT=(THOUR-T1WAVE)/(T2WAVE-T1WAVE)
          DO 8030 I=1,N_CTRD
              HSIG(I)=DWHT(I,1)+FACT*(DWHT(I,2)-DWHT(I,1))
              TSIG(I)=DWPER(I,1)+FACT*(DWPER(I,2)-DWPER(I,1))
              FACT1=DWDIR(I,2)-DWDIR(I,1)
              IF (FACT1.GT.180.) THEN
                  FACT1=FACT1-360.
              ELSEIF(FACT1.LT.-180.) THEN
                  FACT1=FACT1+360.
              ENDIF
              WAVEDIR(I)=DWDIR(I,1)+FACT*FACT1
              IF(WAVEDIR(I).LT.0.) WAVEDIR(I)=WAVEDIR(I)+360
              IF(WAVEDIR(I).GT.360.) WAVEDIR(I)=WAVEDIR(I)-360

              HSIG(I)=RAMP*HSIG(I)
 8030     CONTINUE
#endif
#ifdef TIDE_FLATHER
      IF (NUMEBC.NE.0) THEN
      IF (THOUR.GE.T2F) THEN
          T1F = T2F
          DO N = 1, NUMEBC
          DO K=1,KBM1
              DVFOBC(N,K,1) = DVFOBC(N,K,2)
          ENDDO
          ENDDO
          READ (IUT96,*,END=390) T2F
          DO N=1,NUMEBC
              READ (IUT96,*) (DVFOBC(N,K,2), K=1,KBM1)
          ENDDO
!$ACC UPDATE DEVICE (DVFOBC)
      ENDIF
      FACT = (THOUR-T1F) / (T2F-T1F)
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO N=1,NUMEBC
c      DO K=1,KBM1
      DO concurrent (N=1:NUMEBC,K=1:KBM1)
          VFOBC(N,K)=DVFOBC(N,K,1)+FACT*(DVFOBC(N,K,2)-DVFOBC(N,K,1))
      ENDDO
c      ENDDO
c!$OMP END TARGET
      ENDIF
c!$ACC UPDATE DEVICE (VFOBC)
#endif


!======================================================================   
!                       END: TEMPORAL CYCLING
!======================================================================
      
  
      CASE(10)
!======================================================================   
!                     ELEVATION RADIATION OBC
!======================================================================
!600    CONTINUE
C     THIS CONDITION HAS NOT BE USED IN THE COASTAL TIDAL MODEL CASE SINCE
C     THE SEA LEVEL IS GIVEN AT THE OPEN BOUNDARIES. IT CAN BE OPEN IF WE
C     STUDY THE WAVE PROPAGATION PROBLEM WITH ONE FREE OPEN BOUNDARY. IT
C     SHOULD BE EASY TO CHANGE RADIATION BOUNDARY CONDITION FOR A
C     SPECIFIED BOUNDARY.

       DO 81 N = 1,NUMEBC
         IE = NETA(N)
         IC = NCON(N)
         IF (H(IC).GT.0) THEN
           
#ifdef LBCEL_GRA
          !GRADIENT BOUNDARY CONDITION
          ELF(IE) = ELF(IC)
          
#elif defined LBCEL_CHA
           !CHAPMAN BOUNDARY CONDITION
           CP = SQRT(GRAV*H(IC))*DTI/H1(IC)
           CP1 = 1./(1.+CP)
           ELF(IE) = CP1*(EL(IE)+CP*ELF(IC))
           
#elif defined LBCEL_RAD
           !IMPLICIT UPSTREAM RADIATION CONDITION
           IF (IE.EQ.NIP1(IC)) THEN  !EAST BOUNDARY
             IF (FSM(IE)*FSM(NJP1(IE))*FSM(NJM1(IE)).EQ.0.) THEN
               ELF(IE) = ELF(NIM1(IE))
             ELSE
               DO J = 0,1
                 DO I = -1,0
                   NTEMP = FINDCTRD(IE,I,J)
                   GRAD(NTEMP) = EL(NTEMP)-EL(NJM1(NTEMP))
                 ENDDO
               ENDDO
               
               DUDT = EL(NIM1(IE))-ELF(NIM1(IE))
               DUDX = ELF(NIM1(IE))-ELF(NIM2(IE))
             
               IF ((DUDT*DUDX).LT.0.0) DUDT = 0.0
               IF 
     * ((DUDT*(GRAD(NIM1(IE))+GRAD(FINDCTRD(IE,-1,1)))).GT.0.0) THEN
                 DUDE = GRAD(NIM1(IE))
               ELSE 
                 DUDE = GRAD(FINDCTRD(IE,-1,1))
               ENDIF
               CFF = MAX(DUDX*DUDX+DUDE*DUDE,EPS)
               CX = DUDT*DUDX
               CE = MIN(CFF,MAX(DUDT*DUDE,-CFF))
               
               ELF(IE) = (CFF*EL(IE)+CX*ELF(NIM1(IE)) 
     *                        -MAX(CE,0.0)*GRAD(IE)
     *                        -MIN(CE,0.0)*GRAD(NJP1(IE)))/(CFF+CX)
             ENDIF
             
           ELSEIF (IE.EQ.NIM1(IC)) THEN  !WEST BOUNDARY
             IF (FSM(IE)*FSM(NJP1(IE))*FSM(NJM1(IE)).EQ.0.) THEN
               ELF(IE) = ELF(NIP1(IE))
             ELSE
               DO J = 0,1
                 DO I = 0,1
                   NTEMP = FINDCTRD(IE,I,J)
                   GRAD(NTEMP) = EL(NTEMP)-EL(NJM1(NTEMP))
                 ENDDO
               ENDDO
             
               DVDT = EL(NIP1(IE))-ELF(NIP1(IE))
               DVDX = ELF(NIP1(IE))-ELF(NIP2(IE))
             
               IF ((DVDT*DVDX).LT.0.0) DVDT = 0.0
               IF 
     * ((DVDT*(GRAD(NIP1(IE))+GRAD(FINDCTRD(IE,1,1)))).GT.0.0) THEN
                 DVDE = GRAD(NIP1(IE))
               ELSE 
                 DVDE = GRAD(FINDCTRD(IE,1,1))
               ENDIF
               CFF = MAX(DVDX*DVDX+DVDE*DVDE,EPS)
               CX = DVDT*DVDX
               CE = MIN(CFF,MAX(DVDT*DVDE,-CFF))
               
               ELF(IE) = (CFF*EL(IE)+CX*ELF(NIP1(IE)) 
     *                        -MAX(CE,0.0)*GRAD(IE)
     *                        -MIN(CE,0.0)*GRAD(NJP1(IE)))/(CFF+CX)
             ENDIF
             
           ELSEIF (IE.EQ.NJP1(IC)) THEN  !NORTH BOUNDARY
             IF (FSM(NIM1(IE))*FSM(IE)*FSM(NIP1(IE)).EQ.0.) THEN
               ELF(IE) = ELF(NJM1(IE))
               
             ELSE
               DO J = -1,0
                 DO I = -1,0
                   NTEMP = FINDCTRD(IE,I,J)
                   GRAD(NTEMP) = EL(NIP1(NTEMP))-EL(NTEMP)
                 ENDDO
               ENDDO
             
               DUDT = EL(NJM1(IE))-ELF(NJM1(IE))
               DUDE = ELF(NJM1(IE))-ELF(NJM2(IE))
             
               IF ((DUDT*DUDE).LT.0.0) DUDT = 0.0
               IF 
     * ((DUDT*(GRAD(FINDCTRD(IE,-1,-1))+GRAD(NJM1(IE)))).GT.0.0) THEN
                 DUDX = GRAD(FINDCTRD(IE,-1,-1))
               ELSE 
                 DUDX = GRAD(NJM1(IE))
               ENDIF
               CFF = MAX(DUDX*DUDX+DUDE*DUDE,EPS)
               CX = MIN(CFF,MAX(DUDT*DUDX,-CFF))
               CE = DUDT*DUDE
               
               ELF(IE) = (CFF*EL(IE)+CE*ELF(NJM1(IE)) 
     *                        -MAX(CX,0.0)*GRAD(NIM1(IE))
     *                        -MIN(CX,0.0)*GRAD(IE))/(CFF+CE)
             ENDIF
           
           ELSEIF (IE.EQ.NJM1(IC)) THEN  !SOUTH BOUNDARY
             IF (FSM(NIM1(IE))*FSM(IE)*FSM(NIP1(IE)).EQ.0.) THEN
               ELF(IE) = ELF(NJP1(IE))
               
             ELSE
               DO J = 0,1
                 DO I = -1,0
                   NTEMP = FINDCTRD(IE,I,JE)
                   GRAD(NTEMP) = EL(NIP1(NTEMP))-EL(NTEMP)
                 ENDDO
               ENDDO
             
               DUDT = EL(NJP1(IE))-ELF(NJP1(IE))
               DUDE = ELF(NJP1(IE))-ELF(NJP2(IE))
             
               IF ((DUDT*DUDE).LT.0.0) DUDT = 0.0
               IF 
     * ((DUDT*(GRAD(FINDCTRD(IE,-1,1))+GRAD(NJP1(IE)))).GT.0.0) THEN
                 DUDX = GRAD(FINDCTRD(IE,-1,1))
               ELSE 
                 DUDX = GRAD(NJP1(IE))
               ENDIF
               CFF = MAX(DUDX*DUDX+DUDE*DUDE,EPS)
               CX = MIN(CFF,MAX(DUDT*DUDX,-CFF))
               CE = DUDT*DUDE
               
               ELF(IE) = (CFF*EL(IE)+CE*ELF(NJP1(IE)) 
     *                        -MAX(CX,0.0)*GRAD(NIM1(IE))
     *                        -MIN(CX,0.0)*GRAD(IE))/(CFF+CE)
             ENDIF
             
           ENDIF
           
#endif
           
      ENDIF
81    CONTINUE
       
!======================================================================   
!                     END: ELEVATION RADIATION OBC
!====================================================================== 
      END SELECT
      
      RETURN
    
      
370   WRITE (6,5100) THOUR
      GO TO 420
371   WRITE (6,5100) THOUR
      GO TO 420
380   WRITE (6,5200) THOUR
      GO TO 420
385   WRITE (6,5250) THOUR
      GO TO 420 
390   WRITE (6,5300) THOUR
      GO TO 420
400   WRITE (6,5400) THOUR
      GO TO 420
410   WRITE (6,5500) THOUR
      GO TO 420
5016  WRITE(6,5017) THOUR
      GO TO 420
420   CONTINUE
      CLOSE (IUT90)
      CLOSE (IUT91)
      CLOSE (IUT92)
      CLOSE (IUT93)
      CLOSE (IUT94)
      CLOSE (IUT95)
      
 5000 FORMAT (8E14.7)
 5100 FORMAT (//' THE MODEL HAS RUN OUT OF ELEVATION DATA AT TIME 'F10
     *    .4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)

 5017  FORMAT(//' THE MODEL HAS RUN OUT OF WIND WAVE DATA AT TIME ',
     .    F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)
 5200 FORMAT (//' THE MODEL HAS RUN OUT OF TEMPERATURE-SALINITY DATA AT
     *TIME 'F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT',//)
 5250 FORMAT (//' THE MODEL HAS RUN OUT OF SEDIMENT DATA AT
     *TIME 'F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT',//)
 5300 FORMAT (//' THE MODEL HAS RUN OUT OF DISCHARGE DATA AT TIME 'F10
     *    .4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)
 5400 FORMAT (//' THE MODEL HAS RUN OUT OF DIFFUSER DATA AT TIME 'F10.4,
     *    ' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)
 5500 FORMAT (//' THE MODEL HAS RUN OUT OF METEORLOGICAL DATA AT TIME '
     *    F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)
      
      END SUBROUTINE BCOND
