#include "DEFS.h"
	
	MODULE MOD_SED
	
      USE MOD_GLOBAL
#ifdef MODULE_WAVE
      USE MOD_WAVE
#endif
      SAVE
	CONTAINS
 
!======================================================================
!                                  INITSED
!======================================================================	
      SUBROUTINE INITSED
C     VERSION(3/20/2006)
C COMMON:
C     THIS SOUBROUTINE IF FOR SEDIMENT
C     SED IS GIVEN IN SOME STANDARD LEVEL
C     WE CALCULATE THE RMEAN IN U AND V GRID

      USE MOD_GLOBAL
      INTEGER :: KSL1
      REAL :: DP
      SED = 0.1
      IF (LOG_ISED)THEN
      DO WHILE (.TRUE.)
          READ (IUISED,*,END=10) II,(SED(II,K),K=1,KBM1)
      ENDDO
10    CONTINUE
      CLOSE(IUISED)
      ENDIF
              
      D50 = F_SD50
      TAUE0 = 0.5
      TAUD0 = TAUE0*4./9.
	DO I = 1,N_CTRD
      IF(D(I).GE.0.) THEN
          IF(D50(I).LE.0.0005) DP = 0.0005
          IF(D50(I).GE.0.0005.AND.D50(I).LE.0.01) DP = D50(I)
          IF(D50(I).GE.0.01) DP = 0.01
          USIN(I) = 0.128*((DP/0.01)**(1./6.))
     *    *SQRT(3.6*1.65*9.81*DP+(0.85)**2.5*(1.75E-6/DP+9.81*D(I)*
     *    2.31E-7*SQRT(2.31E-7/DP)/DP))      
	    TAUE0(I) = USIN(I)**2*1025.  !测试表明此值对水深不敏感，计算最小为0.1967,最大为9.7879
          TAUD0(I) = 4.*TAUE0(I)/9.
      ENDIF
      ENDDO
      TAUE = TAUE0 !LXY
      TAUD = TAUD0
      RETURN
	
      END SUBROUTINE INITSED
!======================================================================
!						       END: INITSED
!======================================================================	
      
      
!======================================================================
!							    INIT_TAUE
!======================================================================	
      SUBROUTINE INIT_TAUE
	
      USE MOD_GLOBAL
      INTEGER :: KK,NBOUND,K,IND
      REAL :: COF
      REAL :: XB(10000),YB(10000)
      REAL, ALLOCATABLE :: XBOUND(:)
      REAL, ALLOCATABLE :: YBOUND(:)
      CHARACTER*160 REGION
      
      DO WHILE (.TRUE.)
          READ (IUD50,*,END=117) II,D50(II),TAUE0(II),TAUD0(II) !LXY
      ENDDO
117	CONTINUE
	TAUE0(N_CTRDP1) = -99999.
	TAUD0(N_CTRDP1) = -99999.
      TAUE = TAUE0 
      TAUD = TAUD0
      RETURN
	END SUBROUTINE INIT_TAUE
!======================================================================
!						      END: INIT_TAUE
!======================================================================
            
      
!======================================================================
!							    BOTSTRESS
!======================================================================
	SUBROUTINE BOTSTRESS
	USE MOD_GLOBAL
      INTEGER :: DATE_TIME(8)
      CHARACTER (LEN = 12) REAL_CLOCK (3)
      

	DO I = 1,N_CTRD
	    UBAR(I,KBM1) = 0.5*(U(I,KBM1)+U(NIP1(I),KBM1))!笛卡尔坐标 !CORRECT !LXY
          VBAR(I,KBM1) = 0.5*(V(I,KBM1)+V(NJP1(I),KBM1))
      ENDDO
      DO I = 1, N_CTRD
	    TAU_U(I) = CBC_UV(I)*SQRT(UR(I,KBM1)**2+VR(I,KBM1)**2)
     *    *UR(I,KBM1)*(RHO(I,KBM1)*1000.+1025.)
          TAU_V(I) = CBC_UV(I)*SQRT(UR(I,KBM1)**2+VR(I,KBM1)**2)
     *    *VR(I,KBM1)*(RHO(I,KBM1)*1000.+1025.)
      ENDDO 
      DO I = 1,N_CTRD
	   TAU(I) = SQRT(((TAU_U(I)+TAU_U(NIP1(I)))/2.)**2+
     *  ((TAU_V(I)+TAU_V(NJP1(I)))/2.)**2)
      ENDDO
      
#ifdef MODULE_WAVE      
      TAU_TIDE=TAU
      CALL STRESS_OMP
#else      
	DO I = 1,N_CTRD
      IF (FSM(I).GT.0.0) THEN     
	    TAU_TIDE(I) = TAU(I)
		TAU_WAVE(I) = 0.
      ENDIF
      ENDDO
      DO N = 1,NUMQBC
          IC = NQC(N)
          TAU(IC) = 0.0
      ENDDO         
#endif
      RETURN
	END SUBROUTINE BOTSTRESS
!======================================================================
!							   END: BOTSTRESS
!======================================================================
      
	
	
C**********************************************************************
C********************  SEDDIMENT MODEL SUBROUTINE  *********************
C***********************  EDITED BY LIUGAOFENG  ************************
C***********************  2006.3.17~2006.3.40  *************************

!======================================================================
!							    SEDW
!======================================================================
C-------------计算 泥沙沉降速度 DWS(I,J,K)----------------------------

      SUBROUTINE   SEDW
      USE MOD_GLOBAL
	REAL FTMP,DWS_FREE
      REAL RONGZHONG
      REAL DZZ0
      
	
	RONGZHONG = 1.65  !!! LXY 其实是密度 它乘以重力加速度就是容重
                 !!! 并且这个还是水中实际容重
	DWS = 0.0      
      T=10.!由原来的15改成10，匹配陈曦数据
      DO K = 1,KBM1
	DO I = 1,N_CTRD
      IF(FSM(I).GT.0 .AND. SED(I,K).GT.0) THEN		
      FTMP = SED(I,K)
	!DWS_FREE = 0.0171e-03
      VNIAN(I,K) = 1.792D-6/(1+0.03368*T(I,K)+0.000221*T(I,K)**2)
	DWS_FREE = SQRT(1.09*RONGZHONG*GRAV*D50(I)+
     *     (13.95*VNIAN(I,K)/D50(I))**2)-13.95*VNIAN(I,K)/D50(I)
	DWS(I,K) = DWS_FREE   
     !! IF (FTMP<0.2) THEN
	    !!DWS(I,K) = DWS_FREE
     !! ELSEIF (FTMP >= 0.2 .AND. FTMP<=6.0) THEN	
	    !!DWS(I,K) = (0.012*FTMP**2.2)/(FTMP**2.0+1.7**2.0)**2.8
     !! ELSE
	    !!DWS(I,K) = (0.012*6.0**2.2)/(6.0**2.0+1.7**2.0)**2.8
     !! ENDIF
      IF (FTMP<0.5) THEN
	    DWS(I,K) = DWS_FREE
      ELSEIF (FTMP >= 0.5 .AND. FTMP<=6.0) THEN	
	    DWS(I,K) = (0.002002*FTMP**1.369)/(FTMP**2.0+1.56**2.0)**1.278
      ELSE
	    DWS(I,K) = (0.002002*6.0**1.369)/(6.0**2.0+1.56**2.0)**1.278
      ENDIF
      
      IF (K.LT.KBM1) THEN
          DZZ0=DZZ(K)
      ELSE
          DZZ0=ZZ(K)-(-1.)
      ENDIF
      
      IF (DWS(I,K)*DTI>=D(I)*DZZ0) THEN     !LXY 2014/05/12
	    DWS(I,K) = D(I)*DZZ0/DTI
      ENDIF
      
      ENDIF
      
      IF (K.EQ.KBM1) THEN
          DWS(I,KB) =AMAX1(DWS(I,KBM1),4.0E-4)   !!!DWS(KB)指的是最底边界层沉速                                          !只用来算QDEP !LXY
	ENDIF
	ENDDO
      ENDDO
      
      
      RETURN
	END SUBROUTINE   SEDW
!======================================================================
!							  END: SEDW
!======================================================================


!======================================================================
!							    RESUSPEN
!======================================================================
	SUBROUTINE RESUSPEN(F)
	
      USE MOD_GLOBAL
      !REAL D25,RHOSED  !,DP
      !REAL GAMA,GAMA0
      REAL FTMP
      !REAL AVE_SED,AVE_R
	!REAL F(IM,JM,KB), FF(IM,JM,KB)
      
	REAL F(N_CTRDP1,KB)

	QERO = 0.0
      QDEP = 0.0
      DO I = 1,N_CTRD
      FTMP = F(I,KBM1)
      IF(FSM(I).GT.0.0) THEN
          IF(TAU(I).GE.TAUE(I)) THEN
              QERO(I) = 1*F_MCSXS*(TAU(I)/(TAUE(I)+1.E-30)-1) !冲刷系数
              QDEP(I) = 0.
          ELSEIF (TAU(I).LE.TAUD(I)) THEN
	        QERO(I) = 0.
		    !LXY 沉降概率 !DWS(KB)有单独涵义 
      	    QDEP(I) = 1*F_ALFA*DWS(I,KB)*FTMP*(1-TAU(I)/(TAUD(I)+1.E-30)) 
          !淤积量应不超过底层泥沙总量（包括了上层的下降通量）
	    QDEP(I) = AMIN1(QDEP(I),FTMP*D(I)*DZ(KBM1)/DTI) !LXY 
          
          IF (F(I,KBM1)>=12.3) THEN    !估计不会用到的
            QDEP(I) = AMAX1(QDEP(I),F_MCSXS*FTMP*D(I)*DZ(KBM1)/DTI)
          ENDIF
            
        ELSE
          QERO(I) = 0.
          QDEP(I) = 0.
        ENDIF
      
	  ENDIF     
                         
	ENDDO
      
      !这里的QERO>0，QDEP>0
	
      RETURN
      END  SUBROUTINE RESUSPEN
	  
C+++泥沙的起动速度采用窦国仁的泥沙起动公式——可以适用于粘性，非粘性泥沙的起动+++++ 
!======================================================================
!							  END: RESUSPEN
!======================================================================

        
!======================================================================
!							    RESEABED
!======================================================================
	SUBROUTINE RESEABED 
C     VERSION(10/28/2006)
      USE MOD_GLOBAL
	
      REAL RAMP_E,RAMP_D,RAMPX !!LXY
  
	RGAMA = 0.6    !  THE VOLUMETRIC WEIGHT OF THE RIVER BED
      DO I = 1, N_CTRD
	  IF (FSM(I).EQ.1.0) THEN
        ZBED(I) = (QERO(I)-QDEP(I))/(2650.*RGAMA)*DTI !LXY 2650这个地方是密度
        ZBEDD(I) = ZBEDD(I)+ZBED(I)
	  ENDIF
      ENDDO
      DO I = 1,N_CTRD
	  H(I) = H(I) + ZBED(I)
        D(I) = H(I) + ELF(I)
C       DU(I,J) = AMAX1(0.,HU(I,J)+0.5*(ELF(I,J)+ELF(I-1,J)))
C       DV(I,J) = AMAX1(0.,HV(I,J)+0.5*(ELF(I,J)+ELF(I,J-1)))
      ENDDO
	
      RETURN
	END SUBROUTINE RESEABED
!======================================================================
!							END: RESEABED
!======================================================================
      
      END MODULE 
      
