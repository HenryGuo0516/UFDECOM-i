#include "DEFS.h"    
      SUBROUTINE PROFQ
      USE MOD_GLOBAL
      IMPLICIT NONE
      
      INTEGER I,J,K,KI
      REAL,PARAMETER :: A1=0.92, B1=16.6, A2=0.74, B2=10.1, C1=0.08
      REAL,PARAMETER :: E1=1.8, E2=1.33, E3=1.0, E4=0.25
      REAL,PARAMETER :: KAPPA=0.40, SEF=1.0, GEE=9.806
      REAL,PARAMETER :: SMALL=1.E-12, SHIW=0.0, ZSH=0.01
      REAL,PARAMETER :: CONST1= 16.6 ** 0.6666667 * SEF   !大概是 6

      REAL,Parameter :: COEF1 = A2 * (1.-6.*A1/B1)
      REAL,Parameter :: COEF2 = 3. * A2 * B2 + 18. * A1 * A2
      REAL,Parameter :: COEF3 = A1 * (1.-3.*C1-6.*A1/B1)
      REAL,Parameter :: COEF4 = 18. * A1 * A1 + 9. * A1 * A2
      REAL,Parameter :: COEF5 = 9. * A1 * A2

c      REAL COEF1,COEF2,COEF3,COEF4,COEF5
c      REAL CONST1,WDSS
      REAL WDSS
      REAL BETA,CP,CD10,USTAR,UTAU,UTAU2  !LXY        
c      REAL ZW0(N_CTRD) !CBR: Not Used Now
c  CBR: Moved to MOD_Global:
c      REAL GH(N_CTRD,KB),SM(N_CTRD,KB),SH(N_CTRD,KB)
c      REAL KN(N_CTRD,KB),BOYGR(N_CTRD,KB)
      
c      DATA A1, B1, A2, B2, C1 /0.92, 16.6, 0.74, 10.1, 0.08/
c      DATA E1 /1.8/, E2 /1.33/, E3 /1.0/, E4 /0.25/
c      DATA KAPPA /0.40/, SEF /1./
c      DATA GEE /9.806/
c      DATA SMALL /1.E-12/
c      DATA SHIW /0.0/, ZSH /0.01/    !LXY
      
c      ACB_VARIAT = 0. !CBR: Not Used Now
c      VARIAT1 = 0. !CBR: Not Used Now
c      ZW0 = 0. !CBR: Not Used Now
         
               
C**********************************************************************
C                                                                     *
C        THE FOLLOWING SECTION SOLVES THE EQUATION                    *
C        DTI*(KQ*Q2')' - Q2*(2.*DTI*DTEF+1.) = -Q2B                   *
C                                                                     *
C**********************************************************************
       
c       CONST1 = 16.6 ** 0.6666667 * SEF   !大概是 6
c       VARIAT1 = CONST1 !CBR: Not Used Now
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 1,KB
c      DO I = 1,N_CTRD_AG
      DO concurrent (I=1:N_CTRD_AG,K=1:KB)
          Q2(I,K) = ABS(Q2(I,K))
          Q2L(I,K) = ABS(Q2L(I,K))
      ENDDO
c      ENDDO
c!$OMP END TARGET 
      
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 2,KBM1
c      DO I = 1,N_CTRD_AG
      DO concurrent (I=1:N_CTRD_AG,K=2:KBM1) local_init(UMOL)
      IF(FSMADD(I).GT.0.0) THEN
          A(I,K)=-DTI*(KQ(I,K+1)+KQ(I,K)+2.*UMOL)*.5
     *    /(DZZ(K-1)*DZ(K)*D(I)*D(I))
          C(I,K)=-DTI*(KQ(I,K-1)+KQ(I,K)+2.*UMOL)*.5
     *    /(DZZ(K-1)*DZ(K-1)*D(I)*D(I))
          DTEF(I,K) = Q2(I,K)*SQRT(Q2(I,K))/(B1*Q2L(I,K)+SMALL)
          BOYGR(I,K) = GEE*(RHO(I,K-1)-RHO(I,K))/(DZZ(K-1)*D(I))
c          PROD(I,K) = 
c     * KM(I,K)*000.*(.5*(-BOYGR(I,K)+ABS(BOYGR(I,K))))**1.5
          PROD(I,K) = 0.0 !CBR: 上式中乘以000.，暂时直接设为0.，待需要时恢复。
          PROD(I,K) = PROD(I,K)+KM(I,K)*.25*SEF
     *    *((U(I,K)-U(I,K-1)+U(NIP1(I),K)-U(NIP1(I),K-1))**2
     *    +(V(I,K)-V(I,K-1)+V(NJP1(I),K)-V(NJP1(I),K-1))**2)
     *    /(DZZ(K-1)*D(I))**2
          PROD(I,K) = PROD(I,K)+KH(I,K)*BOYGR(I,K)
      ENDIF  
      ENDDO
c      ENDDO
c!$OMP END TARGET 

c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I = 1,N_CTRD_AG
      DO concurrent (I=1:N_CTRD_AG) local_init(KBM1,DTI) local(K)
      IF(FSMADD(I).GT.0.0) THEN
c ------ surface and bottom:          
c         UTAU2 =0.25*SQRT(WUSURF(I)**2+WVSURF(I)**2)*CONST1
         VHP(I,1) =0.25*SQRT(WUSURF(I)**2+WVSURF(I)**2)*CONST1
         VH(I,1) = 0.0
         UF(I,KB) = SQRT((.5*(WUBOT(I)+WUBOT(NIP1(I))))**2
     * +(.5*(WVBOT(I)+WVBOT(NJP1(I))))**2)*CONST1
c ------ mid layers:
         DO K = 2,KBM1
           VHP(I,K)=1./(A(I,K)+C(I,K)*(1.-VH(I,K-1))-
     * (2.*DTI*DTEF(I,K)+1.))
           VH(I,K)=A(I,K)*VHP(I,K)
           VHP(I,K)=(-2.*DTI*PROD(I,K)+C(I,K)*VHP(I,K-1)-UF(I,K))
     * *VHP(I,K)
         ENDDO
c ------ mid layers:
c         DO K = 1,KBM1
c            KI = KB-K
c            UF(I,KI) = VH(I,KI)*UF(I,KI+1)+VHP(I,KI)
c         ENDDO
         DO K = KBM1,1,-1
            UF(I,K) = VH(I,K)*UF(I,K+1)+VHP(I,K)
         ENDDO
C**********************************************************************
C                                                                     *
C        THE FOLLOWING SECTION SOLVES THE EQUATION                    *
C        DTI*(KQ*Q2L')' - Q2L*(DTI*DTEF+1.) = -Q2LB                   *
C                                                                     *
C**********************************************************************
c ------ surface:          
          VH(I,1) = 0.0
          VHP(I,1) = 0.0 
c ------ mid layers:
          DO K = 2,KBM1
           DTEF(I,K) = DTEF(I,K)*(1.+E2*(L(I,K)/(ABS(Z(K)-Z(KB))*D(I)
     * *KAPPA))**2+E4*(L(I,K)/(ABS(Z(K)-Z(1))*D(I)*KAPPA))**2)
           VHP(I,K)=1./(A(I,K)+C(I,K)*(1.-VH(I,K-1))-(DTI*DTEF(I,K)+1.))
           VH(I,K) = A(I,K)*VHP(I,K)
           VHP(I,K)=(DTI*(-PROD(I,K)*L(I,K)*E1)+C(I,K)*VHP(I,K-1)
     * -VF(I,K))*VHP(I,K)
          ENDDO
c ------ mid layers:
c          DO K = 1,KBM1
c            KI = KB-K
c            VF(I,KI) = VH(I,KI)*VF(I,KI+1)+VHP(I,KI)
c          ENDDO
          DO K = KBM1,1,-1
            VF(I,K) = VH(I,K)*VF(I,K+1)+VHP(I,K)
          ENDDO
      ENDIF
      ENDDO
c!$OMP END TARGET 
       
       
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 2,KBM1
c      DO I = 1,N_CTRD_AG
      DO concurrent (I=1:N_CTRD_AG,K=2:KBM1)
      IF(UF(I,K).LE.SMALL.OR.VF(I,K).LE.SMALL) THEN
             UF(I,K) = SMALL
             VF(I,K) = SMALL
      ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET 
	   
C**********************************************************************
C                                                                     *
C               THE FOLLOWING SECTION SOLVES FOR KM AND KH            *
C                                                                     *
C**********************************************************************
c      COEF1 = A2 * (1.-6.*A1/B1)
c      COEF2 = 3. * A2 * B2 + 18. * A1 * A2
c      COEF3 = A1 * (1.-3.*C1-6.*A1/B1)
c      COEF4 = 18. * A1 * A1 + 9. * A1 * A2
c      COEF5 = 9. * A1 * A2
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 2,KBM1
c      DO I = 1,N_CTRD_AG
      DO concurrent (I=1:N_CTRD_AG,K=2:KBM1)
      IF(FSMADD(I).GT.0.0) THEN
c ----- L, GH:
          L(I,K) = VF(I,K)/UF(I,K)
          IF(L(I,K)**2*BOYGR(I,K).LT.-0.281*UF(I,K)) THEN
              VF(I,K) = UF(I,K)
     *    *SQRT(ABS(-0.281*UF(I,K)/(BOYGR(I,K)+SMALL)))
          L(I,K) = VF(I,K)/UF(I,K)
          ENDIF
c          L(I,K) = AMAX1(L(I,K),KAPPA*ZW0(I)) !CBR: ZW0 Not Used Now
          L(I,K) = AMAX1(L(I,K),0.0) !CBR: ZW0 Not Used Now
          GH(I,K) = L(I,K)**2/UF(I,K)*GEE*(RHO(I,K)-RHO(I,K-1))
     * /(-DZZ(K-1)*D(I))
          GH(I,K) = AMIN1(GH(I,K),.028)
          GH(I,K) = AMAX1(GH(I,K),-.28)
c ----- SH, SM:
          SH(I,K) = COEF1/(1.-COEF2*GH(I,K))
          SM(I,K) = (COEF3+SH(I,K)*COEF4*GH(I,K))/(1.-COEF5*GH(I,K))
c ----- KN, KQ, KH, KM:
          KN(I,K) = L(I,K)*SQRT(ABS(Q2(I,K)))
          KQ(I,K) = (KN(I,K)*.41*SH(I,K)+KQ(I,K))*.5          
          KM(I,K) = (KN(I,K)*SM(I,K)+KM(I,K))*.5
          KH(I,K) = (KN(I,K)*SH(I,K)+KH(I,K))*.5           
          KM(I,K) = KM(I,K)*KM_COF(I)
          KH(I,K) = KH(I,K)*KH_COF(I)
      ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET 
      
      RETURN
      END SUBROUTINE PROFQ