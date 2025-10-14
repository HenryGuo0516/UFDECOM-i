#include "DEFS.h" 
***********************************************************************
*     This version is derived by LXY                                  *
*     And Changed into OMP version by MaRin                           *
***********************************************************************      
      SUBROUTINE PROF_SED(S0)  
	USE MOD_GLOBAL

      REAL LL,LKM1,LK
      REAL KH0(KB),DWSBAR(KB),LSK(KB),SSTAR(KB)
      REAL TA(KB),TB(KB),TC(KB),TD(KB),TE(KB),TF(KB)
      REAL AF(KB),BF(KB),DF(KB),VHF(KB),VHPF(KB)
      REAL S0(N_CTRDP1,KB),DD(N_CTRDP1)
      INTEGER KI,K
      DD=D
#ifdef OMP
!$OMP PARALLEL DO PRIVATE(LL,LKM1,LK,K,KI,KH0,DWSBAR,LSK,SSTAR,
     * TA,TB,TC,TD,TE,TF,AF,BF,DF,VHF,VHPF )
#endif
	DO I = 1,N_CTRD
      IF (FSM(I)*FSM11(I).EQ.1.) THEN
      IF (D(I)<=5.*DMIN) THEN
          DD(I)=5.*DMIN
          DO K=1,KB
          KH0(K)=AMAX1(KH(I,K),0.0001)
          ENDDO
      ELSE
          DO K=1,KB
	    KH0(K)=AMAX1(KH(I,K),0.00001)
          ENDDO
      ENDIF
      DO K = 2,KBM1
	    DWSBAR(K) = 0.5*(DWS(I,K-1)+DWS(I,K))
          LL = DWSBAR(K)*DTI
          LKM1 = 0.5*D(I)*DZ(K-1)
          LK = 0.5*D(I)*DZ(K)
          LL = AMIN1(LL,LKM1) 
          LSK(K) = (LKM1-LL)/(LKM1+LK)
      ENDDO
      TA(1) = -DTI*DWSBAR(2)*(1.-LSK(2))/(DD(I)*DZ(1))
      TB(1) = -DTI*DWSBAR(2)*LSK(2)/(DD(I)*DZ(1)) 
      TC(1) = -DTI*KH0(2)/(DD(I)*DD(I)*DZ(1)*DZZ(1))
      AF(1) = -TB(1)+TC(1)
      BF(1) = 1.-TA(1)-TC(1)
      SSTAR(1) = S0(I,1)
      VHF(1) = -AF(1)/BF(1)
      VHPF(1) = SSTAR(1)/BF(1)
      DO K=2,KBM1
          TA(K) = DTI*DWSBAR(K)*(1.-LSK(K))/(DD(I)*DZ(K))
          TB(K) = DTI*DWSBAR(K)*LSK(K)/(DD(I)*DZ(K))  
          TC(K) = -DTI*KH0(K)/(DD(I)*DD(I)*DZ(K)*DZZ(K-1))
          TD(K) = -DTI*DWSBAR(K+1)/(DD(I)*DZ(K))*(1.-LSK(K+1))
          TE(K) = -DTI*DWSBAR(K+1)/(DD(I)*DZ(K))*LSK(K+1)
          TF(K) = -DTI*KH0(K+1)/(DD(I)*DD(I)*DZ(K)*DZZ(K))
          AF(K) = -TE(K)+TF(K)
          BF(K) = 1.-TB(K)-TC(K)-TD(K)-TF(K)
          DF(K) = -TA(K)+TC(K)
          SSTAR(K)=S0(I,K)
      ENDDO
      DO K = 2, KBM2
          VHPF(K) = (BF(K)+DF(K)*VHF(K-1))
          VHF(K) = - AF(K)/VHPF(K)
          VHPF(K) = (SSTAR(K)-DF(K)*VHPF(K-1))/VHPF(K)
      ENDDO
      !DO K = 1, KBM1
            !SS(I,K) = S0(I,K) !HERE,SS=SSTAR
	!ENDDO
      TA(KBM1) = DTI*DWSBAR(KBM1)*(1.-LSK(KBM1))/(DD(I)*DZ(KBM1))
      TB(KBM1) = DTI*DWSBAR(KBM1)*LSK(KBM1)/(DD(I)*DZ(KBM1))
      TC(KBM1) = -DTI*KH0(KBM1)/(DD(I)*DD(I)*DZ(KBM1)*DZZ(KBM2)) 
      BF(KBM1) = 1.-TB(KBM1)-TC(KBM1)
      DF(KBM1) = -TA(KBM1)+TC(KBM1)
      SSTAR(KBM1)=S0(I,KBM1)-(QDEP(I)-QERO(I))*DTI/(DD(I)*DZ(KBM1))
      LL=SSTAR(KBM1)-DF(KBM1)*VHPF(KBM2)
      LK=BF(KBM1)+DF(KBM1)*VHF(KBM2)
      SSTAR=0.
      SSTAR(KBM1) =LL/LK
      SSTAR(KBM1)=AMIN1(SSTAR(KBM1),12.0)
      DO K=2,KBM1
          KI=KB-K
          SSTAR(KI) = (VHF(KI)*SSTAR(KI+1)+VHPF(KI))
      ENDDO
      DO K=1,KBM1
          S0(I,K) = AMAX1(SSTAR(K), 0.) 
      ENDDO
      ENDIF
      ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif
      RETURN
      END SUBROUTINE PROF_SED