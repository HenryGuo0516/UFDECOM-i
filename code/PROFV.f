#include "DEFS.h"
      
      SUBROUTINE PROFV
      USE MOD_GLOBAL
***********************************************************************
*                                                                     *
*        THE FOLLOWING SECTION SOLVES THE EQUATION                    *
*         DTI*(KM*V')'-V=-VB                                          *
*                                                                     *
***********************************************************************
c      C=0.
c      TPS=0.
c      WVBOT=0.
c      A=0.
c      VH=0.
c      VHP=0.
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 1,KB
c      DO I = 1,N_CTRD
      DO concurrent (I=1:N_CTRD,K=1:KB)
          C(I,K) = (KM(I,K)+KM(NJM1(I),K))*.5
      ENDDO
c      ENDDO
c!$OMP END TARGET


c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 2,KBM1
c      DO I = 1,N_CTRD
      DO concurrent (I=1:N_CTRD,K=2:KBM1) local_init(DTI,UMOL)
      IF(DVM(I).GT.0.0) THEN
c            C0 = -DTI*(C(I,K)+UMOL)/(DZZ(K-1)*DV(I)*DV(I))
c            A(I,K-1) = C0/DZ(K-1)
c            C(I,K) = C0/DZ(K)
          A(I,K-1) = -DTI*(C(I,K)+UMOL)/(DZ(K-1)*DZZ(K-1)*DV(I)*DV(I))
c          C(I,K) = -DTI*(C(I,K)+UMOL)/(DZ(K)*DZZ(K-1)*DV(I)*DV(I))
          C(I,K) = A(I,K-1)*(DZ(K-1)/DZ(K))
      ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET

      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I = 1,N_CTRD
      DO concurrent (I=1:N_CTRD) local(K,TPS0) local_init(DTI,KBM1,KBM2)
      IF(DVM(I).GT.0.0) THEN
          VH(I,1) = A(I,1)/(A(I,1)-1.)
          VHP(I,1) = (-DTI*.5*(WVSURF(I)+WVSURF(NJM1(I)))
     *    /(-DZ(1)*DV(I))-VF(I,1))/(A(I,1)-1.)
          DO K = 2,KBM2
              VHP(I,K) = 1./(A(I,K)+C(I,K)*(1.-VH(I,K-1))-1.)
              VH(I,K) = A(I,K)*VHP(I,K)
              VHP(I,K) = (C(I,K)*VHP(I,K-1)-VF(I,K))*VHP(I,K)
          ENDDO
          TPS0 = CBC_V(I)*SQRT(V(I,KBM1)**2+((
     * U(I,KBM1)*DUM(I)+U(NIP1JM1(I),KBM1)*DUM(NIP1JM1(I))
     * +U(NJM1(I),KBM1)*DUM(NJM1(I))+U(NIP1(I),KBM1)*DUM(NIP1(I)))/
     *(DUM(I)+DUM(NIP1JM1(I))+DUM(NJM1(I))+DUM(NIP1(I))+1.E-10))**2)
          VF(I,KBM1) = (C(I,KBM1)*VHP(I,KBM2)-VF(I,KBM1))
     * /(TPS0*DTI/(-DZ(KBM1)*DV(I))-1.-(VH(I,KBM2)-1.)*C(I,KBM1))
c          DO K = 2,KBM1
c              KI = KB-K
c              VF(I,KI) = (VH(I,KI)*VF(I,KI+1)+VHP(I,KI))
c          ENDDO
          DO K = KBM2,1,-1
              VF(I,K) = (VH(I,K)*VF(I,K+1)+VHP(I,K))
          ENDDO
          WVBOT(I) = -TPS0*VF(I,KBM1)
      ELSE
          WVBOT(I) = 0.0
      ENDIF
      ENDDO
c!$OMP END TARGET

      RETURN
      END  SUBROUTINE PROFV
