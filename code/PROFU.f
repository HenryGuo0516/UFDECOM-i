#include "DEFS.h"
      
      SUBROUTINE PROFU
      USE MOD_GLOBAL

***********************************************************************
*                                                                     *
*        THE FOLLOWING SECTION SOLVES THE EQUATION                    *
*         DTI*(KM*U')'-U=-UB                                          *
*                                                                     *
***********************************************************************
c      C=0.
c      TPS=0.
c      WUBOT=0.
c      A=0.
c      VH=0.
c      VHP=0.

c Useless:
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO I = 1, N_CTRDP1
c      DO K = 1, KB
c          A(I,K)=0.
c          VH(I,K)=0.
c          VHP(I,K)=0.
c      ENDDO
c      ENDDO
c!$OMP END TARGET
      
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 1,KB
c      DO I = 1,N_CTRD
      DO concurrent (I=1:N_CTRD,K=1:KB)
          C(I,K)=(KM(I,K)+KM(NIM1(I),K)) * .5
      ENDDO
c      ENDDO
c!$OMP END TARGET

c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 2,KBM1
c      DO I = 1,N_CTRD
      DO concurrent (I=1:N_CTRD,K=2:KBM1) local_init(DTI,UMOL)
      IF (DUM(I) .GT. 0.0) THEN
c            C0 = -DTI*(C(I,K)+UMOL)/(DZZ(K-1)*DU(I)*DU(I))
c            A(I,K-1) = C0/DZ(K-1)
c            C(I,K) = C0/DZ(K)
            A(I,K-1) = -DTI*(C(I,K)+UMOL)/(DZ(K-1)*DZZ(K-1)*DU(I)*DU(I))
c            C(I,K) = -DTI*(C(I,K)+UMOL)/(DZ(K)*DZZ(K-1)*DU(I)*DU(I))
            C(I,K) = A(I,K-1)*(DZ(K-1)/DZ(K))
      ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET
      
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I = 1,N_CTRD
      DO concurrent (I=1:N_CTRD) local(K,TPS0) local_init(DTI,KBM1,KBM2)
      IF (DUM(I) .GT. 0.0) THEN
          VH(I,1)=A(I,1)/(A(I,1)-1.)
          VHP(I,1)=(-DTI*.5*(WUSURF(I)+WUSURF(NIM1(I)))
     *    /(-DZ(1)*DU(I))-UF(I,1))/(A(I,1)-1.)
          DO K = 2,KBM2
              VHP(I,K) = 1./(A(I,K)+C(I,K)*(1.-VH(I,K-1))-1.)
              VH(I,K) = A(I,K) * VHP(I,K)
              VHP(I,K) = (C(I,K)*VHP(I,K-1)-UF(I,K))*VHP(I,K)
          ENDDO
          TPS0 = CBC_U(I)*SQRT( U(I,KBM1)**2+((
     *V(I,KBM1)*DVM(I)+V(NIM1JP1(I),KBM1)*DVM(NIM1JP1(I))+
     *V(NIM1(I),KBM1)*DVM(NIM1(I))+V(NJP1(I),KBM1)*DVM(NJP1(I)))/
     *(DVM(I)+DVM(NIM1JP1(I))+DVM(NIM1(I))+DVM(NJP1(I))+1.E-10))**2)
          UF(I,KBM1) = (C(I,KBM1)*VHP(I,KBM2)-UF(I,KBM1))
     *    /(TPS0*DTI/(-DZ(KBM1)*DU(I))-1.-(VH(I,KBM2)-1.)*C(I,KBM1))
c          DO K = 2,KBM1
c              KI = KB-K
c              UF(I,KI)=(VH(I,KI)*UF(I,KI+1)+VHP(I,KI))
c          ENDDO
          DO K = KBM2,1,-1
              UF(I,K)=(VH(I,K)*UF(I,K+1)+VHP(I,K))
          ENDDO
          WUBOT(I) = -TPS0*UF(I,KBM1)
      ELSE
          WUBOT(I) = 0.0
      ENDIF
      ENDDO
c!$OMP END TARGET

      RETURN
      END SUBROUTINE PROFU

