#include "DEFS.h"
      
      SUBROUTINE PROFV_OLD
C     VERSION(03/02/90)
C
      USE MOD_GLOBAL

C**********************************************************************
C                                                                     *
C        THE FOLLOWING SECTION SOLVES THE EQUATION                    *
C         DTI*(KM*V')'-V=-VB                                          *
C                                                                     *
C**********************************************************************
     
        DO K = 1,KB
#ifdef OMP
!$OMP PARALLEL DO
#endif
          DO I = 1,N_CTRD_AG
            C(I,K) = (KM(I,K)+KM(NJM1(I),K))*.5
          ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif
        ENDDO
      
      DO K = 2,KBM1
#ifdef OMP
!$OMP PARALLEL DO
#endif
        DO I = 1,N_CTRD_AG
          IF(DVM(I).GT.0.0) THEN
            A(I,K-1) = -DTI*(C(I,K)+UMOL)/(DZ(K-1)*DZZ(K-1)*DV(I)*DV(I))
            C(I,K) = -DTI*(C(I,K)+UMOL)/(DZ(K)*DZZ(K-1)*DV(I)*DV(I))
          ENDIF
        ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif
      ENDDO

#ifdef OMP
!$OMP PARALLEL DO
#endif
      DO I = 1,N_CTRD_AG
        IF(DVM(I).GT.0.0) THEN
          VH(I,1) = A(I,1)/(A(I,1)-1.)
          VHP(I,1) = (-DTI*.5*(WVSURF(I)+WVSURF(NJM1(I)))
     * /(-DZ(1)*DV(I))-VF(I,1))/(A(I,1)-1.)
        ENDIF
      ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif

      DO K = 2,KBM2
#ifdef OMP
!$OMP PARALLEL DO
#endif
        DO I = 1,N_CTRD_AG
          IF(DVM(I).GT.0.0) THEN
            VHP(I,K) = 1./(A(I,K)+C(I,K)*(1.-VH(I,K-1))-1.)
            VH(I,K) = A(I,K)*VHP(I,K)
            VHP(I,K) = (C(I,K)*VHP(I,K-1)-VF(I,K))*VHP(I,K)
          ENDIF
        ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif
      ENDDO
     
#ifdef OMP
!$OMP PARALLEL DO
#endif
      DO I = 1,N_CTRD_AG
        IF(DVM(I).GT.0.0) THEN
          TPS(I) = CBC_V(I)
     * *SQRT((.25*(U(I,KBM1)+U(NIP1(I),KBM1)
     * +U(NJM1(I),KBM1)+U(NIP1JM1(I),KBM1)))**2+V(I,KBM1)**2)
          
          VF(I,KBM1) = (C(I,KBM1)*VHP(I,KBM2)-VF(I,KBM1))
     * /(TPS(I)*DTI/(-DZ(KBM1)*DV(I))-1.-(VH(I,KBM2)-1.)*C(I,KBM1))
        ENDIF
      ENDDO
#ifdef OMP
!$OMP END PARALLEL DO   
#endif

      DO K = 2,KBM1
        KI = KB-K
#ifdef OMP
!$OMP PARALLEL DO
#endif
        DO I = 1,N_CTRD_AG
          IF(DVM(I).GT.0.0) THEN
            VF(I,KI) = (VH(I,KI)*VF(I,KI+1)+VHP(I,KI))
          ENDIF
        ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif
      ENDDO

      
      
      WVBOT = 0.
#ifdef OMP
!$OMP PARALLEL DO
#endif
      DO I = 1,N_CTRD_AG
        IF(DVM(I).GT.0.0) THEN
          WVBOT(I) = -TPS(I)*VF(I,KBM1)
        ENDIF
      ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif
      
      RETURN
      END  SUBROUTINE PROFV_OLD
