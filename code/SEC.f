#include "DEFS.h"
      SUBROUTINE SEC
      USE MOD_GLOBAL
      INTEGER N,I,II
      IF (N_SEC_N.NE.0) THEN
!$ACC UPDATE HOST (S,XMFLUX,YMFLUX) !GPU to CPU fo Sub SEC

          
      DO N=1,N_SEC_N
          
      IF (NUM_SEC(1,N).EQ.NIM1(NUM_SEC(2,N)).OR.
     *NUM_SEC(1,N).EQ.NIP1(NUM_SEC(2,N))) THEN
      FLUXSUM(N)=0.
      SFLUXSUM(N)=0.
      SEDFLUXSUM(N)=0.
      DO I=1,TOLNUM_SEC(N)
          II=NUM_SEC(I,N)
          DO K=1,KBM1
              FLUXSUM(N)=FLUXSUM(N)+YMFLUX(II,K)*DZ(K)
#ifdef MODULE_SAL
              IF (THOUR .GE. S_BEG) THEN  
              SFLUXSUM(N)=SFLUXSUM(N)+S(II,K)*YMFLUX(II,K)*DZ(K)
              ENDIF
#endif
#ifdef MODULE_SED
              IF (THOUR .GE. SED_BEG) THEN  
              SEDFLUXSUM(N)=SEDFLUXSUM(N)+SED(II,K)*YMFLUX(II,K)*DZ(K)
              ENDIF
#endif
          ENDDO
      ENDDO
      FLUX_ACM(N)=FLUX_ACM(N)+DBLE(FLUXSUM(N)*DTI)
#ifdef MODULE_SAL
      IF (THOUR .GE. S_BEG) THEN  
          SFLUX_ACM(N)=SFLUX_ACM(N)+DBLE(SFLUXSUM(N)*DTI)
      ENDIF
#endif
#ifdef MODULE_SED
      IF (THOUR .GE. SED_BEG) THEN  
          SEDFLUX_ACM(N)=SEDFLUX_ACM(N)+DBLE(SEDFLUXSUM(N)*DTI)
      ENDIF
#endif
      ELSEIF(NUM_SEC(1,N).EQ.NJM1(NUM_SEC(2,N)).OR.
     *NUM_SEC(1,N).EQ.NJP1(NUM_SEC(2,N))) THEN
      FLUXSUM(N)=0.
      SFLUXSUM(N)=0.
      SEDFLUXSUM(N)=0.
      DO I=1,TOLNUM_SEC(N)
          II=NUM_SEC(I,N)
          DO K=1,KBM1
              FLUXSUM(N)=FLUXSUM(N)+XMFLUX(II,K)*DZ(K)
#ifdef MODULE_SAL
              IF (THOUR .GE.S_BEG) THEN  
              SFLUXSUM(N)=SFLUXSUM(N)+S(II,K)*XMFLUX(II,K)*DZ(K)
              ENDIF
#endif
#ifdef MODULE_SED
              IF (THOUR .GE. SED_BEG) THEN  
              SEDFLUXSUM(N)=SEDFLUXSUM(N)+SED(II,K)*XMFLUX(II,K)*DZ(K)
              ENDIF
#endif
          ENDDO
      ENDDO
      FLUX_ACM(N)=FLUX_ACM(N)+DBLE(FLUXSUM(N)*DTI)
#ifdef MODULE_SAL
      IF (THOUR .GE. S_BEG) THEN  
          SFLUX_ACM(N)=SFLUX_ACM(N)+DBLE(SFLUXSUM(N)*DTI)
      ENDIF
#endif
#ifdef MODULE_SED
      IF (THOUR .GE. SED_BEG) THEN  
          SEDFLUX_ACM(N)=SEDFLUX_ACM(N)+DBLE(SEDFLUXSUM(N)*DTI)
      ENDIF
#endif
      ENDIF
     
      ENDDO
      ENDIF
      END SUBROUTINE SEC