#include "DEFS.h"      
      SUBROUTINE BAROPG5B
      USE MOD_GLOBAL,DUZ1=>A,DUZ2=>C,DUZ=>VH,RHOZ=>VHP,
     *RHOS1=>UC,RHOS2=>VC
      INTEGER I,K
c      REAL RHOZ(KBM1),RHOS1(KBM1),RHOS2(KBM1),DUZ(KBM1)
c     * ,DUZ1(KBM1),DUZ2(KBM1)
c#ifdef OMP
c!$OMP PARALLEL DO PRIVATE(RHOZ,RHOS1,RHOS2,DUZ,DUZ1,DUZ2,K)
c#endif      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I = 1,N_CTRD
      DO concurrent (I = 1:N_CTRD) local(K,KI,KJ)
     * local_init(RAMP,KBM1)
      IF (DUM(I).GT.0.0) THEN
      DO K=1,KBM1
          DUZ(I,K) = ZZ(K)*0.5*(D(I)+D(NIM1(I)))+0.5*(EL(I)+EL(NIM1(I)))
          DUZ1(I,K) = ZZ(K)*D(NIM1(I))+EL(NIM1(I))
          DUZ2(I,K) = ZZ(K)*D(I)+EL(I)           
          RHOZ(I,K) = 0.5 * (RHO(I,K)+RHO(NIM1(I),K))      
      ENDDO
c      PURE SUBROUTINE SINTER(X,A,Y,B,M,N) Extracted
c      CALL SINTER(DUZ,RHOZ,DUZ1,RHOS1,KBM1,KBM1) Extracted
c ----- rivision 1:   
      KI=1
      KJ=1
      DO WHILE (KI<=KBM1)
          IF (DUZ1(I,KI)>DUZ(I,1)) THEN 
              RHOS1(I,KI)=RHOZ(I,1)+((RHOZ(I,1)-RHOZ(I,2))/
     *      (DUZ(I,1)-DUZ(I,2)))*(DUZ1(I,KI)-DUZ(I,1))
              KI=KI+1
          ELSEIF (KJ>KBM1) THEN
              RHOS1(I,KI)=RHOZ(I,KBM1)
              KI=KI+1
          ELSEIF (DUZ1(I,KI)>DUZ(I,KJ)) THEN
              RHOS1(I,KI) = RHOZ(I,KJ-1) - (RHOZ(I,KJ-1)-RHOZ(I,KJ)) 
     *        * (DUZ(I,KJ-1)-DUZ1(I,KI)) / (DUZ(I,KJ-1)-DUZ(I,KJ))
              KI=KI+1
          ELSE
              KJ=KJ+1
          ENDIF
      ENDDO
c      CALL SINTER(DUZ,RHOZ,DUZ2,RHOS2,KBM1,KBM1) Extracted
      KI=1
      KJ=1
      DO WHILE (KI<=KBM1)
          IF (DUZ2(I,KI)>DUZ(I,1)) THEN 
              RHOS2(I,KI)=RHOZ(I,1)+((RHOZ(I,1)-RHOZ(I,2))/
     *      (DUZ(I,1)-DUZ(I,2)))*(DUZ2(I,KI)-DUZ(I,1))
              KI=KI+1
          ELSEIF (KJ>KBM1) THEN
              RHOS2(I,KI)=RHOZ(I,KBM1)
              KI=KI+1
          ELSEIF (DUZ2(I,KI)>DUZ(I,KJ)) THEN
              RHOS2(I,KI) = RHOZ(I,KJ-1) - (RHOZ(I,KJ-1)-RHOZ(I,KJ))
     *         * (DUZ(I,KJ-1)-DUZ2(I,KI)) / (DUZ(I,KJ-1)-DUZ(I,KJ))
              KI=KI+1
          ELSE
              KJ=KJ+1
          ENDIF
      ENDDO
      
c --- CBR: Using temp variables RHOS1,RHOS2 instead of RHO itself
c      DO K=1,KBM1
c          RHO(I,K) = RHO(I,K) - RHOS2(I,K)
c          RHO(NIM1(I),K) = RHO(NIM1(I),K) - RHOS1(I,K)
c      ENDDO
      DO K=1,KBM1
          RHOS2(I,K) = RHO(I,K) - RHOS2(I,K)
          RHOS1(I,K) = RHO(NIM1(I),K) - RHOS1(I,K)
      ENDDO
      DRHOX(I,1) = .25*GRAV*DZ(1)*(D(I)+D(NIM1(I)))
     **(RHOS2(I,1)-RHOS1(I,1))
      DO K=2,KBM1
          DRHOX(I,K) = DRHOX(I,K-1)
     *    +GRAV*.125*(DZ(K)+DZ(K-1))*(D(I)+D(NIM1(I)))
     *    *(RHOS2(I,K)-RHOS1(I,K)+RHOS2(I,K-1)-RHOS1(I,K-1))
     *    +GRAV*.5*Z(K)*(D(I)-D(NIM1(I)))
     *    *(RHOS2(I,K)+RHOS1(I,K)-RHOS2(I,K-1)-RHOS1(I,K-1))
      ENDDO
      DO K=1,KBM1
          DRHOX(I,K) = .25*(D(I)+D(NIM1(I)))*DRHOX(I,K)
     *    *(H2(I)+H2(NIM1(I)))*RAMP
      ENDDO
c --- CBR: RHO remains unchanged
c      DO K=1,KBM1
c          RHO(I,K) = RHO(I,K) + RHOS2(I,K)
c          RHO(NIM1(I),K) = RHO(NIM1(I),K) + RHOS1(I,K)
c      ENDDO
      ENDIF     
      ENDDO
c!$OMP END TARGET
c#ifdef OMP
c!$OMP END PARALLEL DO
c#endif

c#ifdef OMP
c!$OMP PARALLEL DO PRIVATE(RHOZ,RHOS1,RHOS2,DUZ,DUZ1,DUZ2)
c#endif  
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I= 1,N_CTRD
      DO concurrent (I = 1:N_CTRD) local(K,KI,KJ)
     * local_init(RAMP,KBM1)
      IF (DVM(I).GT.0.0) THEN
      DO K=1,KBM1
          DUZ(I,K) = ZZ(K)*0.5*(D(I)+D(NJM1(I)))+0.5*(EL(I)+EL(NJM1(I)))
          DUZ1(I,K) = ZZ(K)*D(NJM1(I))+EL(NJM1(I))
          DUZ2(I,K) = ZZ(K)*D(I)+EL(I)
          RHOZ(I,K) = 0.5 * (RHO(I,K)+RHO(NJM1(I),K))
      ENDDO
c      PURE SUBROUTINE SINTER(X,A,Y,B,M,N) Extracted
c      CALL SINTER(DUZ,RHOZ,DUZ1,RHOS1,KBM1,KBM1) Extracted
c ----- rivision 1:   
      KI=1
      KJ=1
      DO WHILE (KI<=KBM1)
          IF (DUZ1(I,KI)>DUZ(I,1)) THEN 
              RHOS1(I,KI)=RHOZ(I,1)+((RHOZ(I,1)-RHOZ(I,2))/
     *      (DUZ(I,1)-DUZ(I,2)))*(DUZ1(I,KI)-DUZ(I,1))
              KI=KI+1
          ELSEIF (KJ>KBM1) THEN
              RHOS1(I,KI)=RHOZ(I,KBM1)
              KI=KI+1
          ELSEIF (DUZ1(I,KI)>DUZ(I,KJ)) THEN
              RHOS1(I,KI) = RHOZ(I,KJ-1) - (RHOZ(I,KJ-1)-RHOZ(I,KJ))
     *         * (DUZ(I,KJ-1)-DUZ1(I,KI)) / (DUZ(I,KJ-1)-DUZ(I,KJ))
              KI=KI+1
          ELSE
              KJ=KJ+1
          ENDIF
      ENDDO
c      CALL SINTER(DUZ,RHOZ,DUZ2,RHOS2,KBM1,KBM1) Extracted
      KI=1
      KJ=1
      DO WHILE (KI<=KBM1)
          IF (DUZ2(I,KI)>DUZ(I,1)) THEN 
              RHOS2(I,KI)=RHOZ(I,1)+((RHOZ(I,1)-RHOZ(I,2))/
     *      (DUZ(I,1)-DUZ(I,2)))*(DUZ2(I,KI)-DUZ(I,1))
              KI=KI+1
          ELSEIF (KJ>KBM1) THEN
              RHOS2(I,KI)=RHOZ(I,KBM1)
              KI=KI+1
          ELSEIF (DUZ2(I,KI)>DUZ(I,KJ)) THEN
              RHOS2(I,KI) = RHOZ(I,KJ-1) - (RHOZ(I,KJ-1)-RHOZ(I,KJ))
     *         * (DUZ(I,KJ-1)-DUZ2(I,KI)) / (DUZ(I,KJ-1)-DUZ(I,KJ))
              KI=KI+1
          ELSE
              KJ=KJ+1
          ENDIF
      ENDDO
c --- CBR: Using temp variables RHOS1,RHOS2 instead of RHO itself
c      DO K=1,KBM1
c          RHO(I,K) = RHO(I,K) - RHOS2(I,K)
c          RHO(NJM1(I),K) = RHO(NJM1(I),K) - RHOS1(I,K)
c      ENDDO
      DO K=1,KBM1
          RHOS2(I,K) = RHO(I,K) - RHOS2(I,K)
          RHOS1(I,K) = RHO(NJM1(I),K) - RHOS1(I,K)
      ENDDO
      DRHOY(I,1) = .25*GRAV*DZ(1)*(D(I)+D(NJM1(I)))
     **(RHOS2(I,1)-RHOS1(I,1))
      DO K=2,KBM1
          DRHOY(I,K) = DRHOY(I,K-1)
     *    +.125*GRAV*(DZ(K)+DZ(K-1))*(D(I)+D(NJM1(I)))
     *    *(RHOS2(I,K)-RHOS1(I,K)+RHOS2(I,K-1)-RHOS1(I,K-1))
     *    +.5*GRAV*Z(K)*(D(I)-D(NJM1(I)))
     *    *(RHOS2(I,K)+RHOS1(I,K)-RHOS2(I,K-1)-RHOS1(I,K-1))
      ENDDO
      DO K=1,KBM1
          DRHOY(I,K) = .25*(D(I)+D(NJM1(I)))
     *    *DRHOY(I,K)*(H1(I)+H1(NJM1(I)))*RAMP
      ENDDO
c --- CBR: RHO remains unchanged
c      DO K=1,KBM1
c          RHO(I,K) = RHO(I,K) + RHOS2(I,K)
c          RHO(NJM1(I),K) = RHO(NJM1(I),K) + RHOS1(I,K)
c      ENDDO
      ENDIF
      ENDDO
c!$OMP END TARGET
c#ifdef OMP
c!$OMP END PARALLEL DO
c#endif
#if defined TIDE_EL || defined TIDE_FLUX  || defined TIDE_FLATHER      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP PARALLEL DO COLLAPSE(2)
c      DO I=1,NUMEBC
c      DO K=1,KBM1
      DO concurrent (I=1:NUMEBC,K=1:KBM1) local(IC,IE)
      IC = NETA(I)
      IE = NCON(I)
      IF (IC==NIM1(IE)) THEN
          DRHOX(IE,K)=0
      ELSEIF (IC==NIP1(IE)) THEN
          DRHOX(IC,K)=0
      ELSEIF (IC==NJM1(IE)) THEN
          DRHOY(IE,K)=0
      ELSE
          DRHOY(IC,K)=0
      ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET
#endif      
      RETURN
      END SUBROUTINE BAROPG5B
