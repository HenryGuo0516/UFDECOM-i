#include "DEFS.h"

c  CBR: NCON1==1 for T; Else for S, Sed, etc. 
      SUBROUTINE PROFT_GPU(NCON1,F,WFSURF)
C     VERSION(03/02/90)
C
c      USE MOD_GLOBAL,AF=>A,FF=>A,CF=>PROD,VHF=>VH,VHPF=>VHP
      USE MOD_GLOBAL,AF=>A,CF=>PROD,VHF=>VH,VHPF=>VHP
      Parameter (UMOLPR = 0.000001)

c      REAL AF(N_CTRD_AG,KB), CF(N_CTRD_AG,KB)
c      REAL VHF(N_CTRD_AG,KB), VHPF(N_CTRD_AG,KB)
c      REAL FF(N_CTRD_AG,KB), RAD(N_CTRD_AG,KB)
      DIMENSION F(N_CTRDP1,KB), WFSURF(N_CTRD)

c      UMOLPR = 0.000001
      
C**********************************************************************
C                                                                     *
C        THE FOLLOWING SECTION SOLVES THE EQUATION                    *
C         DTI*(KH*F')'-F=-FB                                          *
C                                                                     *
C**********************************************************************
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c       DO K = 2,KBM1
c         DO I = 1,N_CTRD_AG
      DO concurrent (I=1:N_CTRD_AG,K=2:KBM1) local_init(DTI) 
           IF(FSM(I)*FSM11(I).EQ.1.) THEN
          AF(I,K-1) = -DTI*(KH(I,K)+UMOLPR)/(DZ(K-1)*DZZ(K-1)*D(I)*D(I))
          CF(I,K) = AF(I,K-1)*(DZ(K-1)/DZ(K))
c          CF(I,K) = -DTI*(KH(I,K)+UMOLPR)/(DZ(K)*DZZ(K-1)*D(I)*D(I))
           ENDIF
         ENDDO
c       ENDDO
c!$OMP END TARGET
       
      IF (NCON1==1) THEN !CBR: For T
C
C     THE NET HEAT FLUX INPUT
C
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c         DO K = 1,KB
c           DO I = 1,N_CTRD_AG
      DO concurrent (I=1:N_CTRD_AG,K=1:KB) local(ZDEP) 
             IF (FSM(I)*FSM11(I).EQ.1.) THEN
               ZDEP = Z(K)*D(I)
               RAD(I,K) = SWRAD(I)
     * *(RHEAT*EXP(ZDEP/ZETA1)+(1.-RHEAT)*EXP(ZDEP/ZETA2))
             ENDIF
           ENDDO
c         ENDDO
c!$OMP END TARGET

C----- SURFACE BC'S; WFSURF ----------------------------------
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I = 1,N_CTRD_AG
      DO concurrent (I=1:N_CTRD_AG) local(K) local_init(KBM2,DTI) 
        IF(FSM(I)*FSM11(I).EQ.1.) THEN
          VHF(I,1) = AF(I,1) / (AF(I,1)-1.)
          VHPF(I,1) = -DTI*(WFSURF(I)-SWRAD(I)+RAD(I,1)-RAD(I,2))
     * /(-DZ(1)*D(I))-F(I,1)
          VHPF(I,1) = VHPF(I,1)/(AF(I,1)-1.)
          DO K = 2,KBM2
            VHPF(I,K) = 1./(AF(I,K)+CF(I,K)*(1.-VHF(I,K-1))-1.)
            VHF(I,K) = AF(I,K)*VHPF(I,K)
            VHPF(I,K) = (CF(I,K)*VHPF(I,K-1)-F(I,K)+DTI*(RAD(I,K)
     * -RAD(I,K+1))/(D(I)*DZ(K)))*VHPF(I,K)
          ENDDO
        ENDIF
      ENDDO
c!$OMP END TARGET


      ELSE !CBR: For S, Sed, etc.

C----- SURFACE BC'S; WFSURF ----------------------------------
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c       DO I = 1,N_CTRD_AG
      DO concurrent (I=1:N_CTRD_AG) local(K) local_init(KBM2,DTI) 
       IF(FSM(I)*FSM11(I).EQ.1.) THEN
          VHF(I,1) = AF(I,1) / (AF(I,1)-1.)
          VHPF(I,1) = -DTI*WFSURF(I)
     *            /(-DZ(1)*D(I))-F(I,1)
          VHPF(I,1) = VHPF(I,1)/(AF(I,1)-1.)
          DO K = 2,KBM2
            VHPF(I,K) = 1./(AF(I,K)+CF(I,K)*(1.-VHF(I,K-1))-1.)
            VHF(I,K) = AF(I,K)*VHPF(I,K)
            VHPF(I,K) = (CF(I,K)*VHPF(I,K-1)-F(I,K))*VHPF(I,K)
          ENDDO
       ENDIF
       ENDDO
c!$OMP END TARGET

      ENDIF !CBR: NCON1

       
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c       DO K = 1,KBM1
c         DO I = 1,N_CTRD_AG
c           IF(FSM(I)*FSM11(I).EQ.1.) THEN
c              FF(I,K) = F(I,K)
c           ENDIF
c         ENDDO
c       ENDDO
c!$OMP END TARGET

       
      IF (NCON1==1) THEN !CBR: For T
       
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c       DO I = 1,N_CTRD_AG
      DO concurrent (I=1:N_CTRD_AG) local(K) 
     * local_init(KB,KBM1,KBM2,DTI) 
         IF (FSM(I)*FSM11(I).EQ.1.) THEN
           F(I,KBM1) = ((CF(I,KBM1)*VHPF(I,KBM2)-F(I,KBM1)
     * +DTI*(RAD(I,KBM1)-RAD(I,KB))/(D(I)*DZ(KBM1)))
     * /(CF(I,KBM1)*(1.-VHF(I,KBM2))-1.))
c          DO K = 2,KBM1
c            KI = KB-K
c            F(I,KI) = (VHF(I,KI)*F(I,KI+1)+VHPF(I,KI))
c          ENDDO
           DO K = KBM2,1,-1
             F(I,K) = (VHF(I,K)*F(I,K+1)+VHPF(I,K))
           ENDDO
         ENDIF
       ENDDO
c!$OMP END TARGET

      ELSE !CBR: For S, Sed, etc.

c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c       DO I = 1,N_CTRD_AG
      DO concurrent (I=1:N_CTRD_AG) local(K) 
     * local_init(KBM1,KBM2,DTI) 
         IF (FSM(I)*FSM11(I).EQ.1.) THEN
           F(I,KBM1) = (CF(I,KBM1)*VHPF(I,KBM2)-F(I,KBM1))
     * /(CF(I,KBM1)*(1.-VHF(I,KBM2))-1.)
c          DO K = 2,KBM1
c            KI = KB-K
c            F(I,KI) = (VHF(I,KI)*F(I,KI+1)+VHPF(I,KI))
c          ENDDO
           DO K = KBM2,1,-1
             F(I,K) = (VHF(I,K)*F(I,K+1)+VHPF(I,K))
           ENDDO
         ENDIF
       ENDDO
c!$OMP END TARGET
          
      ENDIF !CBR: NCON1


c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 1,KBM1
c        DO I = 1,N_CTRD_AG
c          IF(FSM(I)*FSM11(I).EQ.1.) THEN
c            F(I,K) = FF(I,K)
c          ENDIF
c        ENDDO
c      ENDDO
c      ENDDO
c!$OMP END TARGET
      
      RETURN
      
      END SUBROUTINE PROFT_GPU

c ---------------------------------------------------------------------
c ---------------------------------------------------------------------
c ---------------------------------------------------------------------
c ---------------------------------------------------------------------
c ---------------------------------------------------------------------
c ---------------------------------------------------------------------

c  CBR: NCON1==1 for T; Else for S, Sed, etc. 
      SUBROUTINE PROFT(NCON1,F,WFSURF)
C     VERSION(03/02/90)
C
      USE MOD_GLOBAL,AF=>A,FF=>A,CF=>PROD,VHF=>VH,VHPF=>VHP
      Parameter (UMOLPR = 0.000001)

c      REAL AF(N_CTRD_AG,KB), CF(N_CTRD_AG,KB)
c      REAL VHF(N_CTRD_AG,KB), VHPF(N_CTRD_AG,KB)
c      REAL FF(N_CTRD_AG,KB), RAD(N_CTRD_AG,KB)
      DIMENSION F(N_CTRDP1,KB), WFSURF(N_CTRD)

c      UMOLPR = 0.000001
      
C**********************************************************************
C                                                                     *
C        THE FOLLOWING SECTION SOLVES THE EQUATION                    *
C         DTI*(KH*F')'-F=-FB                                          *
C                                                                     *
C**********************************************************************
       DO K = 2,KBM1
#ifdef OMP
!$OMP PARALLEL DO
#endif
         DO I = 1,N_CTRD_AG
           IF(FSM(I)*FSM11(I).EQ.1.) THEN
          AF(I,K-1) = -DTI*(KH(I,K)+UMOLPR)/(DZ(K-1)*DZZ(K-1)*D(I)*D(I))
          CF(I,K) = -DTI*(KH(I,K)+UMOLPR)/(DZ(K)*DZZ(K-1)*D(I)*D(I))
           ENDIF
         ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif
       ENDDO

       
      IF (NCON1==1) THEN !CBR: For T
C
C     THE NET HEAT FLUX INPUT
C
         DO K = 1,KB
         !DO K=1,KBM1
#ifdef OMP
!$OMP PARALLEL DO PRIVATE(ZDEP)
#endif
           DO I = 1,N_CTRD_AG
             IF (FSM(I)*FSM11(I).EQ.1.) THEN
               ZDEP = Z(K)*D(I)
               RAD(I,K) = SWRAD(I)
     * *(RHEAT*EXP(ZDEP/ZETA1)+(1-RHEAT)*EXP(ZDEP/ZETA2))
             ENDIF
           ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif
         ENDDO


C----- SURFACE BC'S; WFSURF ----------------------------------
#ifdef OMP
!$OMP PARALLEL DO
#endif
      DO I = 1,N_CTRD_AG
        IF(FSM(I)*FSM11(I).EQ.1.) THEN
          VHF(I,1) = AF(I,1) / (AF(I,1)-1.)
          VHPF(I,1) = -DTI*(WFSURF(I)-SWRAD(I)+RAD(I,1)-RAD(I,2))
     * /(-DZ(1)*D(I))-F(I,1)
          VHPF(I,1) = VHPF(I,1)/(AF(I,1)-1.)
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
          IF(FSM(I)*FSM11(I).EQ.1.) THEN
            VHPF(I,K) = 1./(AF(I,K)+CF(I,K)*(1.-VHF(I,K-1))-1.)
            VHF(I,K) = AF(I,K)*VHPF(I,K)
c            VHPF(I,K) = (CF(I,K)*VHPF(I,K-1)-DBLE(F(I,K))+DTI*(RAD(I,K)
            VHPF(I,K) = (CF(I,K)*VHPF(I,K-1)-F(I,K)+DTI*(RAD(I,K)
     * -RAD(I,K+1))/(D(I)*DZ(K)))*VHPF(I,K)
          ENDIF
        ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif
      ENDDO

      ELSE !CBR: For S, Sed, etc.

C----- SURFACE BC'S; WFSURF ----------------------------------
          DO I = 1,N_CTRD_AG
            IF(FSM(I)*FSM11(I).EQ.1.) THEN
              VHF(I,1) = AF(I,1) / (AF(I,1)-1.)
              VHPF(I,1) = -DTI*WFSURF(I)
     *                     /(-DZ(1)*D(I))-F(I,1)
              VHPF(I,1) = VHPF(I,1)/(AF(I,1)-1.)
            ENDIF
          ENDDO

        DO K = 2,KBM2
        DO I = 1,N_CTRD_AG
          IF(FSM(I)*FSM11(I).EQ.1.) THEN
            VHPF(I,K) = 1./(AF(I,K)+CF(I,K)*(1.-VHF(I,K-1))-1.)
            VHF(I,K) = AF(I,K)*VHPF(I,K)
c            VHPF(I,K) = (CF(I,K)*VHPF(I,K-1)-DBLE(F(I,K)))*VHPF(I,K)
            VHPF(I,K) = (CF(I,K)*VHPF(I,K-1)-F(I,K))*VHPF(I,K)
          ENDIF
        ENDDO
        ENDDO
          
      ENDIF !CBR: NCON1

       
       DO K = 1,KBM1
         DO I = 1,N_CTRD_AG
           IF(FSM(I)*FSM11(I).EQ.1.) THEN
              FF(I,K) = F(I,K)
           ENDIF
         ENDDO
       ENDDO

       
      IF (NCON1==1) THEN !CBR: For T
       
#ifdef OMP
!$OMP PARALLEL DO
#endif
       DO I = 1,N_CTRD_AG
         IF (FSM(I)*FSM11(I).EQ.1.) THEN
           FF(I,KBM1) = ((CF(I,KBM1)*VHPF(I,KBM2)-FF(I,KBM1)
     * +DTI*(RAD(I,KBM1)-RAD(I,KB))/(D(I)*DZ(KBM1)))
     * /(CF(I,KBM1)*(1.-VHF(I,KBM2))-1.))
         ENDIF
       ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif

      ELSE !CBR: For S, Sed, etc.

       DO I = 1,N_CTRD_AG
         IF (FSM(I)*FSM11(I).EQ.1.) THEN
           FF(I,KBM1) = (CF(I,KBM1)*VHPF(I,KBM2)-FF(I,KBM1))
     * /(CF(I,KBM1)*(1.-VHF(I,KBM2))-1.)
         ENDIF
       ENDDO
          
          
      ENDIF !CBR: NCON1


      DO K = 2,KBM1
        KI = KB-K
#ifdef OMP
!$OMP PARALLEL DO
#endif
        DO I = 1,N_CTRD_AG
          IF(FSM(I)*FSM11(I).EQ.1.) THEN
            FF(I,KI) = (VHF(I,KI)*FF(I,KI+1)+VHPF(I,KI))
          ENDIF
        ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif
      ENDDO

      DO K = 1,KBM1
        DO I = 1,N_CTRD_AG
          IF(FSM(I)*FSM11(I).EQ.1.) THEN
            F(I,K) = FF(I,K)
          ENDIF
        ENDDO
      ENDDO

      
      RETURN
      
      END SUBROUTINE PROFT
