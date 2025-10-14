#include "DEFS.h"
	SUBROUTINE SOLVER_SOR_TEST
      USE MOD_GLOBAL 
      USE OMP_LIB
      IMPLICIT NONE
      integer,parameter :: ITMAX=1000
      REAL,parameter :: OMEGA=1.41,TOL=1.E-5
      
      INTEGER I,J,K,IDX_S,IDX_E,ITERR !,ITMAX
      REAL ERR_NRM,ERR_TMP,X_NRM,RELERR !,OMEGA,TOL
!      REAL DPOI(NPOI)  !CBR: scalarized
      REAL DPOII  !CBR: scalar replacing DPOI(I)
!      ITMAX=1000
!      OMEGA=1.41
!      TOL=1.E-5
      
!      DPOI=0.  !no use:
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I=1,NPOI
c          DPOI(I)=0.0
c      ENDDO
c!$OMP END TARGET 

!      if (NStep==1) then
!          open (997,file="IDX_SE.dat")
!          DO I=1,NPOI
!              IDX_S=IAPOI(I)
!              IDX_E=IAPOI(I+1)-1
!              write (997,9001) I,IDX_S,IDX_E,(JAPOI(K),K=IDX_S,IDX_E)
!          ENDDO
!          close (997)          
!      endif
!9001  format(20I8)      
!      pause
      
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I=1,NPOI
      DO concurrent (I=1:NPOI) local(DPOII)
          DPOII=1./APOI(KAPOI(I))
          DO K=IAPOI(I),IAPOI(I+1)-1
              APOI(K)=DPOII*APOI(K)
          ENDDO
          RHSPOI(I)=RHSPOI(I)*DPOII
c          IDX_S=IAPOI(I)
c          IDX_E=IAPOI(I+1)-1
c          DO K=IDX_S,IDX_E
c              IF (JAPOI(K).EQ.I) THEN !don't have to for every step
c                  DPOI(I)=1./APOI(K)
c                  EXIT
c              ENDIF
c          ENDDO
      ENDDO
c!$OMP END TARGET 
      
C --- CBR: Merged:
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I=1,NPOI
c          IDX_S=IAPOI(I)
c          IDX_E=IAPOI(I+1)-1
c          DO K=IDX_S,IDX_E
c              APOI(K)=DPOI(I)*APOI(K)
c          ENDDO
c      ENDDO
c!$OMP END TARGET 

C --- CBR: Merged:
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I=1,NPOI
c          RHSPOI(I)=RHSPOI(I)*DPOI(I)
c      ENDDO
c!$OMP END TARGET 

!$ACC UPDATE HOST (APOI,RHSPOI,XPOI) !Pass to CPU
#ifdef TIDE_EL
      APOI(ANUM)=1. !V2412: moved from before "call SOLVER_SOR_TEST"
      RHSPOI(N_CTRD_AG+1)=0 !V2412: moved from before "call SOLVER_SOR_TEST"
#endif
      
      DO ITERR=1,ITMAX
          
          ERR_NRM=0.
          X_NRM=0. !CBR
          
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO REDUCTION (+:ERR_NRM,X_NRM)
          DO I=1,NPOI
              IDX_S=IAPOI(I)
              IDX_E=IAPOI(I+1)-1
              ERR_TMP=RHSPOI(I)
              DO K=IDX_S,IDX_E
                  J=JAPOI(K)
                  ERR_TMP=ERR_TMP-APOI(K)*XPOI(J)
              ENDDO
              ERR_NRM=ERR_NRM+ERR_TMP*ERR_TMP
              XPOI(I)=XPOI(I)+OMEGA*ERR_TMP
              X_NRM=X_NRM+XPOI(I)**2  !CBR
          ENDDO
c!$OMP END TARGET 
          
!          X_NRM=DOT_PRODUCT(XPOI,XPOI)
!          RELERR=OMEGA*SQRT(ERR_NRM)/SQRT(X_NRM)
          RELERR=OMEGA*SQRT(ERR_NRM/X_NRM)
          
!          print*,ITERR,RELERR,ERR_NRM,X_NRM
          
          IF (RELERR.LT.TOL)  THEN
!$ACC UPDATE DEVICE (XPOI) !Pass to GPU
              RETURN
          ENDIF
      ENDDO
      
      PRINT*,'ERROR IN ELEVATION SOLVER (ITERR=ITMAX)'
      PRINT*,'NSTEP = ',NSTEP
      PAUSE
      STOP
      
      RETURN
	END 