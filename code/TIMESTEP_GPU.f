#include "DEFS.h"      
      SUBROUTINE TIMESTEP_GPU
! --- CBR: For "defined IMP_EL" only  
      
	USE MOD_GLOBAL
      PARAMETER (UVSCALE = 1.0)
      
      
      INTEGER N(4)
      REAL DT(4)
	REAL DT1,DT2
      DTIMAX = DTIMAX1
 
      
c      UVSCALE = 1.0  !CBR: Moved to Parameter
      
c	N = 0
c	N_DTI = 0
c	DT = 1000000000
      DTT = 1000000000.
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO reduction(MIN:DTT)
c      DO I = 1,N_CTRD_AG
      DO concurrent (I = 1:N_CTRD_AG) reduce(MIN:DTT)
      IF (DUM(I).EQ.1..AND.((U(I,1).NE.0.).OR.(V(I,1).NE.0.))) THEN
          DTT = MIN(
     *    H1(I)/SQRT(U(I,1)**2+(UVSCALE*V(I,1))**2),
     *    H2(I)/SQRT(V(I,1)**2+(UVSCALE*U(I,1))**2) )
      ENDIF
      ENDDO
c!$OMP END TARGET

      DTI=MIN(FACTOR*DTT,DTIMAX)
      
c	IF (DTI.GE.DTIMAX) DTI = DTIMAX
c      IF (N_DTI .GE. 1) THEN
      IF (DTT < 900000000.) THEN

      ELSE
          IF (NSTEP.GT.500) THEN 
              PRINT*, 'NSTEP = ',NSTEP
              PRINT*, 'CFL criterion failed!'
#ifdef MODEL_WARNING
              CALL SYSTEM('python model_warning.py')
#endif
              PAUSE
              STOP
          ENDIF
      ENDIF
	IF (DTI .LE. DTIMAX/100.) THEN
          WRITE (*,*) 'TIME STEP IS ABNORMAL, CHECKING AGAIN...'
!$ACC UPDATE HOST (DUM,DVM,U,V,DU,DV)
          CALL TIMESTEP
      ENDIF
      
	RETURN
	END SUBROUTINE TIMESTEP_GPU
