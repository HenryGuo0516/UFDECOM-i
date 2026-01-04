#include "DEFS.h"      
      SUBROUTINE TIMESTEP
	USE MOD_GLOBAL
      PARAMETER (UVSCALE = 1.0)
      
      
      INTEGER N(4)
      REAL DT(4)
	REAL DT1,DT2
      DTIMAX = DTIMAX1
 
      
c      UVSCALE = 1.0  !CBR: Moved to Parameter
	N = 0
	N_DTI = 0
	DT = 1000000000
      DO I = 1,N_CTRD_AG
      IF (DUM(I).EQ.1..AND.((U(I,1).NE.0.).OR.(V(I,1).NE.0.))) THEN
	    DTT = FACTOR*ABS(H1(I)/
     *    SQRT(U(I,1)**2+(UVSCALE*V(I,1))**2))
          IF (DTT.LE.DT(1))THEN
              DT(1) = DTT
              N(1) = I
          ENDIF
          DTT = FACTOR*ABS(H2(I)/
     *    SQRT(V(I,1)**2+(UVSCALE*U(I,1))**2))
          IF (DTT.LE.DT(2))THEN
	        DT(2) = DTT
	        N(2) = I
          ENDIF
      ENDIF
      
#ifdef EXP_EL
      IF (DUM(I).EQ.1.) THEN
      DTT = FACTOR*ABS(H1(I)/SQRT(GRAV*DU(I)))
      IF (DTT.LE.DT(3))THEN
	      DT(3) = DTT
	      N(3) = I
      ENDIF
      ENDIF
      IF (DVM(I).EQ.1.) THEN
      DTT = FACTOR*ABS(H2(I)/SQRT(GRAV*DV(I)))
      IF (DTT.LE.DT(4))THEN
          DT(4) = DTT
          N(4) = I
      ENDIF
      ENDIF
#endif  
      ENDDO
#ifdef EXP_EL
      DTI = MINVAL(DT)
      N_DTI = N(MINLOC(DT,1))
#elif defined IMP_EL
      IF (DT(1).GT.DT(2)) THEN
          DTI=DT(2)
          N_DTI=N(2)
      ELSE
          DTI=DT(1)
          N_DTI=N(1)
      ENDIF
#endif  
	IF (DTI.GE.DTIMAX) DTI = DTIMAX
      IF (N_DTI .GE. 1) THEN
          
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
      WRITE (*,*) 'TIME STEP IS ABNORMAL, MAYBE IT HAS BEEN OVERFLOW!'
      PRINT*,N_DTI,U(N_DIT,1),V(N_DTI,1)
#ifdef MODEL_WARNING
          !CALL SYSTEM('python model_warning.py')
#endif
	    PAUSE
          STOP
      ENDIF
      
	RETURN
	END SUBROUTINE TIMESTEP
