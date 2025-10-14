#include "DEFS.h"
	
      SUBROUTINE CAL_SED
	
	USE MOD_GLOBAL
	USE ADVT
#ifdef MODULE_SED
      USE MOD_SED
#endif
#ifdef MODULE_WAVE
      USE MOD_WAVE
#endif
	
	IMPLICIT NONE
	
	INTEGER NCON2
	INTEGER I,J,K,N,II
	
#ifdef MODULE_SED
	CALL SEDW   !³Á½µËÙ¶È 
#endif
	
#ifdef HSIMT
	NCON2 = 2
      CALL ADVT_HSIMT(SED,WFF)  !WU HUI
#else
      PRINT*, 'Advection scheme undefined!'
	PAUSE
	STOP
#endif

#ifdef MODULE_SED
	CALL RESUSPEN(WFF)  ! ÆðÐü
#endif

#ifdef BED
	!IF (LOG_BED) THEN
        IF(THOUR .GE. BED_BEG)THEN
          CALL RESEABED
	  ENDIF	    
	!ENDIF  
#endif
     
      CALL STEP2_SED2(WFF)
	
      CALL BCOND(7)		 

	DO K = 1, KB
        DO I = 1, N_CTRD_AG
	    IF (FSM(I) .GE. 0.5) THEN                !WUHUI
	      SED(I,K) = WFF(I,K)
	    ENDIF        
        ENDDO
	ENDDO   
		
	DO N = 1, NUMEBC
	  II = NETA(N)
	  DO K = 1,KB
		SED(II,K) = WFF(II,K)
	  ENDDO
	ENDDO
		
	DO N = 1, NUMQBC
        II = NQC(N)
        DO K = 1,KB
          SED(II,K) = WFF(II,K)
        ENDDO
      ENDDO

	RETURN
	
	END SUBROUTINE CAL_SED