#include "DEFS.h"	
      SUBROUTINE INIT_CBCCOF
      USE MOD_GLOBAL
      
      INTEGER :: NCOF
      INTEGER :: KK,NBOUND,K,IND
      REAL :: COF1,COF2,COF3,COF4
      REAL, ALLOCATABLE :: XB(:)
      REAL, ALLOCATABLE :: YB(:)      
      REAL, ALLOCATABLE :: XBOUND(:)
      REAL, ALLOCATABLE :: YBOUND(:)
      CHARACTER*160 REGION
#ifdef Z0BTYP_USERDEF
	FN="./"//TRIM(IN_DIRE)//TRIM(XG)//TRIM(CASENAME)//"_lba.dat"  
	OPEN (IULBA,FILE=FN,STATUS='OLD')
#endif
      NCOF = 0
      DO WHILE(.TRUE.)
          READ (IULBA,*,END=300)
          NCOF = NCOF+1
      ENDDO
300   CONTINUE
      REWIND(IULBA)      
      
      READ (IULBA,*)
      READ (IULBA,*) REGION,COF1,COF2,COF3,COF4
      
      CBC_COF_A1 = COF1
      CBC_COF_B1 = COF2
      CBC_COF_A2 = COF3
      CBC_COF_B2 = COF4
      
      DO WHILE(.TRUE.) 
          READ (IULBA,*,END=200) REGION,COF1,COF2,COF3,COF4
          KK = 1
          ALLOCATE (XB(NCOF))
          ALLOCATE (YB(NCOF))
          DO WHILE (.TRUE.) 
              READ (IULBA,*,ERR=100,END=100) XB(KK),YB(KK)
              KK = KK+1
          ENDDO
100       CONTINUE
          BACKSPACE(IULBA)
          NBOUND = KK-1
          ALLOCATE(XBOUND(NBOUND))
          ALLOCATE(YBOUND(NBOUND))
          DO K = 1,NBOUND
              XBOUND(K) = XB(K)
              YBOUND(K) = YB(K)
          ENDDO
          DEALLOCATE (XB)
          DEALLOCATE (YB)
          DO I = 1,N_CTRD
          IF (H(I).GE.-10) THEN
              CALL INSIDE(XR(I),YR(I),XBOUND,YBOUND,NBOUND,IND)
              IF (IND.EQ.1) THEN
                  CBC_COF_A1(I) = COF1
                  CBC_COF_B1(I) = COF2
                  CBC_COF_A2(I) = COF3
                  CBC_COF_B2(I) = COF4
              ENDIF
          ENDIF
          ENDDO
          DEALLOCATE(XBOUND)
          DEALLOCATE(YBOUND)
      ENDDO
200   CONTINUE
      CLOSE(IULBA)
      RETURN
      END SUBROUTINE INIT_CBCCOF
      