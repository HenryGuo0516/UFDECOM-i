#INCLUDE "DEFS.h"
	
	SUBROUTINE BCOND_TWO(IDX)
      
      USE MOD_GLOBAL
      
      INTEGER, INTENT(IN) :: IDX
      
      SELECT CASE (IDX)  
          
      CASE(1)
!======================================================================   
!                     MATERIAL   OBC
!======================================================================
#ifdef MODULE_MATERIAL
      DO N = 1,NUMQBC
      ID = NQD(N)
      IC = NQC(N)
      DO K = 1,KBM1
      IF (VQDIST(N,K).NE.0.0.AND.QDIS(N).GE.0.0) THEN     
      IF (ID.EQ.NIM1(IC)) THEN  ! EAST
          CONCENTRATION(IC,K)=CONCENTRATION(ID,K)
      ELSEIF (ID.EQ.NIP1(IC)) THEN  !WEST
          CONCENTRATION(IC,K)=MDIS(N)
      ELSEIF (ID.EQ.NJM1(IC)) THEN  !NORTH
          CONCENTRATION(IC,K)=CONCENTRATION(ID,K)
      ELSEIF (ID.EQ.NJP1(IC)) THEN  !SORTH
          CONCENTRATION(IC,K)=MDIS(N)
      ENDIF
      ELSEIF (VQDIST(N,K).NE.0.0.AND.QDIS(N).LT.0.0) THEN
      IF (ID.EQ.NIM1(IC)) THEN  ! EAST
          CONCENTRATION(IC,K)=MDIS(N)
      ELSEIF (ID.EQ.NIP1(IC)) THEN  !WEST
          CONCENTRATION(IC,K)=CONCENTRATION(ID,K)
      ELSEIF (ID.EQ.NJM1(IC)) THEN  !NORTH
          CONCENTRATION(IC,K)=MDIS(N)
      ELSEIF (ID.EQ.NJP1(IC)) THEN  !SORTH
          CONCENTRATION(IC,K)=CONCENTRATION(ID,K)
      ENDIF
      ENDIF
      ENDDO
      ENDDO
#endif    
      END SELECT
      RETURN
      END  SUBROUTINE BCOND_TWO