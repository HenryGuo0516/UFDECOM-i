#include "DEFS.h"	
      SUBROUTINE INIT_COF   !!XK
      
      USE MOD_GLOBAL
      INTEGER :: NNCOF
      INTEGER :: KKK,NNBOUND,IK,IIND
      REAL :: NCOFH,NCOFM,COF_KM,COF_KH  !!XK
      REAL :: WEI_VER(N_CTRDP1)
      REAL :: KM_TEMP(N_CTRD)
      REAL :: KH_TEMP(N_CTRD)
      REAL, ALLOCATABLE :: XXB(:)
      REAL, ALLOCATABLE :: YYB(:)      
      REAL, ALLOCATABLE :: XXBOUND(:)
      REAL, ALLOCATABLE :: YYBOUND(:)
      CHARACTER*160 NREGION
      
      FN="./"//TRIM(IN_DIRE)//TRIM(XG)//TRIM(CASENAME)//"_ver.dat"  !XK
	OPEN (IUVER,FILE=FN,STATUS='OLD')
      
      NNCOF = 0
      DO WHILE(.TRUE.)
          READ (IUVER,*,END=301)
          NNCOF = NNCOF+1
      ENDDO
301   CONTINUE
      REWIND(IUVER)      
      READ (IUVER,*)
      READ (IUVER,*) NREGION,NCOFM,NCOFH
      DO I = 1,N_CTRD
	    KM_COF(I) = NCOFM    !!XK
          KH_COF(I) = NCOFH    !!XK
      ENDDO
      DO WHILE(.TRUE.)
          READ (IUVER,*,END=201) NREGION,COF_KM,COF_KH  !!XK
          KKK = 1
          ALLOCATE (XXB(NNCOF))
          ALLOCATE (YYB(NNCOF))
          DO WHILE (.TRUE.) 
              READ (IUVER,*,ERR=101,END=101) XXB(KKK),YYB(KKK)
              KKK = KKK+1
          ENDDO
101       CONTINUE
          BACKSPACE(IUVER)
          NNBOUND = KKK-1
          ALLOCATE (XXBOUND(NNBOUND))
          ALLOCATE (YYBOUND(NNBOUND))
          DO IK = 1,NNBOUND
              XXBOUND(IK) = XXB(IK)
              YYBOUND(IK) = YYB(IK)
          ENDDO
          DEALLOCATE (XXB)
          DEALLOCATE (YYB) 
          DO I = 1,N_CTRD
          IF (H(I).GE.-10) THEN
        	CALL INSIDE(XR(I),YR(I),XXBOUND,YYBOUND,NNBOUND,IIND)
        	IF (IIND.EQ.1) THEN
	  	    KM_COF(I) = COF_KM  !!XK
	  	    KH_COF(I) = COF_KH  !!XK
        	ENDIF
          ENDIF
          ENDDO
          DEALLOCATE(XXBOUND)
          DEALLOCATE(YYBOUND)
      ENDDO
201   CONTINUE
      REWIND(IUVER)
      CLOSE(IUVER)
      
********************************* SMOOTH *********************************
      WEI_VER = 0.
      DO I = 1,N_CTRD
          IF (H(I).GE.-10) WEI_VER(I) = 1.
      ENDDO
      
      KM_TEMP = 0.
      KH_TEMP = 0.
 !     DO M = 1,7
 !     DO I = 1,N_CTRD_AG
 !     IF (H(I).GE. -10) THEN 
 !         KM_TEMP(I) = (KM_COF(I)*WEI_VER(I)   !XK
 !    *    +KM_COF(NIP1(I))*WEI_VER(NIP1(I))
 !    *    +KM_COF(NIM1(I))*WEI_VER(NIM1(I))
 !    *    +KM_COF(NJM1(I))*WEI_VER(NJM1(I))
 !    *    +KM_COF(NJP1(I))*WEI_VER(NJP1(I)))
 !    *    /(WEI_VER(I)+WEI_VER(NIP1(I))
 !    *    +WEI_VER(NIM1(I))+WEI_VER(NJP1(I))
 !    *    +WEI_VER(NJM1(I)))
 !
 !         KH_TEMP(I) = (KH_COF(I)*WEI_VER(I)    !XK
 !    *    +KH_COF(NIP1(I))*WEI_VER(NIP1(I))
 !    *    +KH_COF(NIM1(I))*WEI_VER(NIM1(I))
 !    *    +KH_COF(NJM1(I))*WEI_VER(NJM1(I))
 !    *    +KH_COF(NJP1(I))*WEI_VER(NJP1(I)))
 !    *    /(WEI_VER(I)+WEI_VER(NIP1(I))
 !    *    +WEI_VER(NIM1(I))+WEI_VER(NJP1(I))
 !    *    +WEI_VER(NJM1(I)))
 !     ENDIF
 !     ENDDO
 !     
	!VTP = 0
 !     
 !     CALL INFO_EXCH(VTP,1,KM_TEMP,FSMADD(1:N_CTRD)
 !    * ,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
 !     CALL INFO_EXCH(VTP,1,KH_TEMP,FSMADD(1:N_CTRD)
 !    * ,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
 !     DO I=1,N_CTRD
 !     IF (H(I) .GE. -10) THEN
 !         KM_COF(I) = KM_TEMP(I)  !XK
 !         KH_COF(I) = KH_TEMP(I)  !XK
 !     ENDIF
 !     ENDDO
 !     ENDDO
      VTP = 0
      
      CALL INFO_EXCH(VTP,1,KM_COF,FSMADD(1:N_CTRD)
     * ,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
      CALL INFO_EXCH(VTP,1,KH_COF,FSMADD(1:N_CTRD)
     * ,NAGC_EVG,WTC_EVG,NAGC_IVG,WTC_IVG)
1211  FORMAT(2I5,2F7.2) 
      RETURN
      END SUBROUTINE INIT_COF