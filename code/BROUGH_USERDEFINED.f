#include "DEFS.h"		
      SUBROUTINE BROUGH_USERDEFINED_GPU
	USE MOD_GLOBAL
      REAL Z0,CD_TEST,VEL1,VEL2
      
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I = 1,N_CTRD
      DO concurrent (I=1:N_CTRD) local(VEL1,VEL2,CD_TEST)
     * local_init(BFRIC,Z0B)
          VEL1=AMIN1(0.64,SQRT(U(I,KBM1)**2+V(I,KBM1)**2))  
          VEL2=SQRT(WINDU(I)**2+WINDV(I)**2)
          VEL2=AMIN1(10.,VEL2)
          VEL2=AMAX1(5.,VEL2)
          VEL2=(VEL2-5.)**2
          
          CD_TEST=(CBC_COF_A1(I)+CBC_COF_B1(I)*(VEL1/0.64))
     *    *(CBC_COF_A2(I)+CBC_COF_B2(I)*(VEL2/25.))
      
      IF(DU(I).GE.3.0) THEN
          CBC_U(I) = AMAX1(BFRIC,
     *    .16/ALOG((ZZ(KBM1)-Z(KB))*DU(I)/Z0B)**2)*CD_TEST
      ELSEIF(DU(I).GE.1..AND.DU(I).LT.3.0) THEN
          CBC_U(I) = AMAX1(BFRIC,
     *    .16/ALOG((ZZ(KBM1)-Z(KB))*3.0/Z0B)**2)*CD_TEST*1.3
      ELSEIF(HU(I).GE.-10..AND.DU(I).LT.1.0) THEN
          CBC_U(I) = AMAX1(BFRIC,
     *    .16/ALOG((ZZ(KBM1)-Z(KB))*1.0/Z0B)**2)
          CBC_U(I) = CBC_U(I)*2
      ENDIF
      
      IF(DV(I).GE.3.0) THEN
          CBC_V(I) = AMAX1(BFRIC,
     *    .16/ALOG((ZZ(KBM1)-Z(KB))*DV(I)/Z0B)**2)*CD_TEST
      ELSEIF(DV(I).GE.1..AND.DV(I).LT.3.0) THEN
          CBC_V(I) = AMAX1(BFRIC,
     *    .16/ALOG((ZZ(KBM1)-Z(KB))*3.0/Z0B)**2)*CD_TEST*1.3
      ELSEIF(HV(I).GE.-10..AND.DV(I).LT.1.0) THEN
          CBC_V(I) = AMAX1(BFRIC,
     *    .16/ALOG((ZZ(KBM1)-Z(KB))*1.0/Z0B)**2)
          CBC_V(I) = CBC_V(I)*2
      ENDIF
#ifdef MODULE_SED
      IF (D(I).GE.3.0) THEN  
          CBC_UV(I) = AMAX1(BFRIC,.16/
     *    ALOG((ZZ(KBM1)-Z(KB))*D(I)/Z0B)**2)*CD_TEST
      ELSEIF (D(I).GE.1..AND.D(I).LT.3.0) THEN
          CBC_UV(I) = AMAX1(BFRIC,.16/                 
     *    ALOG((ZZ(KBM1)-Z(KB))*3.0/Z0B)**2)*CD_TEST*1.3
      ELSEIF (H(I).GE.-10..AND.D(I).LT.1.0) THEN  
          CBC_UV(I) = AMAX1(BFRIC,.16/                 
     *    ALOG((ZZ(KBM1)-Z(KB))*1.0/Z0B)**2)*2    
      ENDIF
#endif
      ENDDO
c!$OMP END TARGET      


c	IF(LOG_CBCADJ)THEN
c      DO I=1,NUMCBCADJ
c      II=CBCADJ(I)
c      IF (FSM(II).GT.0.) THEN
     !! IF(DU(II).LE.2.0.AND.DU(II).GT.0.0) THEN
     !!     CBCMIN = 0.00001
     !!     CBC_U(II) = AMAX1(CBCMIN,(9.8*(0.055)**2.
     !!*    /(DU(II)**(1./3.))))*CBC_COF(II)  
     !! ELSEIF(DU(II).LE.3.0.AND.DU(II).GT.2.0) THEN
     !!     CBCMIN = 0.00001
     !!     CBC_U(II) = AMAX1(CBCMIN,(9.8*(0.040)**2.
     !!*    /(DU(II)**(1./3.))))*CBC_COF(II)  
     !! ELSEIF(DU(II).LE.4.AND.DU(II).GT.3.0) THEN
     !!     CBCMIN = 0.00001
     !!     CBC_U(II) = AMAX1(CBCMIN,(9.8*(0.032)**2.
     !!*    /(DU(II)**(1./3.))))*CBC_COF(II)  
     !! ELSEIF(DU(II).LE.7.0.AND.DU(II).GT.4)THEN    
     !!     CBCMIN = 0.00001
     !!     CBC_U(II) = AMAX1(CBCMIN,(9.8*(0.027)**2.
     !!*    /(DU(II)**(1./3.))))*CBC_COF(II)
     !! ELSEIF(DU(II).LE.10.0.AND.DU(II).GT.7.0)THEN    
     !!     CBCMIN = 0.00001
     !!     CBC_U(II) = AMAX1(CBCMIN,(9.8*(0.014)**2.
     !!*    /(DU(II)**(1./3.))))*CBC_COF(II)
     !! ELSEIF(DU(II).GT.10.0)THEN    
     !!     CBCMIN = 0.00001
     !!     CBC_U(II) = AMAX1(CBCMIN,(9.8*(0.01)**2.
     !!*/(DU(II)**(1./3.))))*CBC_COF(II)
     !! ENDIF
     !! 
     !! IF(DV(II).LE.2.0.AND.DV(II).GT.0.0) THEN
     !!     CBCMIN = 0.00001
     !!     CBC_V(II) = AMAX1(CBCMIN,(9.8*(0.055)**2.
     !!*    /(DV(II)**(1./3.))))*CBC_COF(II)  
     !! ELSEIF(DV(II).LE.3.0.AND.DV(II).GT.2.0) THEN
     !!     CBCMIN = 0.00001
     !!     CBC_V(II) = AMAX1(CBCMIN,(9.8*(0.040)**2.
     !!*    /(DV(II)**(1./3.))))*CBC_COF(II)  
     !! ELSEIF(DV(II).LE.4.AND.DV(II).GT.3.0) THEN
     !!     CBCMIN = 0.00001
     !!     CBC_V(II) = AMAX1(CBCMIN,(9.8*(0.032)**2.
     !!*    /(DV(II)**(1./3.))))*CBC_COF(II)  
     !! ELSEIF(DV(II).LE.7.0.AND.DV(II).GT.4)THEN    
     !!     CBCMIN = 0.00001
     !!     CBC_V(II) = AMAX1(CBCMIN,(9.8*(0.027)**2.
     !!*    /(DV(II)**(1./3.))))*CBC_COF(II)
     !! ELSEIF(DV(II).LE.10.0.AND.DV(II).GT.7.0)THEN    
     !!     CBCMIN = 0.00001
     !!     CBC_V(II) = AMAX1(CBCMIN,(9.8*(0.014)**2.
     !!*    /(DV(II)**(1./3.))))*CBC_COF(II)
     !! ELSEIF(DV(II).GT.10.0)THEN    
     !!     CBCMIN = 0.00001
     !!     CBC_V(II) = AMAX1(CBCMIN,(9.8*(0.01)**2.
     !!*    /(DV(II)**(1./3.))))*CBC_COF(II)
     !! ENDIF
     !! IF(D(II).LE.2.0.AND.D(II).GT.0.0) THEN
     !!     CBCMIN = 0.00001
     !!     CBC_UV(II) = AMAX1(CBCMIN,(9.8*(0.055)**2.
     !!* /(D(II)**(1./3.))))*CBC_COF(II)  
     !! ELSEIF(D(II).LE.3.0.AND.D(II).GT.2.0) THEN
     !!     CBCMIN = 0.00001
     !!     CBC_UV(II) = AMAX1(CBCMIN,(9.8*(0.040)**2.
     !!*    /(D(II)**(1./3.))))*CBC_COF(II)  
     !! ELSEIF(D(II).LE.4.AND.D(II).GT.3.0) THEN
     !!     CBCMIN = 0.00001
     !!     CBC_UV(II) = AMAX1(CBCMIN,(9.8*(0.032)**2.
     !!*    /(D(II)**(1./3.))))*CBC_COF(II)  
     !! ELSEIF(D(II).LE.7.0.AND.D(II).GT.6)THEN    
     !!     CBCMIN = 0.00001
     !!     CBC_UV(II) = AMAX1(CBCMIN,(9.8*(0.027)**2.
     !!*    /(D(II)**(1./3.))))*CBC_COF(II)
     !! ELSEIF(D(II).LE.10.0.AND.D(II).GT.7.0)THEN    
     !!     CBCMIN = 0.00001
     !!     CBC_UV(II) = AMAX1(CBCMIN,(9.8*(0.014)**2.
     !!*    /(D(II)**(1./3.))))*CBC_COF(II)
     !! ELSEIF(D(II).GT.10.0)THEN    
     !!     CBCMIN = 0.00001
     !!     CBC_UV(II) = AMAX1(CBCMIN,(9.8*(0.01)**2.
     !!*    /(D(II)**(1./3.))))*CBC_COF(II)
     !! ENDIF
c      ENDIF
c      ENDDO
c      ENDIF
      
   	RETURN
      END SUBROUTINE BROUGH_USERDEFINED_GPU

      
C  ========================================================= CPU version:      
      SUBROUTINE BROUGH_USERDEFINED
	USE MOD_GLOBAL
      REAL Z0,CD_TEST,VEL1,VEL2
      
      
      DO I = 1,N_CTRD
          VEL1=AMIN1(0.64,SQRT(U(I,KBM1)**2+V(I,KBM1)**2))  
          VEL2=SQRT(WINDU(I)**2+WINDV(I)**2)
          VEL2=AMIN1(10.,VEL2)
          VEL2=AMAX1(5.,VEL2)
          VEL2=(VEL2-5.)**2
          
          CD_TEST=(CBC_COF_A1(I)+CBC_COF_B1(I)*(VEL1/0.64))
     *    *(CBC_COF_A2(I)+CBC_COF_B2(I)*(VEL2/25.))
      
      IF(DU(I).GE.3.0) THEN
          CBC_U(I) = AMAX1(BFRIC,
     *    .16/ALOG((ZZ(KBM1)-Z(KB))*DU(I)/Z0B)**2)*CD_TEST
      ELSEIF(DU(I).GE.1..AND.DU(I).LT.3.0) THEN
          CBC_U(I) = AMAX1(BFRIC,
     *    .16/ALOG((ZZ(KBM1)-Z(KB))*3.0/Z0B)**2)*CD_TEST*1.3
      ELSEIF(HU(I).GE.-10..AND.DU(I).LT.1.0) THEN
          CBC_U(I) = AMAX1(BFRIC,
     *    .16/ALOG((ZZ(KBM1)-Z(KB))*1.0/Z0B)**2)
          CBC_U(I) = CBC_U(I)*2
      ENDIF
      
      IF(DV(I).GE.3.0) THEN
          CBC_V(I) = AMAX1(BFRIC,
     *    .16/ALOG((ZZ(KBM1)-Z(KB))*DV(I)/Z0B)**2)*CD_TEST
      ELSEIF(DV(I).GE.1..AND.DV(I).LT.3.0) THEN
          CBC_V(I) = AMAX1(BFRIC,
     *    .16/ALOG((ZZ(KBM1)-Z(KB))*3.0/Z0B)**2)*CD_TEST*1.3
      ELSEIF(HV(I).GE.-10..AND.DV(I).LT.1.0) THEN
          CBC_V(I) = AMAX1(BFRIC,
     *    .16/ALOG((ZZ(KBM1)-Z(KB))*1.0/Z0B)**2)
          CBC_V(I) = CBC_V(I)*2
      ENDIF
#ifdef MODULE_SED
      IF (D(I).GE.3.0) THEN  
          CBC_UV(I) = AMAX1(BFRIC,.16/
     *    ALOG((ZZ(KBM1)-Z(KB))*D(I)/Z0B)**2)*CD_TEST
      ELSEIF (D(I).GE.1..AND.D(I).LT.3.0) THEN
          CBC_UV(I) = AMAX1(BFRIC,.16/                 
     *    ALOG((ZZ(KBM1)-Z(KB))*3.0/Z0B)**2)*CD_TEST*1.3
      ELSEIF (H(I).GE.-10..AND.D(I).LT.1.0) THEN  
          CBC_UV(I) = AMAX1(BFRIC,.16/                 
     *    ALOG((ZZ(KBM1)-Z(KB))*1.0/Z0B)**2)*2    
      ENDIF
#endif
      ENDDO


c	IF(LOG_CBCADJ)THEN
c      DO I=1,NUMCBCADJ
c      II=CBCADJ(I)
c      IF (FSM(II).GT.0.) THEN
     !! IF(DU(II).LE.2.0.AND.DU(II).GT.0.0) THEN
     !!     CBCMIN = 0.00001
     !!     CBC_U(II) = AMAX1(CBCMIN,(9.8*(0.055)**2.
     !!*    /(DU(II)**(1./3.))))*CBC_COF(II)  
     !! ELSEIF(DU(II).LE.3.0.AND.DU(II).GT.2.0) THEN
     !!     CBCMIN = 0.00001
     !!     CBC_U(II) = AMAX1(CBCMIN,(9.8*(0.040)**2.
     !!*    /(DU(II)**(1./3.))))*CBC_COF(II)  
     !! ELSEIF(DU(II).LE.4.AND.DU(II).GT.3.0) THEN
     !!     CBCMIN = 0.00001
     !!     CBC_U(II) = AMAX1(CBCMIN,(9.8*(0.032)**2.
     !!*    /(DU(II)**(1./3.))))*CBC_COF(II)  
     !! ELSEIF(DU(II).LE.7.0.AND.DU(II).GT.4)THEN    
     !!     CBCMIN = 0.00001
     !!     CBC_U(II) = AMAX1(CBCMIN,(9.8*(0.027)**2.
     !!*    /(DU(II)**(1./3.))))*CBC_COF(II)
     !! ELSEIF(DU(II).LE.10.0.AND.DU(II).GT.7.0)THEN    
     !!     CBCMIN = 0.00001
     !!     CBC_U(II) = AMAX1(CBCMIN,(9.8*(0.014)**2.
     !!*    /(DU(II)**(1./3.))))*CBC_COF(II)
     !! ELSEIF(DU(II).GT.10.0)THEN    
     !!     CBCMIN = 0.00001
     !!     CBC_U(II) = AMAX1(CBCMIN,(9.8*(0.01)**2.
     !!*/(DU(II)**(1./3.))))*CBC_COF(II)
     !! ENDIF
     !! 
     !! IF(DV(II).LE.2.0.AND.DV(II).GT.0.0) THEN
     !!     CBCMIN = 0.00001
     !!     CBC_V(II) = AMAX1(CBCMIN,(9.8*(0.055)**2.
     !!*    /(DV(II)**(1./3.))))*CBC_COF(II)  
     !! ELSEIF(DV(II).LE.3.0.AND.DV(II).GT.2.0) THEN
     !!     CBCMIN = 0.00001
     !!     CBC_V(II) = AMAX1(CBCMIN,(9.8*(0.040)**2.
     !!*    /(DV(II)**(1./3.))))*CBC_COF(II)  
     !! ELSEIF(DV(II).LE.4.AND.DV(II).GT.3.0) THEN
     !!     CBCMIN = 0.00001
     !!     CBC_V(II) = AMAX1(CBCMIN,(9.8*(0.032)**2.
     !!*    /(DV(II)**(1./3.))))*CBC_COF(II)  
     !! ELSEIF(DV(II).LE.7.0.AND.DV(II).GT.4)THEN    
     !!     CBCMIN = 0.00001
     !!     CBC_V(II) = AMAX1(CBCMIN,(9.8*(0.027)**2.
     !!*    /(DV(II)**(1./3.))))*CBC_COF(II)
     !! ELSEIF(DV(II).LE.10.0.AND.DV(II).GT.7.0)THEN    
     !!     CBCMIN = 0.00001
     !!     CBC_V(II) = AMAX1(CBCMIN,(9.8*(0.014)**2.
     !!*    /(DV(II)**(1./3.))))*CBC_COF(II)
     !! ELSEIF(DV(II).GT.10.0)THEN    
     !!     CBCMIN = 0.00001
     !!     CBC_V(II) = AMAX1(CBCMIN,(9.8*(0.01)**2.
     !!*    /(DV(II)**(1./3.))))*CBC_COF(II)
     !! ENDIF
     !! IF(D(II).LE.2.0.AND.D(II).GT.0.0) THEN
     !!     CBCMIN = 0.00001
     !!     CBC_UV(II) = AMAX1(CBCMIN,(9.8*(0.055)**2.
     !!* /(D(II)**(1./3.))))*CBC_COF(II)  
     !! ELSEIF(D(II).LE.3.0.AND.D(II).GT.2.0) THEN
     !!     CBCMIN = 0.00001
     !!     CBC_UV(II) = AMAX1(CBCMIN,(9.8*(0.040)**2.
     !!*    /(D(II)**(1./3.))))*CBC_COF(II)  
     !! ELSEIF(D(II).LE.4.AND.D(II).GT.3.0) THEN
     !!     CBCMIN = 0.00001
     !!     CBC_UV(II) = AMAX1(CBCMIN,(9.8*(0.032)**2.
     !!*    /(D(II)**(1./3.))))*CBC_COF(II)  
     !! ELSEIF(D(II).LE.7.0.AND.D(II).GT.6)THEN    
     !!     CBCMIN = 0.00001
     !!     CBC_UV(II) = AMAX1(CBCMIN,(9.8*(0.027)**2.
     !!*    /(D(II)**(1./3.))))*CBC_COF(II)
     !! ELSEIF(D(II).LE.10.0.AND.D(II).GT.7.0)THEN    
     !!     CBCMIN = 0.00001
     !!     CBC_UV(II) = AMAX1(CBCMIN,(9.8*(0.014)**2.
     !!*    /(D(II)**(1./3.))))*CBC_COF(II)
     !! ELSEIF(D(II).GT.10.0)THEN    
     !!     CBCMIN = 0.00001
     !!     CBC_UV(II) = AMAX1(CBCMIN,(9.8*(0.01)**2.
     !!*    /(D(II)**(1./3.))))*CBC_COF(II)
     !! ENDIF
c      ENDIF
c      ENDDO
c      ENDIF
      
   	RETURN
      END SUBROUTINE BROUGH_USERDEFINED
      