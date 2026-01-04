#include "DEFS.h"

	SUBROUTINE INFO_EXCH3D_GPU(VA,MASK,NAG_EVG,WT_EVG,NAG_IVG,WT_IVG)
      
!======================================================================
! COMMENTS
!
! VTYPE = 0
! 1: ELEVATION
! 2: U VARIABLE
! 3: V VARIABLE
! 4: OTHERS
!======================================================================
  
      USE MOD_GLOBAL
  
      IMPLICIT NONE

      
      INTEGER I,J,K
!	INTEGER VTYPE
      INTEGER NAG_EVG(NMAX,N_CTRD_EVG),NAG_IVG(NMAX,N_CTRD_IVG)
      REAL WEIALL, SUMALL
      REAL VAMAX
      REAL VA(N_CTRDP1,KB)
!      REAL VA_AG(N_CTRD_AG) !,VA_EVG(N_CTRD_EVG),VA_IVG(N_CTRD_IVG)
      REAL MASK(N_CTRDP1),L4IE !, L4IE(N_CTRD)
c  CBR: Moved to MOD_Global:
c      REAL, ALLOCATABLE :: MASK4VAR(:,:)
      REAL, PARAMETER :: EPSILON = 1E-10
      REAL WT_EVG(NMAX,N_CTRD_EVG),WT_IVG(NMAX,N_CTRD_IVG)
      REAL WTM1,WTM2
      REAL VA_INT1,VA_INT2,VA_INT3,VA_INT4
      REAL MASK_INT1,MASK_INT2,MASK_INT3,MASK_INT4
      
      INTEGER N_CTRD_AEG
	
#ifdef INT_UAEI
      LOGICAL FLAG1, FLAG2, FLAG3
#endif

#ifdef INT_HPI
      LOGICAL FLAG1, FLAG2, FLAG3
      REAL MASK_INT(4), VA_INT(4) !, COF(3)
      REAL COF1,COF2,COF3 
      INTEGER, PARAMETER ::  IND(4)=(/3,5,1,7/)
c      INTEGER IND(4)
#endif

# ifdef INT_HAEI
      INTEGER, PARAMETER ::  IND(5)=(/3,5,1,7,9/)
c      INTEGER IND(5)
      REAL MASK_INT(5), VA_INT(5)
      REAL PHIL, PHIR, EPS, RL, RTR, A1, B1, BETAL, BETAR, KAPPA
      REAL, PARAMETER :: EPSON = 0.0001
      LOGICAL FLAG1, FLAG2, FLAG3
# endif


      N_CTRD_AEG=N_CTRD_AG+N_CTRD_EVG
      
c      VA_AG = VA(1:N_CTRD_AG)
c      VA_EVG = VA(N_CTRD_AG+1:N_CTRD_AG+N_CTRD_EVG)
c	VA_IVG = VA(N_CTRD_AG+N_CTRD_EVG+1:N_CTRD)
	

!======================================================================
!                         INTERPOLATION
!======================================================================
      
!**********************************************************************
! 0TH ORDER UNIFORM
#ifdef INT_UNI
	!ELSEIF (INTSCH.EQ.'UNI') THEN
	  
		
		DO I = 1,N_CTRD_EVG
		  IF (MASK(N_CTRD_AG+I).GT.0.0) THEN
			VA(N_CTRD_AG+I) = VA(NAG_EVG(1,I))
		  ENDIF
		ENDDO
		
        
	
!**********************************************************************
! INVERSE DISTANCE WEIGHTING INTERPOLATION
! / INVERSE BILINEAR INTERPOLATION
#elif defined INT_IDW || defined INT_LINEAR
	!IF (INTSCH.EQ.'IDW') THEN
	  
	
        DO I = 1,N_CTRD_EVG
          IF (MASK(N_CTRD_AG+I) .GT. 0.0) THEN
            WEIALL = WT_EVG(1,I)*MASK(NAG_EVG(1,I))
     *               +WT_EVG(2,I)*MASK(NAG_EVG(2,I))
     *               +WT_EVG(3,I)*MASK(NAG_EVG(3,I))
     *               +WT_EVG(4,I)*MASK(NAG_EVG(4,I))
!            IF (WEIALL .EQ. 0.) CYCLE
            IF (WEIALL /= 0.)
       VA(N_CTRD_AG+I) = (WT_EVG(1,I)*MASK(NAG_EVG(1,I))*VA(NAG_EVG(1,I))
     *              +WT_EVG(2,I)*MASK(NAG_EVG(2,I))*VA(NAG_EVG(2,I))
     *              +WT_EVG(3,I)*MASK(NAG_EVG(3,I))*VA(NAG_EVG(3,I))
     *              +WT_EVG(4,I)*MASK(NAG_EVG(4,I))*VA(NAG_EVG(4,I)))
     *              /WEIALL

          ENDIF
	  ENDDO
	  
	  

!**********************************************************************
! QUADRATIC INTERPOLATION
#elif defined INT_QDT
        DO I = 1,N_CTRD_EVG
          IF ((MASK(N_CTRD_AG+I) .GT. 0.0)) THEN
            WEIALL = 
     * WT_EVG(1,I)*WT_EVG(7,I)*MASK(NAG_EVG(1,I))
     * +WT_EVG(2,I)*WT_EVG(7,I)*MASK(NAG_EVG(2,I))
     * +WT_EVG(3,I)*WT_EVG(8,I)*MASK(NAG_EVG(3,I))
     * +WT_EVG(4,I)*WT_EVG(8,I)*MASK(NAG_EVG(4,I))
     * +WT_EVG(5,I)*WT_EVG(9,I)*MASK(NAG_EVG(5,I))
     * +WT_EVG(6,I)*WT_EVG(9,I)*MASK(NAG_EVG(6,I))
!            IF (WEIALL.LT.0.2) CYCLE
            IF (WEIALL>=0.2) THEN
                VA(N_CTRD_AG+I) = 
     * (WT_EVG(1,I)*WT_EVG(7,I)*MASK(NAG_EVG(1,I))*VA(NAG_EVG(1,I))
     * +WT_EVG(2,I)*WT_EVG(7,I)*MASK(NAG_EVG(2,I))*VA(NAG_EVG(2,I))
     * +WT_EVG(3,I)*WT_EVG(8,I)*MASK(NAG_EVG(3,I))*VA(NAG_EVG(3,I))
     * +WT_EVG(4,I)*WT_EVG(8,I)*MASK(NAG_EVG(4,I))*VA(NAG_EVG(4,I))
     * +WT_EVG(5,I)*WT_EVG(9,I)*MASK(NAG_EVG(5,I))*VA(NAG_EVG(5,I))
     * +WT_EVG(6,I)*WT_EVG(9,I)*MASK(NAG_EVG(6,I))*VA(NAG_EVG(6,I)))
     * /WEIALL

              !LIMIT VA_EVG
                VAMAX = MAX(ABS(VA(NAG_EVG(1,I))),
     * ABS(VA(NAG_EVG(2,I))), ABS(VA(NAG_EVG(3,I))),
     * ABS(VA(NAG_EVG(4,I))), ABS(VA(NAG_EVG(5,I))),
     * ABS(VA(NAG_EVG(6,I))))
                IF (ABS(VA(N_CTRD_AG+I)).GT.VAMAX) THEN
                  VA(N_CTRD_AG+I) = VA(NAG_EVG(1,I))
                ENDIF
             ENDIF
          ENDIF
	  ENDDO
	  

        
!**********************************************************************
! HSIMT PARABOLIC INTERPOLATION
#elif defined INT_HPI

      !---------------------------------------------------
      !LOGICAL SWITCHES FOR INFO EXCHANGE & VARIABLE MASKS
c      L4IE = 0.
c  CBR: Moved to ALLOC_VARS:
c      ALLOCATE (MASK4VAR(8,N_CTRD_EVG))
c      MASK4VAR = 0. !CBR
      
      !OTHER VARIABLES
c        L4IE = MASK(1:N_CTRD)
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c        DO I = 1,N_CTRD_EVG
c          DO J = 1,8
        doconcurrent (I = 1:N_CTRD_EVG,J=1:8)
            MASK4VAR(J,I) = MASK(NAG_EVG(J,I))
          ENDDO
c        ENDDO
c!$OMP END TARGET
      
      !------------------------------------
      !INTERPOLATION
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) !private(MASK_INT, VA_INT)
c      DO K=1,KB
c      DO I = 1,N_CTRD_EVG
      doconcurrent (K=1:KB,I=1:N_CTRD_EVG) 
     * local(FLAG1,FLAG2,FLAG3,WTM1,WTM2,VA_INT1,MASK_INT1,
     * VA_INT2,MASK_INT2,VA_INT3,MASK_INT3,VA_INT4,MASK_INT4,
     * COF1,COF2,COF3,WEIALL,VAMAX)
!          L4IE=MASK(N_CTRD_AG+I)
          
c          IF (L4IE(N_CTRD_AG+I).GT.0.0) THEN
!          IF (L4IE.GT.0.0) THEN
          IF (MASK(N_CTRD_AG+I).GT.0.0) THEN
          !FLAGS
          !1 FOR WESTERN/EASTERN BOUNDARY 
          !O FOR SOUTHERN/NORTHERN BOUNDARY
          FLAG1 = (NIP2(NAG_EVG(1,I)).GT.N_CTRD_AG) .OR. 
     * (NIM2(NAG_EVG(1,I)).GT.N_CTRD_AG)
          !U DIRECTION
          FLAG2 = U(NAG_EVG(1,I),K).GT.0.
          !V DIRECTION
          FLAG3 = V(NAG_EVG(1,I),K).GT.0.
            
          !DETERMINE MASK & VARIABLES FOR INTERPOLATION
c          MASK_INT = 0.
c          IND = (/3,5,1,7/)
          !I-2 TO I+1
c --- J=1:       
              WTM1=WT_EVG(3,I)*MASK4VAR(3,I)
              WTM2=WT_EVG(4,I)*MASK4VAR(4,I)
            VA_INT1 = 
     * (VA(NAG_EVG(3,I),K)*WTM1
     * +VA(NAG_EVG(4,I),K)*WTM2)
     * /(WTM1+WTM2+EPSILON)
            IF 
     * ((MASK4VAR(3,I)+MASK4VAR(4,I)).NE.0.) THEN
              MASK_INT1 = 1.
            ELSE
              MASK_INT1 = 0.
            ENDIF
c --- J=2:       
              WTM1=WT_EVG(5,I)*MASK4VAR(5,I)
              WTM2=WT_EVG(6,I)*MASK4VAR(6,I)
            VA_INT2 = 
     * (VA(NAG_EVG(5,I),K)*WTM1
     * +VA(NAG_EVG(6,I),K)*WTM2)
     * /(WTM1+WTM2+EPSILON)
            IF 
     * ((MASK4VAR(5,I)+MASK4VAR(6,I)).NE.0.) THEN
              MASK_INT2 = 1.
            ELSE
              MASK_INT2 = 0.
            ENDIF
c --- J=3:       
              WTM1=WT_EVG(1,I)*MASK4VAR(1,I)
              WTM2=WT_EVG(2,I)*MASK4VAR(2,I)
            VA_INT3 = 
     * (VA(NAG_EVG(1,I),K)*WTM1
     * +VA(NAG_EVG(2,I),K)*WTM2)
     * /(WTM1+WTM2+EPSILON)
            IF 
     * ((MASK4VAR(1,I)+MASK4VAR(2,I)).NE.0.) THEN
              MASK_INT3 = 1.
            ELSE
              MASK_INT3 = 0.
            ENDIF
c --- J=4:       
              WTM1=WT_EVG(7,I)*MASK4VAR(7,I)
              WTM2=WT_EVG(8,I)*MASK4VAR(8,I)
            VA_INT4 = 
     * (VA(NAG_EVG(7,I),K)*WTM1
     * +VA(NAG_EVG(8,I),K)*WTM2)
     * /(WTM1+WTM2+EPSILON)
            IF 
     * ((MASK4VAR(7,I)+MASK4VAR(8,I)).NE.0.) THEN
              MASK_INT4 = 1.
            ELSE
              MASK_INT4 = 0.
            ENDIF
            
            
          IF ((FLAG1.AND.FLAG2).OR.((.NOT.FLAG1).AND.FLAG3)) THEN
          !U/V > 0.
              
            IF ((MASK_INT1*MASK_INT2*MASK_INT3).NE.0.) THEN
            !HSIMT PARABOLIC INTERPOLATION
                
              !COEFFICIENTS
              COF1 = 0.5*(VA_INT3+VA_INT1)-VA_INT2
              COF2 = VA_INT3-VA_INT2
              COF3 = (1./3.)*VA_INT3+(5./6.)*VA_INT2-(1./6.)*VA_INT1
              !COEFFICIENTS
c              COF1 = 
c     * 0.5*(VA_INT3-VA_INT2)-0.5*(VA_INT2-VA_INT1)
c              COF2 = VA_INT3-VA_INT2
c              COF3 = 1./3.*(VA_INT3-VA_INT2)
c     * +1./6.*(VA_INT2-VA_INT1)+VA_INT2
                
              !INTERPOLATE
              VA(N_CTRD_AG+I,K) = COF1*XIC(I)**2+COF2*XIC(I)+COF3
                
            ELSE
            !UPWIND AEI
              !OTHER VARIABLES
                WEIALL = 
     * WT_EVG(1,I)*WT_EVG(11,I)+WT_EVG(2,I)*WT_EVG(11,I)
     * +WT_EVG(5,I)*WT_EVG(12,I)+WT_EVG(6,I)*WT_EVG(12,I)
                IF (WEIALL .EQ. 0.) THEN
!                  CYCLE
                ELSE
                  VA(N_CTRD_AG+I,K) = 
     * (WT_EVG(1,I)*WT_EVG(11,I)*VA(NAG_EVG(1,I),K)
     * +WT_EVG(2,I)*WT_EVG(11,I)*VA(NAG_EVG(2,I),K)
     * +WT_EVG(5,I)*WT_EVG(12,I)*VA(NAG_EVG(5,I),K)
     * +WT_EVG(6,I)*WT_EVG(12,I)*VA(NAG_EVG(6,I),K))
     * /WEIALL
                ENDIF
   
            ENDIF
                 
          ELSE 
          !U/V<0.
              
            IF ((MASK_INT2*MASK_INT3*MASK_INT4).NE.0.) THEN
            !HSIMT PARABOLIC INTERPOLATION
                
              !COEFFICIENTS
              COF1 = 0.5*(VA_INT4+VA_INT2)-VA_INT3
              COF2 = VA_INT3-VA_INT2
              COF3 = (1./3.)*VA_INT2+(5./6.)*VA_INT3-(1./6.)*VA_INT4
              !COEFFICIENTS
c              COF1 = 
c     * 0.5*(VA_INT4-VA_INT3)-0.5*(VA_INT3-VA_INT2)
c              COF2 = VA_INT3-VA_INT2
c              COF3 = -1./6.*(VA_INT4-VA_INT3)
c     * -1./3.*(VA_INT3-VA_INT2)+VA_INT3
                
              !INTERPOLATE
                VA(N_CTRD_AG+I,K) = COF1*XIC(I)**2+COF2*XIC(I)+COF3
                
            ELSE
            !UPWIND AEI
              
              !OTHER VARIABLES
                WEIALL = 
     * WT_EVG(1,I)*WT_EVG(9,I)*MASK4VAR(1,I)
     * +WT_EVG(2,I)*WT_EVG(9,I)*MASK4VAR(2,I)
     * +WT_EVG(7,I)*WT_EVG(10,I)*MASK4VAR(7,I)
     * +WT_EVG(8,I)*WT_EVG(10,I)*MASK4VAR(8,I)
                IF (WEIALL .EQ. 0.) THEN
!                  CYCLE
                ELSE
                  VA(N_CTRD_AG+I,K) = 
     * (WT_EVG(1,I)*WT_EVG(9,I)*MASK4VAR(1,I)*VA(NAG_EVG(1,I),K)
     * +WT_EVG(2,I)*WT_EVG(9,I)*MASK4VAR(2,I)*VA(NAG_EVG(2,I),K)
     * +WT_EVG(7,I)*WT_EVG(10,I)*MASK4VAR(7,I)*VA(NAG_EVG(7,I),K)
     * +WT_EVG(8,I)*WT_EVG(10,I)*MASK4VAR(8,I)*VA(NAG_EVG(8,I),K))
     * /WEIALL
                ENDIF
                
   
            ENDIF
                
          ENDIF

          !LIMIT VA_EVG
          VAMAX = MAX(ABS(VA(NAG_EVG(1,I),K)),
     * ABS(VA(NAG_EVG(2,I),K)), ABS(VA(NAG_EVG(3,I),K)),
     * ABS(VA(NAG_EVG(4,I),K)), ABS(VA(NAG_EVG(5,I),K)),
     * ABS(VA(NAG_EVG(6,I),K)), ABS(VA(NAG_EVG(7,I),K)),
     * ABS(VA(NAG_EVG(8,I),K)))
          IF (ABS(VA(N_CTRD_AG+I,K)).GT.VAMAX) THEN
            VA(N_CTRD_AG+I,K) = VA(NAG_EVG(1,I),K)
          ENDIF
        ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET


!**********************************************************************
! UPWIND AEI (UPWIND ADVECTION-EQUIVALENT INTERPOLATION)
#elif defined INT_UAEI

	
        DO I = 1,N_CTRD_EVG
          IF ((MASK(N_CTRD_AG+I) .GT. 0.0)) THEN
            FLAG1 = (NIP2(NAG_EVG(1,I)).GT.N_CTRD_AG) .OR. 
     * (NIM2(NAG_EVG(1,I)).GT.N_CTRD_AG)
            FLAG2 = U(NAG_EVG(1,I),KIN).GT.0.
            FLAG3 = V(NAG_EVG(1,I),KIN).GT.0.
			
            IF ((FLAG1.AND.FLAG2).OR.((.NOT.FLAG1).AND.FLAG3)) THEN
			  
              WEIALL = 
     * WT_EVG(1,I)*WT_EVG(9,I)+WT_EVG(2,I)*WT_EVG(9,I)
     * +WT_EVG(5,I)*WT_EVG(10,I)+WT_EVG(6,I)*WT_EVG(10,I)
              IF (WEIALL .EQ. 0.) THEN
!                CYCLE
              ELSE
                VA(N_CTRD_AG+I) = 
     * (WT_EVG(1,I)*WT_EVG(9,I)*VA(NAG_EVG(1,I))
     * +WT_EVG(2,I)*WT_EVG(9,I)*VA(NAG_EVG(2,I))
     * +WT_EVG(5,I)*WT_EVG(10,I)*VA(NAG_EVG(5,I))
     * +WT_EVG(6,I)*WT_EVG(10,I)*VA(NAG_EVG(6,I)))
     * /WEIALL
              ENDIF
			  
            ELSE
			  
              WEIALL = 
     * WT_EVG(1,I)*WT_EVG(7,I)+WT_EVG(2,I)*WT_EVG(7,I)
     * +WT_EVG(3,I)*WT_EVG(8,I)+WT_EVG(4,I)*WT_EVG(8,I)
              IF (WEIALL .EQ. 0.) THEN
!                CYCLE
              ELSE
                VA(N_CTRD_AG+I) = 
     * (WT_EVG(1,I)*WT_EVG(7,I)*VA(NAG_EVG(1,I))
     * +WT_EVG(2,I)*WT_EVG(7,I)*VA(NAG_EVG(2,I))
     * +WT_EVG(3,I)*WT_EVG(8,I)*VA(NAG_EVG(3,I))
     * +WT_EVG(4,I)*WT_EVG(8,I)*VA(NAG_EVG(4,I)))
     * /WEIALL
                !LIMIT VA_EVG
                VAMAX = MAX(ABS(VA(NAG_EVG(1,I)))
     * ,ABS(VA(NAG_EVG(2,I))),ABS(VA(NAG_EVG(3,I)))
     * ,ABS(VA(NAG_EVG(4,I))))
                IF (ABS(VA(N_CTRD_AG+I)).GT.VAMAX) THEN
                  VA(N_CTRD_AG+I) = VA(NAG_EVG(1,I))
                ENDIF
                
              ENDIF
			
            ENDIF

          ENDIF
	  ENDDO
	  

        
!**********************************************************************
! HSIMT AEI (HSIMT ADVECTION-EQUIVALENT INTERPOLATION)
#elif defined INT_HAEI

      !-----------------------------------------------------!
      ! LOGICAL SWITCHES FOR INFO EXCHANGE & VARIABLE MASKS !
      !-----------------------------------------------------!
c      L4IE = 0.
c  CBR: Moved to ALLOC_VARS:
c      ALLOCATE (MASK4VAR(10,N_CTRD_EVG))
c      MASK4VAR = 0. !CBR
      
      !OTHER VARIABLES
c        L4IE = MASK(1:N_CTRD)
        DO I = 1,N_CTRD_EVG
          DO J = 1,10
            MASK4VAR(J,I) = MASK(NAG_EVG(J,I))
          ENDDO
        ENDDO

      !--------------------!
      ! HAEI INTERPOLATION !
      !--------------------!
      DO I = 1,N_CTRD_EVG
!              L4IE=MASK(N_CTRD_AG+I)
          
c          IF (L4IE(N_CTRD_AG+I).GT.0.0) THEN
!          IF (L4IE.GT.0.0) THEN
          IF (MASK(N_CTRD_AG+I).GT.0.0) THEN
          !DEFINE VALUES & MASKS OF REFERENCE POINTS
          MASK_INT = 0.
c          IND = (/3,5,1,7,9/)
          
          !J-2 TO J+2
          DO J = 1,5
            VA_INT(J) = 
     * (VA(NAG_EVG(IND(J),I))*WT_EVG(IND(J),I)
     * *MASK4VAR(IND(J),I)
     * +VA(NAG_EVG(IND(J)+1,I))*WT_EVG(IND(J)+1,I)
     * *MASK4VAR(IND(J)+1,I))
     * /(WT_EVG(IND(J),I)*MASK4VAR(IND(J),I)
     * +WT_EVG(IND(J)+1,I)*MASK4VAR(IND(J)+1,I)+EPSILON)
            IF 
     * ((MASK4VAR(IND(J),I)+MASK4VAR(IND(J)+1,I)).NE.0.) THEN
              MASK_INT(J) = 1.
            ENDIF
          ENDDO
          
          !EPS & KAPPA
          ! OTHER VARIABLES
            EPS = EPSC(I)
          KAPPA = 1-ABS(EPS)
          
          !--------------------!
          ! CALCULATE LEFT PHI !
          !--------------------!
          ! RL & RTR
          IF (ABS(VA_INT(3)-VA_INT(2)).LE.EPSON) THEN
            RL = 0.
            RTR = 0.
          ELSE
            RL = (VA_INT(2)-VA_INT(1))/(VA_INT(3)-VA_INT(2))
            RTR = (VA_INT(4)-VA_INT(3))/(VA_INT(3)-VA_INT(2))
          ENDIF
          
          !FLAGS
          !1 FOR WESTERN/EASTERN BOUNDARY 
          !O FOR SOUTHERN/NORTHERN BOUNDARY
          FLAG1 = (NIP2(NAG_EVG(1,I)).GT.N_CTRD_AG) .OR. 
     * (NIM2(NAG_EVG(1,I)).GT.N_CTRD_AG)
          !U DIRECTION
          FLAG2 = U(NAG_EVG(1,I),KIN).GT.0.
          !V DIRECTION
          FLAG3 = V(NAG_EVG(1,I),KIN).GT.0.
          
          IF ((FLAG1.AND.FLAG2).OR.((.NOT.FLAG1).AND.FLAG3)) THEN
          !U/V > 0.
            IF (MASK_INT(1).EQ.0. .OR. MASK_INT(3).EQ.0.) THEN
              PHIL = VA_INT(2)
            ELSE
              A1 = 1./4.*KAPPA+1./2.-1./(12.*KAPPA)
              B1 = -1./4.*KAPPA+1./2.+1./(12*KAPPA)
              BETAL = A1+B1*RL
              PHIL = VA_INT(2)+0.5*BETAL*KAPPA*(VA_INT(3)-VA_INT(2))
            ENDIF
            
          ELSE
          !U/V < 0.
            IF (MASK_INT(4).EQ.0. .OR. MASK_INT(2).EQ.0.) THEN
              PHIL = VA_INT(3)
            ELSE
              A1 = 1./4.*KAPPA+1./2.-1./(12.*KAPPA)
              B1 = -1./4.*KAPPA+1./2.+1./(12*KAPPA)
              BETAR = A1+B1*RTR
              PHIL = VA_INT(3)-0.5*BETAR*KAPPA*(VA_INT(3)-VA_INT(2))
            ENDIF
            
          ENDIF
          
          !---------------------!
          ! CALCULATE RIGHT PHI !
          !---------------------!
          ! RL & RTR
          IF (ABS(VA_INT(4)-VA_INT(3)).LE.EPSON) THEN
            RL = 0.
            RTR = 0.
          ELSE
            RL = (VA_INT(3)-VA_INT(2))/(VA_INT(4)-VA_INT(3))
            RTR = (VA_INT(5)-VA_INT(4))/(VA_INT(4)-VA_INT(3))
          ENDIF
          
          !FLAGS
          !1 FOR WESTERN/EASTERN BOUNDARY 
          !O FOR SOUTHERN/NORTHERN BOUNDARY
          FLAG1 = (NIP2(NAG_EVG(1,I)).GT.N_CTRD_AG) .OR. 
     * (NIM2(NAG_EVG(1,I)).GT.N_CTRD_AG)
          !U DIRECTION
          FLAG2 = U(NAG_EVG(7,I),KIN).GT.0.
          !V DIRECTION
          FLAG3 = V(NAG_EVG(7,I),KIN).GT.0.
          
          IF ((FLAG1.AND.FLAG2).OR.((.NOT.FLAG1).AND.FLAG3)) THEN
          !U/V > 0.
            IF (MASK_INT(2).EQ.0. .OR. MASK_INT(4).EQ.0.) THEN
              PHIR = VA_INT(3)
            ELSE
              A1 = 1./4.*KAPPA+1./2.-1./(12.*KAPPA)
              B1 = -1./4.*KAPPA+1./2.+1./(12*KAPPA)
              BETAL = A1+B1*RL
              PHIR = VA_INT(3)+0.5*BETAL*KAPPA*(VA_INT(4)-VA_INT(3))
            ENDIF
            
          ELSE
          !U/V < 0.
            IF (MASK_INT(5).EQ.0. .OR. MASK_INT(3).EQ.0.) THEN
              PHIR = VA_INT(4)
            ELSE
              A1 = 1./4.*KAPPA+1./2.-1./(12.*KAPPA)
              B1 = -1./4.*KAPPA+1./2.+1./(12*KAPPA)
              BETAR = A1+B1*RTR
              PHIR = VA_INT(4)-0.5*BETAR*KAPPA*(VA_INT(4)-VA_INT(3))
            ENDIF
            
          ENDIF
          
          VA(N_CTRD_AG+I) = VA_INT(3)+EPS*(PHIR-PHIL)
        ENDIF
      ENDDO



        
!**********************************************************************	  
	!ENDIF
#else
      PRINT*, 'interpolation scheme undefined!'
	PAUSE
	STOP
#endif

!======================================================================
!                             END: INTERPOLATION
!======================================================================
  
  

!======================================================================
!                                 UPDATE
!======================================================================
	
!**********************************************************************
! DIRECTLY-REPLACING
#ifdef UD_DR
	!ELSEIF (UDSCH.EQ.'DR') THEN
		
		DO I = 1,N_CTRD_IVG
		  IF (MASK(N_CTRD_AEG+I).GT.0.0) THEN
			VA(N_CTRD_AEG+I) = VA(NAG_IVG(1,I))
		  ENDIF
		ENDDO
		

!**********************************************************************
! INVERSE DISTANCE WEIGHTING INTERPOLATION
#elif defined UD_IDW
      !IF (UDSCH.EQ.'IDW') THEN
	  
	
        DO I = 1,N_CTRD_IVG
          IF (MASK(N_CTRD_AEG+I) .GT. 0.0) THEN
            WEIALL = WT_IVG(1,I)*MASK(NAG_IVG(1,I))
     *               +WT_IVG(2,I)*MASK(NAG_IVG(2,I))
     *               +WT_IVG(3,I)*MASK(NAG_IVG(3,I))
     *               +WT_IVG(4,I)*MASK(NAG_IVG(4,I))
!            IF (WEIALL .EQ. 0.) CYCLE
            IF (WEIALL /= 0.) THEN
       VA(N_CTRD_AEG+I) = (WT_IVG(1,I)*MASK(NAG_IVG(1,I))*VA(NAG_IVG(1,I))
     *              +WT_IVG(2,I)*MASK(NAG_IVG(2,I))*VA(NAG_IVG(2,I))
     *              +WT_IVG(3,I)*MASK(NAG_IVG(3,I))*VA(NAG_IVG(3,I))
     *              +WT_IVG(4,I)*MASK(NAG_IVG(4,I))*VA(NAG_IVG(4,I)))
     *              /WEIALL
            ENDIF

          ENDIF
	  ENDDO
	  
	  
!**********************************************************************
!AREA-AVERAGING
#elif defined UD_AVE
	  DO I = 1,N_CTRD_IVG
		IF (MASK(N_CTRD_AEG+I).GT.0.0) THEN
		
		  WEIALL = 0
		  SUMALL = 0
		  DO J = 1,NMAX
		    !IF (NAG_IVG(J,I).NE.0) THEN
			  !IF (MASK(NAG_IVG(J,I)).GT.0.0) THEN
		  	  SUMALL = SUMALL+VA(NAG_IVG(J,I))*WT_IVG(J,I)
			    WEIALL = WEIALL+WT_IVG(J,I)
			  !ENDIF
		    !ELSE
		  	 ! EXIT
		    !ENDIF
		  ENDDO 
		
!		  IF (WEIALL.EQ.0) CYCLE
		  IF (WEIALL/=0) VA(N_CTRD_AEG+I) = SUMALL/WEIALL
		  
		ENDIF
	  ENDDO
	  

	
!**********************************************************************
!9-POINT SHAPIRO FILTERING
#elif defined UD_SF
	  DO I = 1,N_CTRD_IVG
		IF (MASK(N_CTRD_AEG+I).GT.0.0) THEN
		
		  WEIALL = 0
		  SUMALL = 0
		  DO J = 1,NMAX
		    IF (NAG_IVG(J,I).NE.1) THEN
			  !IF (MASK(NAG_IVG(J,I)).GT.0.0) THEN
              SUMALL = SUMALL+VA(NAG_IVG(J,I))*WT_IVG(J,I)
			    WEIALL = WEIALL+WT_IVG(J,I)
			  !ENDIF
		    ELSE
		  	  EXIT
		    ENDIF
		  ENDDO
		
!		  IF (WEIALL.EQ.0) CYCLE
		  IF (WEIALL/=0) VA(N_CTRD_AEG+I) = SUMALL
		  
		ENDIF
	  ENDDO
	  

      
!**********************************************************************
! FULL-WEIGHTING OPERATOR (25-POINT)
#elif defined UD_FWO
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c        DO K=1,KB
c	  DO I = 1,N_CTRD_IVG
	  DO concurrent (I=1:N_CTRD_IVG,K=1:KB) local(WEIALL,SUMALL,J)
		IF (MASK(N_CTRD_AEG+I).GT.0.0) THEN
		
		  WEIALL = 0
		  SUMALL = 0
		  DO J = 1,NMAX
		    IF (NAG_IVG(J,I).NE.1) THEN
			  !IF (MASK(NAG_IVG(J,I)).GT.0.0) THEN
              SUMALL = SUMALL+VA(NAG_IVG(J,I),K)*WT_IVG(J,I)
			    WEIALL = WEIALL+WT_IVG(J,I)
			  !ENDIF
		    ELSE
		  	  EXIT
		    ENDIF
		  ENDDO
		
!		  IF (WEIALL.EQ.0) CYCLE
		  IF (WEIALL/=0) VA(N_CTRD_AEG+I,K) = SUMALL
		  
		ENDIF
	  ENDDO
c	  ENDDO
c!$OMP END TARGET
	  

!**********************************************************************
	!ENDIF
      
#else
      PRINT*, 'Update scheme undefined!'
	PAUSE
	STOP
#endif

!======================================================================
!							 END: UPDATE
!======================================================================
	
	
c	VA(N_CTRD_AG+1:N_CTRD_AG+N_CTRD_EVG) = VA_EVG
c	VA(N_CTRD_AG+N_CTRD_EVG+1:N_CTRD) = VA_IVG
  
      RETURN
  
      END SUBROUTINE INFO_EXCH3D_GPU
  