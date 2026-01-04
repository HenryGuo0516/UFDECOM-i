      SUBROUTINE INTERP_POINT(VA,VAM,MASK,NAGREF,WEIREF)
  
      USE MOD_GLOBAL
  
      IMPLICIT NONE
  
c      INTEGER I
      INTEGER NAGREF(NMAX)
      REAL WEIALL
      REAL VA
      REAL VAM(N_CTRDP1) !,VA_AG(N_CTRD_AG)
      REAL MASK(N_CTRDP1)
      REAL WEIREF(NMAX)

c      VA = 0.0
c      VA_AG = VAM(1:N_CTRD_AG)
  
!==============================================================================
! INTERPOLATION
      WEIALL = WEIREF(1)*MASK(NAGREF(1))+WEIREF(2)*MASK(NAGREF(2))
     * +WEIREF(3)*MASK(NAGREF(3))+WEIREF(4)*MASK(NAGREF(4))
      IF (WEIALL /= 0.) THEN
      VA = (WEIREF(1)*MASK(NAGREF(1))*VAM(NAGREF(1))
     * +WEIREF(2)*MASK(NAGREF(2))*VAM(NAGREF(2))
     * +WEIREF(3)*MASK(NAGREF(3))*VAM(NAGREF(3))
     * +WEIREF(4)*MASK(NAGREF(4))*VAM(NAGREF(4)))/WEIALL
      ELSE
      VA = 0.0
      ENDIF
      
! END: INTERPOLATION
!==============================================================================

      RETURN

      END SUBROUTINE INTERP_POINT

      
      
      
      SUBROUTINE INTERP_POINT3D(VA,VAM,K,MASK,NAGREF,WEIREF)
  
      USE MOD_GLOBAL
  
      IMPLICIT NONE
  
      INTEGER K
      INTEGER NAGREF(NMAX)
      REAL WEIALL
      REAL VA
      REAL VAM(N_CTRDP1,KB) !,VA_AG(N_CTRD_AG)
      REAL MASK(N_CTRDP1)
      REAL WEIREF(NMAX)

c      VA = 0.0
c      VA_AG = VAM(1:N_CTRD_AG)
  
!==============================================================================
! INTERPOLATION
      WEIALL = WEIREF(1)*MASK(NAGREF(1))+WEIREF(2)*MASK(NAGREF(2))
     * +WEIREF(3)*MASK(NAGREF(3))+WEIREF(4)*MASK(NAGREF(4))
      IF (WEIALL /= 0.) THEN
      VA = (WEIREF(1)*MASK(NAGREF(1))*VAM(NAGREF(1),K)
     * +WEIREF(2)*MASK(NAGREF(2))*VAM(NAGREF(2),K)
     * +WEIREF(3)*MASK(NAGREF(3))*VAM(NAGREF(3),K)
     * +WEIREF(4)*MASK(NAGREF(4))*VAM(NAGREF(4),K))/WEIALL
      ELSE
      VA = 0.0
      ENDIF
      
! END: INTERPOLATION
!==============================================================================

      RETURN

      END SUBROUTINE INTERP_POINT3D
      
      