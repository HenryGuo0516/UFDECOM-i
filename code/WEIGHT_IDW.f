
      SUBROUTINE WEIGHT_IDW(XINT,YINT,X,Y,WEIREF,NAGREF)

      ! INVERSE DISTANCE WEIGHT METHOD
      USE MOD_GLOBAL

      IMPLICIT NONE

      INTEGER I
      INTEGER NAG
      INTEGER NAGREF(NMAX)
      REAL, PARAMETER :: PI = 3.1415926/180.
      REAL ENMIN,WNMIN,ESMIN,WSMIN
      REAL DIS,DX,DY
      REAL XINT, YINT
      REAL WEIREF(NMAX)
      REAL X(N_CTRD),Y(N_CTRD)
      LOGICAL ATEAST,ATNORTH

      !================================================================
      ! DETERMINE WEIGHT
      WEIREF = 0.0
      NAGREF = 1

      ENMIN = SEARCHRADIUS
      WNMIN = SEARCHRADIUS
      ESMIN = SEARCHRADIUS
      WSMIN = SEARCHRADIUS
      
      DO NAG = 1,N_CTRD_AG
      IF (COORD .EQ. 'XY') THEN
          DIS = SQRT((X(NAG)-XINT)**2+(Y(NAG)-YINT)**2)
      ELSEIF (COORD .EQ. 'BL') THEN
          DX = (XINT-X(NAG))*PI*ER*COS(0.5*(YINT+Y(NAG))*PI)
          DY = (YINT-Y(NAG))*ER*PI
          DIS = SQRT(DX**2+DY**2)
      ELSE 
          PRINT*, 'WRONG COORDINATE SYSTEM!'
      ENDIF
      IF (DIS<SEARCHRADIUS) THEN
      IF (DIS.LE.0.1) THEN
          DIS = 0.1
          WEIREF = 1.
          NAGREF = NAG
          EXIT  
      ELSE
          ATEAST = X(NAG)>XINT
          ATNORTH = Y(NAG)>YINT
          IF (ATEAST .AND. ATNORTH) THEN
              IF (DIS<ENMIN) THEN
                  ENMIN = DIS
                  WEIREF(1) = 1000.0/DIS**NPOW
                  NAGREF(1) = NAG
              ENDIF
          ELSEIF ((.NOT.ATEAST).AND.ATNORTH) THEN
              IF (DIS<WNMIN) THEN
                  WNMIN = DIS
                  WEIREF(2) = 1000.0/DIS**NPOW
                  NAGREF(2) = NAG
              ENDIF
          ELSEIF (ATEAST.AND.(.NOT.ATNORTH)) THEN
              IF (DIS<ESMIN) THEN
                  ESMIN = DIS
                  WEIREF(3) = 1000.0/DIS**NPOW
                  NAGREF(3) = NAG
              ENDIF
          ELSE 
              IF (DIS<WSMIN) THEN
                  WSMIN = DIS
                  WEIREF(4) = 1000.0/DIS**NPOW
                  NAGREF(4) = NAG
              ENDIF
          ENDIF
      ENDIF
      ENDIF
      ENDDO
      ! END: DETERMINE WEIGHT
      !================================================================

      RETURN
      END SUBROUTINE WEIGHT_IDW
