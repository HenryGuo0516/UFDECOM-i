      INTEGER FUNCTION FINDCTRD(CTRDB, STEPI, STEPJ)
*Compile by Sylvers Ding  
*Find the point(i+stepi,j+stepj) by following i direction before j
*or j direction before i.   ---- Notes by MaRin
      USE MOD_GLOBAL
      IMPLICIT NONE
      
      INTEGER CTRDB, STEPI, STEPJ
      INTEGER CTRDTEMP
      INTEGER FINDCTRDI, FINDCTRDJ
      
      CTRDTEMP = FINDCTRDI(CTRDB,STEPI)
      
      IF (CTRDTEMP.NE.N_CTRDP1) THEN
        FINDCTRD = FINDCTRDJ(CTRDTEMP,STEPJ)
      ELSE
        CTRDTEMP = FINDCTRDJ(CTRDB,STEPJ)
        IF (CTRDTEMP.NE.N_CTRDP1) THEN
          FINDCTRD = FINDCTRDI(CTRDTEMP,STEPI)
        ELSE
          FINDCTRD = 1
        ENDIF
      ENDIF
      
      RETURN
      
      END
      
      
      INTEGER FUNCTION FINDCTRDI(CTRDB,STEPI)
      
      USE MOD_GLOBAL
      
      IMPLICIT NONE
      
      INTEGER CTRDB, STEPI
      
      SELECT CASE (STEPI)
        
      CASE (-4)
        FINDCTRDI = NIM4(CTRDB)
      CASE (-3)
        FINDCTRDI = NIM3(CTRDB)
      CASE (-2)
        FINDCTRDI = NIM2(CTRDB)
      CASE (-1)
        FINDCTRDI = NIM1(CTRDB)
      CASE (0)
        FINDCTRDI = CTRDB
      CASE (1)
        FINDCTRDI = NIP1(CTRDB)
      CASE (2)
        FINDCTRDI = NIP2(CTRDB)
      CASE (3)
        FINDCTRDI = NIP3(CTRDB)
      CASE (4)
        FINDCTRDI = NIP4(CTRDB)
      CASE DEFAULT
        PRINT*, "WRONG STEPI IN FINDCTRDI!"
        PAUSE
        STOP
      END SELECT
      
      RETURN
      
      END
      
      
      INTEGER FUNCTION FINDCTRDJ(CTRDB,STEPJ)
      
      USE MOD_GLOBAL
      
      IMPLICIT NONE
      
      INTEGER CTRDB, STEPJ
      
      SELECT CASE (STEPJ)
        
      CASE (-4)
        FINDCTRDJ = NJM4(CTRDB)
      CASE (-3)      
        FINDCTRDJ = NJM3(CTRDB)
      CASE (-2)      
        FINDCTRDJ = NJM2(CTRDB)
      CASE (-1)      
        FINDCTRDJ = NJM1(CTRDB)
      CASE (0)      
        FINDCTRDJ = CTRDB
      CASE (1)       
        FINDCTRDJ = NJP1(CTRDB)
      CASE (2)       
        FINDCTRDJ = NJP2(CTRDB)
      CASE (3)       
        FINDCTRDJ = NJP3(CTRDB)
      CASE (4)       
        FINDCTRDJ = NJP4(CTRDB)
      CASE DEFAULT
        PRINT*, "WRONG STEPJ IN FINDCTRDJ!"
        PAUSE
        STOP
      END SELECT
      
      RETURN
      
      END