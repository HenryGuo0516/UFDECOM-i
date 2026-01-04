      SUBROUTINE FIELD_CHECK
  
      USE MOD_GLOBAL
  
      IMPLICIT NONE 
  
      INTEGER FID,I,NLAYER
      REAL BI1,BI2,BJ1,BJ2
      CHARACTER*20 TSTR

      NLAYER = 1 !THE LAYER TO BE VISUALIZED
      FID = 201
      BI1 = -586182.867627
      BI2 = 2748928.94
      BJ1 = 1392748.28404
      BJ2 = 4156064.6
  
      WRITE(TSTR,'(I6.6)') INT(THOUR)
  
      !H
      FN = '.\MODGEN'//TRIM(XG)//'FIELD_CHECK'//TRIM(XG)
     *//'FIELD_CHECK_'//TRIM(TSTR)//'_H.MIF'
      OPEN(FID,FILE=FN,STATUS='REPLACE')
      WRITE (FID,*) 'VERSION 300'
      WRITE (FID,*) 'CHARSET "WINDOWSSIMPCHINESE"'
      WRITE (FID,*) 'DELIMITER ","'
      WRITE (FID,2000) BI1,BI2,BJ1,BJ2
      WRITE (FID,*) 'COLUMNS 1'
      WRITE (FID,*) 'ID FLOAT'
      WRITE (FID,*) 'DATA'
      DO I=1,N_CTRD
        IF (FSM(I).EQ.1.) THEN
          WRITE (FID,3000) XR(I),YR(I)
          WRITE (FID,*) '    SYMBOL (35,128,12) '
        ENDIF
      ENDDO
      CLOSE(FID)
  
      FN = '.\MODGEN'//TRIM(XG)//'FIELD_CHECK'//TRIM(XG)
     *//'FIELD_CHECK_'//TRIM(TSTR)//'_H.MID'
      OPEN(FID,FILE=FN,STATUS='REPLACE')
      DO I=1,N_CTRD
        IF (FSM(I).EQ.1.) THEN
          WRITE(FID,*) H(I)
        ENDIF
      ENDDO
      CLOSE(FID)
  
      !SALINITY
      FN = '.\MODGEN'//TRIM(XG)//'FIELD_CHECK'//TRIM(XG)
     *//'FIELD_CHECK_'//TRIM(TSTR)//'_SAL.MIF'
      OPEN(FID,FILE=FN,STATUS='REPLACE')
      WRITE (FID,*) 'VERSION 300'
      WRITE (FID,*) 'CHARSET "WINDOWSSIMPCHINESE"'
      WRITE (FID,*) 'DELIMITER ","'
      WRITE (FID,2000) BI1,BI2,BJ1,BJ2
      WRITE (FID,*) 'COLUMNS 1'
      WRITE (FID,*) 'ID FLOAT'
      WRITE (FID,*) 'DATA'
      DO I=1,N_CTRD
        IF (FSM(I).EQ.1.) THEN
          WRITE (FID,3000) XR(I),YR(I)
          WRITE (FID,*) '    SYMBOL (35,128,12) '
        ENDIF
      ENDDO
      CLOSE(FID)
  
      FN = '.\MODGEN'//TRIM(XG)//'FIELD_CHECK'//TRIM(XG)
     *//'FIELD_CHECK_'//TRIM(TSTR)//'_SAL.MID'
      OPEN(FID,FILE=FN,STATUS='REPLACE')
      DO I=1,N_CTRD
        IF (FSM(I).EQ.1.) THEN
          WRITE(FID,*) S(I,NLAYER)
        ENDIF
      ENDDO
      CLOSE(FID)
  
      !TEMPERATURE
      FN = '.\MODGEN'//TRIM(XG)//'FIELD_CHECK'//TRIM(XG)
     *//'FIELD_CHECK_'//TRIM(TSTR)//'_TEMP.MIF'
      OPEN(FID,FILE=FN,STATUS='REPLACE')
      WRITE (FID,*) 'VERSION 300'
      WRITE (FID,*) 'CHARSET "WINDOWSSIMPCHINESE"'
      WRITE (FID,*) 'DELIMITER ","'
      WRITE (FID,2000) BI1,BI2,BJ1,BJ2
      WRITE (FID,*) 'COLUMNS 1'
      WRITE (FID,*) 'ID FLOAT'
      WRITE (FID,*) 'DATA'
      DO I=1,N_CTRD
        IF (FSM(I).EQ.1.) THEN
          WRITE (FID,3000) XR(I),YR(I)
          WRITE (FID,*) '    SYMBOL (35,128,12) '
        ENDIF
      ENDDO
      CLOSE(FID)
  
      FN = '.\MODGEN'//TRIM(XG)//'FIELD_CHECK'//TRIM(XG)
     *//'FIELD_CHECK_'//TRIM(TSTR)//'_TEMP.MID'
      OPEN(FID,FILE=FN,STATUS='REPLACE')
      DO I=1,N_CTRD
        IF (FSM(I).EQ.1.) THEN
          WRITE(FID,*) T(I,NLAYER)
        ENDIF
      ENDDO
      CLOSE(FID)

  
2000  FORMAT('COORDSYS NONEARTH UNITS "M" BOUNDS 
     * (',F16.6,',',F16.6,') (',F16.6,','F16.6,')')   
3000  FORMAT ('POINT   ',F14.5,'     ',F14.5)
     
      RETURN
  
      END SUBROUTINE FIELD_CHECK