      SUBROUTINE VGA_CHK(VA,X,Y,MASK,VARN)

      USE MOD_GLOBAL

      IMPLICIT NONE


      INTEGER FID,I,NLAYER
      REAL BI1,BI2,BJ1,BJ2
      REAL VA(N_CTRD),X(N_CTRD),Y(N_CTRD),MASK(N_CTRD)
      CHARACTER*10 VARN
      CHARACTER*20 TSTR

      FID = 202
      !BI1 = -586182.867627
      !BI2 = 2748928.94
      !BJ1 = 1392748.28404
      !BJ2 = 4156064.6
      BI1 = MINVAL(X)-100000
      BI2 = MINVAL(Y)-100000
      BJ1 = MAXVAL(X)+100000
      BJ2 = MAXVAL(Y)+100000

      WRITE(TSTR,'(I6.6)') INT(THOUR)

      !VG
      FN = '.\MODGEN'//TRIM(XG)//'VGA_CHECK'//TRIM(XG)//TRIM(TSTR)
     *//'_'//TRIM(VARN)//'_VG.MIF'
      OPEN(FID,FILE=FN,STATUS='REPLACE')
      WRITE (FID,*) 'VERSION 300'
      WRITE (FID,*) 'CHARSET "WINDOWSSIMPCHINESE"'
      WRITE (FID,*) 'DELIMITER ","'
      WRITE (FID,2000) BI1,BI2,BJ1,BJ2
      WRITE (FID,*) 'COLUMNS 1'
      WRITE (FID,*) 'ID FLOAT'
      WRITE (FID,*) 'DATA'
  
      DO I = N_CTRD_AG+1,N_CTRD
        IF (MASK(I).EQ.1.) THEN
          WRITE (FID,3000) X(I),Y(I)
          WRITE (FID,*) '    SYMBOL (35,128,12) '
        ENDIF
      ENDDO
      CLOSE(FID)
  
      FN = '.\MODGEN'//TRIM(XG)//'VGA_CHECK'//TRIM(XG)//TRIM(TSTR)
     *//'_'//TRIM(VARN)//'_VG.MID'
      OPEN(FID,FILE=FN,STATUS='REPLACE')
      DO I = N_CTRD_AG+1,N_CTRD
        IF (MASK(I).EQ.1.) THEN
          WRITE(FID,*) VA(I)
        ENDIF
      ENDDO
      CLOSE(FID)


      !AG
      FN = '.\MODGEN'//TRIM(XG)//'VGA_CHECK'//TRIM(XG)//TRIM(TSTR)
     *//'_'//TRIM(VARN)//'_AG.MIF'
      OPEN(FID,FILE=FN,STATUS='REPLACE')
      WRITE (FID,*) 'VERSION 300'
      WRITE (FID,*) 'CHARSET "WINDOWSSIMPCHINESE"'
      WRITE (FID,*) 'DELIMITER ","'
      WRITE (FID,2000) BI1,BI2,BJ1,BJ2
      WRITE (FID,*) 'COLUMNS 1'
      WRITE (FID,*) 'ID FLOAT'
      WRITE (FID,*) 'DATA'
  
      DO I = 1,N_CTRD_AG
        IF (MASK(I).EQ.1.) THEN
          WRITE (FID,3000) X(I),Y(I)
          WRITE (FID,*) '    SYMBOL (35,128,12) '
        ENDIF
      ENDDO
      CLOSE(FID)
  
      FN = '.\MODGEN'//TRIM(XG)//'VGA_CHECK'//TRIM(XG)//TRIM(TSTR)
     *//'_'//TRIM(VARN)//'_AG.MID'
      OPEN(FID,FILE=FN,STATUS='REPLACE')
      DO I = 1,N_CTRD_AG
        IF (MASK(I).EQ.1.) THEN
          WRITE(FID,*) VA(I)
        ENDIF
      ENDDO
      CLOSE(FID)


2000  FORMAT('COORDSYS NONEARTH UNITS "M" BOUNDS 
     * (',F16.6,',',F16.6,') (',F16.6,','F16.6,')')   
3000  FORMAT ('POINT   ',F14.5,'     ',F14.5)

      RETURN

      END SUBROUTINE VGA_CHK