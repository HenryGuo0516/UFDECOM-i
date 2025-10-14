      SUBROUTINE PRE_TSR
	
!======================================================================
!  This program is used to prepare for the TSR interpolation process
!
!  by Sylvers Ding, Dec. 31, 2019	
!======================================================================
	
	USE MOD_GLOBAL
	
	IMPLICIT NONE
	
	INTEGER I
	
	NPOW = 1
	SEARCHRADIUS = MAX(MAXVAL(H1),MAXVAL(H2))*1.5
	
	
	IF (COORD .EQ. 'XY') THEN
        
      DO I = 1,N_TSR
        CALL WEIGHT_IDW
     * (X_TSR(I),Y_TSR(I),XR,YR,WTC_TSR(:,I),NAGC_TSR(:,I))
        CALL WEIGHT_IDW
     * (X_TSR(I),Y_TSR(I),XU,YU,WTU_TSR(:,I),NAGU_TSR(:,I))
        CALL WEIGHT_IDW
     * (X_TSR(I),Y_TSR(I),XV,YV,WTV_TSR(:,I),NAGV_TSR(:,I))
      ENDDO
      
      ELSE
        
      DO I = 1,N_TSR
        CALL WEIGHT_IDW
     * (X_TSR(I),Y_TSR(I),LON,LAT,WTC_TSR(:,I),NAGC_TSR(:,I))
        CALL WEIGHT_IDW
     * (X_TSR(I),Y_TSR(I),LONU,LATU,WTU_TSR(:,I),NAGU_TSR(:,I))
        CALL WEIGHT_IDW
     * (X_TSR(I),Y_TSR(I),LONV,LATV,WTV_TSR(:,I),NAGV_TSR(:,I))
      ENDDO
 
	ENDIF
	RETURN
	
	END SUBROUTINE PRE_TSR