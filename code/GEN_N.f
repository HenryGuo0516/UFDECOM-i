	SUBROUTINE GEN_N
  
      USE MOD_GLOBAL
      PI=3.1415926/180.
      IF (N_TSR.NE.0) THEN
      DO II = 1,N_TSR
          DPMIN = 1000000
          XX = X_TSR(II)
          YY = Y_TSR(II)
          DO I = 1,N_CTRD_AG
          IF (H(I).GE.-10.) THEN
          IF (COORD .EQ. 'XY') THEN
              DM = SQRT((XX-XR(I))**2+(YY-YR(I))**2)
          ELSE
              DX = (XX-LON(I))*PI*ER*COS(0.5*(YY+LAT(I))*PI)
              DY = (YY-LAT(I))*ER*PI
              DM = SQRT(DX**2+DY**2)
          ENDIF
          IF (DM.LE.DPMIN) THEN
              DPMIN = DM
              I0 = I
          ENDIF
          ENDIF
          ENDDO
          NN_TSR(II) = I0
      ENDDO
      ENDIF

! NXY_SEC
      !IF (N_SEC_XY.NE.0) THEN
      !DO II = 1,N_SEC_XY
      !    DO JJ = 1,NN_SEC(II)
      !    DPMIN = 1000000
      !    XX = XD_SEC(II,JJ)
      !    YY = YD_SEC(II,JJ)
      !    DO I = 1,N_CTRD_AG
      !    IF (H(I).GE.-10.) THEN
      !    IF (COORD .EQ. 'XY') THEN
      !        DM = SQRT((XX-XR(I))**2+(YY-YR(I))**2)
      !    ELSE
      !        DX = (XX-XR(I))*PI*ER*COS(0.5*(YY+YR(I))*PI)
      !        DY = (YY-YR(I))*ER*PI
      !        DM = SQRT(DX**2+DY**2)
      !    ENDIF
      !    IF (DM.LE.DPMIN) THEN
      !        DPMIN = DM
      !        I0 = I
      !    ENDIF
      !    ENDIF
      !    ENDDO
	     !   NXY_SEC(II,JJ) = I0
      !    ENDDO
      !ENDDO
      !ENDIF
      RETURN
	END SUBROUTINE GEN_N
