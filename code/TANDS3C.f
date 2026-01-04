      SUBROUTINE TANDS3C
C     VERSION(7/27/98)
C COMMON:
C     THIS IS FOR BAROPG3 AND BAROPG4
C     S AND T IS GIVEN IN SOME STANDARD LEVEL
C     WE CALCULATE THE RMEAN IN U AND V GRID
C ANOTHER COMMON: T OF CHANGJIANG RIVER MOUTH IS BE VERIFIED

      USE MOD_GLOBAL
      INTEGER :: KSL1
      REAL, ALLOCATABLE :: RHZ(:)
      REAL, ALLOCATABLE :: TZZ(:)
      REAL, ALLOCATABLE :: SZZ(:)
      REAL, ALLOCATABLE :: DEP(:)
      REAL, ALLOCATABLE :: TZ(:,:)      
      REAL, ALLOCATABLE :: SZ(:,:) 
      DIMENSION RHS(KB),TI(KB),SI(KB)
      DIMENSION DEPS(KB),TA(KSL),SA(KSL)
      KSL1 = KSL+1
      ALLOCATE(RHZ(KSL1),TZZ(KSL1),SZZ(KSL1),DEP(KSL1))
      ALLOCATE(TZ(N_CTRD,KSL1),SZ(N_CTRD,KSL1))
      DI = 200.0
      
      READ (IUITS,*)
      DO I = 1,N_CTRD
	    READ (IUITS,*) II,(TZ(II,K),K=1,KSL1)
	    READ (IUITS,*) II,(SZ(II,K),K=1,KSL1)
	ENDDO
	CLOSE (IUITS)
    
      DO I=1,N_CTRD
      IF (H(I).GE.-10.AND.H(I).LT.DMIN) THEN
          DO K = 1,KBM1
	      S(I,K) = SZ(I,1)
	      T(I,K) = TZ(I,1)
          ENDDO
      ELSEIF(H(I).GT.DMIN.AND.H(I).LT.DI) THEN
          DO K=1,KSL
          IF(H(I).GE.-DPTHSL(K)) THEN
              DEP(K) = DPTHSL(K)
              TZZ(K) = TZ(I,K)
              SZZ(K) = SZ(I,K)
              T0ZZ(I,K) = TZ(I,K)
              S0ZZ(I,K) = SZ(I,K)
              K1 = K
              IF(H(I).EQ.-DPTHSL(K)) GOTO 15
          ELSE
              K1 = K
              GOTO 10
          ENDIF
          ENDDO
10        DEP(K1) = -H(I)
          TZZ(K1) = TZ(I,KSL1)
          SZZ(K1) = SZ(I,KSL1)
          T0ZZ(I,K1) = TZ(I,KSL1)
          S0ZZ(I,K1) = SZ(I,KSL1)
15        CONTINUE
          DO K = 1,KBM1
              DEPS(K) = ZZ(K)*H(I)
          ENDDO
          CALL SINTER(DEP,TZZ,DEPS,RHS,K1,KBM1)
          DO K = 1,KBM1
                T(I,K) = RHS(K)
          ENDDO
          CALL SINTER(DEP,SZZ,DEPS,RHS,K1,KBM1)
          DO K = 1,KBM1
              S(I,K) = RHS(K)
          ENDDO
      ELSEIF (H(I).GT.DI) THEN
          DO K = 1,KSL
          IF(H(I).GE.-DPTHSL(K)) THEN
              DEP(K) = DPTHSL(K)
              TZZ(K) = TZ(I,K)
              SZZ(K) = SZ(I,K)
              T0ZZ(I,K) = TZ(I,K)
              S0ZZ(I,K) = SZ(I,K)
              K1 = K
              IF(H(I).EQ.-DPTHSL(K)) GOTO 20
          ELSE
              IF(H(I).GE.-DPTHSL(KSL-1)) THEN
                  T0ZZ(I,K) = TZ(I,KSL-1)
                  S0ZZ(I,K) = SZ(I,KSL-1)
              ELSE
                  T0ZZ(I,K) = TZ(I,K-1)+(TZ(I,K-1)-TZ(I,K-2))*
     *            (H(I)+DPTHSL(K-1))/(-DPTHSL(K-1)+DPTHSL(K-2))
                  S0ZZ(I,K) = SZ(I,K-1)+(SZ(I,K-1)-SZ(I,K-2))*
     *            (H(I)+DPTHSL(K-1))/(-DPTHSL(K-1)+DPTHSL(K-2))
              ENDIF
              GOTO 20
          ENDIF
          ENDDO
20        CONTINUE
          DO K = 1,KBM1
              DEPS(K) = ZZ(K)*H(I)
          ENDDO
          CALL SINTER(DEP,TZZ,DEPS,RHS,K1,KBM1)
          DO K = 1,KBM1
              T(I,K) = RHS(K)
          ENDDO
          CALL SINTER(DEP,SZZ,DEPS,RHS,K1,KBM1)
          DO K = 1,KBM1
              S(I,K) = RHS(K)
          ENDDO
          DO K = 1,KBM1
          IF(-DEPS(K).GT.-DPTHSL(KSL-1)) THEN
              T(I,K) = TZZ(KSL-1)
              S(I,K) = SZZ(KSL-1)
          ENDIF
          ENDDO
      ENDIF
      ENDDO
      
      CALL DENS
      RETURN
      
      END SUBROUTINE TANDS3C

