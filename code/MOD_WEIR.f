      MODULE MOD_WEIR  !LTY
      
      USE MOD_GLOBAL
      
      SAVE
      CONTAINS
      
      
      SUBROUTINE ADDWEIR(IDX)
      USE MOD_GLOBAL
      INTEGER :: M,N,I,J
      REAL :: WH0,WP0
      REAL :: OVERFLOW1,SELFV,TOL
      REAL :: VC1
      INTEGER :: KB_TEMP
      INTEGER, INTENT(IN) :: IDX
      REAL :: XISHU,DU00,DV00,DELTA5
      
      GO TO (100,200,250), IDX
      
 100  CONTINUE

      OCSL=0.         !!!! ONLY CASULLI
      DO M=1,ND
        
      IF (TRIM(WEIR_TYP(M)).EQ.'J_DIRECT') THEN
        
        DO N=1,WEIR_N(M)
      
          I = WEIR_I(M,N)
          WH0 = -WEIR_H(M,N)
      
          IF ((DUM(I).EQ.0.)
     * .OR. (H(NIM1(I)).LE.0.5*WEIR_H(M,N))
     * .OR. (H(I).LE.0.5*WEIR_H(M,N)))THEN
            OCSL(M,N)=1.
            CYCLE
          ENDIF
      
          DELTA0(M,N) = ELF(NIM1(I))-WH0
          DELTA1(M,N) = ELF(I)-WH0
          DELTA5 = 0.5*(ELF(NIM1(I))+ELF(I))-WH0
      
          KB_WEIR_0(M,N) = CEILING(DELTA0(M,N)/D(NIM1(I))*KBM1)   !原来上个版本就是向上取整
          KB_WEIR_0(M,N) = MAX(KB_WEIR_0(M,N),0)
          !KB_WEIR_0(M,N)=MIN(KB_WEIR_0(M,N),5)
      
          KB_WEIR_1(M,N) = CEILING(DELTA1(M,N)/D(I)*KBM1)
          KB_WEIR_1(M,N) = MAX(KB_WEIR_1(M,N),0)
          !KB_WEIR_1(M,N)=MIN(KB_WEIR_1(M,N),5)
      
          DU00 = 0.5*(D(I)+D(NIM1(I)))
          KB_WEIR_5(M,N) = CEILING(DELTA5/DU00*KBM1)
          KB_WEIR_5(M,N) = MAX(KB_WEIR_5(M,N),0)
          !KB_WEIR_5(M,N)=MIN(KB_WEIR_5(M,N),5)
      
      
          IF((DUM(I).EQ.1.) 
     * .AND. (KB_WEIR_0(M,N)*KB_WEIR_1(M,N).GT.0)) THEN
            DUM(I) = 1.
          ELSE
            IF (OCSL(M,N).EQ.0.)THEN
              DUM(I) = 0.
            ENDIF
          ENDIF
      
        ENDDO   
      
      
      ELSE          !!!!! I DIRECTION
        DO N=1,WEIR_N(M)
          
          I = WEIR_I(M,N)
          WH0 = -WEIR_H(M,N)
      
          IF ((DVM(I).EQ.0.)
     * .OR.(H(NJM1(I)).LE.0.5*WEIR_H(M,N))
     * .OR.(H(I).LE.0.5*WEIR_H(M,N))) THEN
            OCSL(M,N) = 1.
            CYCLE
          ENDIF
      
          DELTA0(M,N) = ELF(NJM1(I))-WH0
          DELTA1(M,N) = ELF(I)-WH0
          DELTA5 = 0.5*(ELF(NJM1(I))+ELF(I))-WH0
      
          KB_WEIR_0(M,N) = CEILING(DELTA0(M,N)/D(NJM1(1))*KBM1)   !
          KB_WEIR_0(M,N) = MAX(KB_WEIR_0(M,N),0)
          !KB_WEIR_0(M,N) = MIN(KB_WEIR_0(M,N),5)
      
          KB_WEIR_1(M,N) = CEILING(DELTA1(M,N)/D(I)*KBM1)
          KB_WEIR_1(M,N) = MAX(KB_WEIR_1(M,N),0)
          !KB_WEIR_1(M,N) = MIN(KB_WEIR_1(M,N),5)
      
          DV00 = 0.5*(D(I)+D(NJM1(I)))
          KB_WEIR_5(M,N) = CEILING(DELTA5/DV00*KBM1)
          KB_WEIR_5(M,N) = MAX(KB_WEIR_5(M,N),0)
          !KB_WEIR_5(M,N)=MIN(KB_WEIR_5(M,N),5)
      
      
          IF((DVM(I).EQ.1.) 
     * .AND. (KB_WEIR_0(M,N)*KB_WEIR_1(M,N).GT.0)) THEN
            DVM(I) = 1.
          ELSE
            IF (OCSL(M,N).EQ.0.)THEN
              DVM(I) = 0.
            ENDIF
          ENDIF
      
        ENDDO   
      
      ENDIF
      
      ENDDO
      
      RETURN
      
      
 200  CONTINUE
      DO M=1,ND
      IF (TRIM(WEIR_TYP(M)).EQ.'J_DIRECT') THEN
        
      DO N=1,WEIR_N(M)
        I = WEIR_I(M,N)
      
        WH0 = -WEIR_H(M,N)
        WP0 = 0.5*(H(I)+H(NIM1(I)))+WH0
      
        IF (OCSL(M,N).EQ.1.) THEN
          CYCLE
        ENDIF
      
        IF(KB_WEIR_0(M,N).EQ.0 .AND. KB_WEIR_1(M,N).EQ.0) THEN
          UF(I,:) = 0.    
          XMFLUX(I,:) = 0.
          CYCLE !!! 两边水位低于堰高
        ENDIF
      
        IF(KB_WEIR_0(M,N).GT.0 .AND. KB_WEIR_1(M,N).EQ.0) THEN !!! 仅左面水位高于堰高
          CALL WEIR_FLOW(I,1,WH0,WP0,OVERFLOW1)
          OVERFLOW1 = OVERFLOW1/DTI
          !PRINT*,I,OVERFLOW1,1
          DO K = 1, KBM1
            IF(K.LE.KB_WEIR_0(M,N))THEN
              XMFLUX(I,K) = OVERFLOW1/DZ(K)/KB_WEIR_0(M,N)
            ELSE
              UF(I,K) = 0.    !应该加上吧
              XMFLUX(I,K) = 0.
            ENDIF
          ENDDO
        ENDIF
      
        IF(KB_WEIR_0(M,N).EQ.0 .AND. KB_WEIR_1(M,N).GT.0) THEN !!! 仅右面水位高于堰高
          CALL WEIR_FLOW(I,2,WH0,WP0,OVERFLOW1)
          OVERFLOW1 = OVERFLOW1/DTI
          !PRINT*,I,J,OVERFLOW1,2
          DO K = 1, KBM1
            IF(K.LE.KB_WEIR_1(M,N))THEN
              XMFLUX(I,K) = -OVERFLOW1/DZ(K)/KB_WEIR_1(M,N)
            ELSE
              UF(I,K) = 0.    !应该加上吧
              XMFLUX(I,K) = 0.
            ENDIF
          ENDDO
        ENDIF
      
        IF(KB_WEIR_0(M,N).GT.0 .AND. KB_WEIR_1(M,N).GT.0) THEN !!! 两边水位均高于堰高
          !KB_TEMP=MIN(KB_WEIR_0(M,N),KB_WEIR_1(M,N))      !!! 想想是取较小值好
          !KB_TEMP=MAX(KB_WEIR_0(M,N),KB_WEIR_1(M,N))      !!! 
          KB_TEMP = KB_WEIR_5(M,N)
          
          DO K = 1, KBM1
            IF(K.LE.KB_TEMP)THEN
            
              DU00 = 0.5*(D(I)+D(NIM1(I)))
              
              IF (K.LT.KB_TEMP)THEN
                XISHU = 1.
              ELSE
                XISHU = AMIN1( (1+Z(K))*DU00-WH0 , DU00*DZ(K) )
                XISHU = AMAX1(XISHU,0.)
                XISHU = XISHU/(DU00*DZ(K))
              ENDIF
            
              XMFLUX(I,K) = XMFLUX(I,K)*XISHU
            ELSE
              UF(I,K) = 0.
              XMFLUX(I,K) = 0.
            ENDIF
          ENDDO
        ENDIF
      
      ENDDO 
      
      ENDIF
      ENDDO
      
      RETURN
      
      
      
      
 250  CONTINUE
      DO M=1,ND
        IF (TRIM(WEIR_TYP(M)).EQ.'J_DIRECT') THEN
        ELSE        !! I_DIRECT
          DO N=1,WEIR_N(M)
            I = WEIR_I(M,N)
            
            WH0 = -WEIR_H(M,N)
            WP0 = 0.5*(H(I)+H(NJM1(I)))+WH0
            
            IF (OCSL(M,N).EQ.1.)THEN
              CYCLE
            ENDIF
            
            IF(KB_WEIR_0(M,N).EQ.0 .AND. KB_WEIR_1(M,N).EQ.0) THEN
              VF(I,:) = 0.    
              YMFLUX(I,:) = 0.
              CYCLE !!! 两边水位低于堰高
            ENDIF
            
            IF(KB_WEIR_0(M,N).GT.0 .AND. KB_WEIR_1(M,N).EQ.0) THEN !!! 仅下面水位高于堰高
              CALL WEIR_FLOW(I,3,WH0,WP0,OVERFLOW1)
              OVERFLOW1 = OVERFLOW1/DTI
              DO K = 1, KBM1
                IF(K.LE.KB_WEIR_0(M,N))THEN
                  !PRINT*,OVERFLOW1,K,3
                  YMFLUX(I,K) = OVERFLOW1/DZ(K)/KB_WEIR_0(M,N)
                ELSE
                  VF(I,K) = 0.    !应该加上吧
                  YMFLUX(I,K) = 0.
                ENDIF
              ENDDO
            ENDIF
            
            IF(KB_WEIR_0(M,N).EQ.0 .AND. KB_WEIR_1(M,N).GT.0) THEN !!! 仅上面水位高于堰高
              CALL WEIR_FLOW(I,4,WH0,WP0,OVERFLOW1)
              OVERFLOW1 = OVERFLOW1/DTI
              DO K = 1, KBM1
                IF(K.LE.KB_WEIR_1(M,N))THEN
                  !PRINT*,OVERFLOW1,K,4
                  YMFLUX(I,K) = -OVERFLOW1/DZ(K)/KB_WEIR_1(M,N)
                ELSE
                  VF(I,K) = 0.    !应该加上吧
                  YMFLUX(I,K) = 0.
                ENDIF
              ENDDO
            ENDIF
            
            IF(KB_WEIR_0(M,N).GT.0 .AND. KB_WEIR_1(M,N).GT.0) THEN !!! 两边水位均高于堰高
              !KB_TEMP=MIN(KB_WEIR_0(M,N),KB_WEIR_1(M,N))      !!! 想想是取较小值好
              !KB_TEMP=MAX(KB_WEIR_0(M,N),KB_WEIR_1(M,N))      !!! 
              KB_TEMP = KB_WEIR_5(M,N)
            
              DO K = 1, KBM1
                IF(K.LE.KB_TEMP)THEN
            
                  DV00 = 0.5*(D(I)+D(NJM1(I)))
            
                  IF (K.LT.KB_TEMP)THEN
                    XISHU = 1.
                  ELSE
                    XISHU = AMIN1( (1+Z(K))*DV00-WH0 , DV00*DZ(K) )
                    XISHU = AMAX1(XISHU,0.)
                    XISHU = XISHU/(DV00*DZ(K))
                  ENDIF
            
                  YMFLUX(I,K) = YMFLUX(I,K)*XISHU
                ELSE
                  VF(I,K) = 0.
                  YMFLUX(I,K) = 0.
                ENDIF
              ENDDO
            ENDIF
          
          ENDDO   
        ENDIF
      ENDDO
      
      RETURN
      
      END SUBROUTINE ADDWEIR
      
      
      
      
      SUBROUTINE ADDWEIR1(F)
      USE MOD_GLOBAL
      INTEGER :: M,N,I,J
      REAL :: WH0,WP0
      REAL :: OVERFLOW1,SELFV,TOL
      REAL :: VC1
      INTEGER :: KB_TEMP
      REAL :: XISHU,DU00,DV00,DELTA5
      
      REAL F(N_CTRDP1,KB)
      
      
      DO M=1,ND
        IF (TRIM(WEIR_TYP(M)).EQ.'J_DIRECT') THEN
          DO N=1,WEIR_N(M)
            I = WEIR_I(M,N)
      
            WH0 = -WEIR_H(M,N)
            WP0 = 0.5*(H(I)+H(NIM1(I)))+WH0
      
            IF (OCSL(M,N).EQ.1. .AND. DUM(I).NE.0.)THEN
              CYCLE
            ENDIF
      
            IF(KB_WEIR_0(M,N).EQ.0 .AND. KB_WEIR_1(M,N).EQ.0) THEN
              XFLUX(I,:) = 0.
            ENDIF
      
            IF(KB_WEIR_0(M,N).GT.0 .AND. KB_WEIR_1(M,N).EQ.0) THEN !!! 仅左面水位高于堰高
              DO K = 1, KBM1
                IF(K.LE.KB_WEIR_0(M,N))THEN                    
                  XFLUX(I,K) = F(NIM1(I),K)*XMFLUX(I,K)
                ELSE
                  XFLUX(I,K) = 0.
                ENDIF
              ENDDO
            ENDIF
      
            IF(KB_WEIR_0(M,N).EQ.0 .AND. KB_WEIR_1(M,N).GT.0) THEN !!! 仅右面水位高于堰高
              DO K = 1, KBM1
                IF(K.LE.KB_WEIR_1(M,N))THEN
                  XFLUX(I,K) = F(I,K)*XMFLUX(I,K)      
                ELSE
                  XFLUX(I,K) = 0.
                ENDIF
              ENDDO
            ENDIF
      
            IF(KB_WEIR_0(M,N).GT.0 .AND. KB_WEIR_1(M,N).GT.0) THEN !!! 两边水位均高于堰高
              !KB_TEMP=MIN(KB_WEIR_0(M,N),KB_WEIR_1(M,N))
              !KB_TEMP=MAX(KB_WEIR_0(M,N),KB_WEIR_1(M,N))      !!! 
              KB_TEMP = KB_WEIR_5(M,N)

              DO K = 1, KBM1
                IF(K.LE.KB_TEMP)THEN
      
                  IF(XMFLUX(I,K)>=0)THEN
                    XFLUX(I,K) = F(NIM1(I),K)*XMFLUX(I,K)
                  ELSE
                    XFLUX(I,K) = F(I,K)*XMFLUX(I,K)
                  ENDIF
      
                ELSE
                  XFLUX(I,K) = 0.              
                ENDIF
              ENDDO
            ENDIF
      
          ENDDO 
          
        ELSE        !! I_DIRECT
            
          DO N=1,WEIR_N(M)
            I = WEIR_I(M,N)
      
            WH0 = -WEIR_H(M,N)
            WP0 = 0.5*(H(I)+H(NJM1(I)))+WH0
      
      
            IF (OCSL(M,N).EQ.1. .AND. DVM(I).NE.0.)THEN
              CYCLE
            ENDIF
      
            IF(KB_WEIR_0(M,N).EQ.0 .AND. KB_WEIR_1(M,N).EQ.0) THEN
              YFLUX(I,:) = 0.
            ENDIF

            IF(KB_WEIR_0(M,N).GT.0 .AND. KB_WEIR_1(M,N).EQ.0) THEN !!! 仅下面水位高于堰高
              DO K = 1, KBM1
                IF(K.LE.KB_WEIR_0(M,N))THEN
                  YFLUX(I,K) = F(NJM1(I),K)*YMFLUX(I,K)
                ELSE
                  YFLUX(I,K) = 0.
                ENDIF
              ENDDO
            ENDIF
      
            IF(KB_WEIR_0(M,N).EQ.0 .AND. KB_WEIR_1(M,N).GT.0) THEN !!! 仅上面水位高于堰高
              DO K = 1, KBM1
                IF(K.LE.KB_WEIR_0(M,N))THEN
                  YFLUX(I,K) = F(I,K)*YMFLUX(I,K)
                ELSE
                  YFLUX(I,K) = 0.
                ENDIF
              ENDDO
            ENDIF
      
            IF(KB_WEIR_0(M,N).GT.0 .AND. KB_WEIR_1(M,N).GT.0) THEN !!! 两边水位均高于堰高
              !KB_TEMP=MIN(KB_WEIR_0(M,N),KB_WEIR_1(M,N))
              !KB_TEMP=MAX(KB_WEIR_0(M,N),KB_WEIR_1(M,N))      !!! 
              KB_TEMP = KB_WEIR_5(M,N)
      
              DO K = 1, KBM1
                IF(K.LE.KB_TEMP)THEN
                  IF(YMFLUX(I,K)>=0)THEN
                    YFLUX(I,K) = F(NJM1(I),K)*YMFLUX(I,K)
                  ELSE
                    YFLUX(I,K) = F(I,K)*YMFLUX(I,K)
                  ENDIF
      
                ELSE
                  YFLUX(I,K)=0.              
                ENDIF
              ENDDO
            ENDIF
      
          ENDDO   
        ENDIF
      ENDDO
      
      RETURN
      
      END SUBROUTINE ADDWEIR1
      
      
      
      SUBROUTINE WEIR_FLOW(I,TYPE0,WH0,WP0,OVERFLOW1)  !LTY

      USE MOD_GLOBAL
      REAL :: OVERFLOW1
      REAL :: Q_TOL
      REAL :: T_TOL
      REAL :: DTII
      REAL :: SEC_LEN,EL_HIGH,EL_LOW,S_HIGH,S_LOW
      REAL :: WH0,WP0,Q
      REAL :: M0,BH
      INTEGER :: TYPE0
          !!! TYPE0 = 1 :I-1,J --> I,J
          !!! TYPE0 = 2 :I,J --> I-1,J
          !!! TYPE0 = 3 :I,J-1 --> I,J
          !!! TYPE0 = 4 :I,J --> I,J-1
      INTEGER :: MM,MMM
      INTEGER :: I
            

      IF (TYPE0.EQ.1)THEN
        SEC_LEN = H2(I)
        EL_HIGH = EL(NIM1(I))
        EL_LOW = EL(I)
        S_HIGH = DJ(NIM1(I))
        S_LOW = DJ(I)
      ELSEIF (TYPE0.EQ.2)THEN
        SEC_LEN = H2(I)
        EL_LOW = EL(NIM1(I))
        EL_HIGH = EL(I)
        S_LOW = DJ(NIM1(I))
        S_HIGH = DJ(I)
      ELSEIF (TYPE0.EQ.3)THEN
        SEC_LEN = H1(I)
        EL_HIGH = EL(NJM1(I))
        EL_LOW = EL(I)
        S_HIGH = DJ(NJM1(I))
        S_LOW = DJ(I)
      ELSE
        SEC_LEN = H1(I)
        EL_LOW = EL(NJM1(I))
        EL_HIGH = EL(I)
        S_LOW = DJ(NJM1(I))
        S_HIGH = DJ(I)
      ENDIF

      MM = 500
      DTII = REAL(DTI/REAL(MM))     !LXY 20150507
            
      Q_TOL = 0.
      T_TOL = 0.
            
      DO MMM=1,MM
        IF(EL_HIGH.LE.WH0 .OR. EL_HIGH.LE.EL_LOW)THEN 
          GOTO 111
        ENDIF
        BH = EL_HIGH-WH0
                                
        M0 = 2./3.*(0.605+0.001/BH+0.08*BH/WP0) !LXY, T.REHBOCK
        M0 = AMIN1(M0,2.0)

        Q = SEC_LEN*M0*(SQRT(2.*9.806))*(BH**1.5)
        EL_HIGH = EL_HIGH-DTII*Q/S_HIGH
        EL_LOW = EL_LOW+DTII*Q/S_LOW
        Q_TOL = Q_TOL+Q*DTII
        T_TOL = T_TOL+DTII
      ENDDO
 111  CONTINUE
            
      OVERFLOW1 = Q_TOL  
      OVERFLOW1 = OVERFLOW1*0.8 !!!LXY FOR TEST

      RETURN 
      END SUBROUTINE WEIR_FLOW
      
      
      END MODULE MOD_WEIR
      