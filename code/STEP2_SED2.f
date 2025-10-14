      SUBROUTINE STEP2_SED2(S0)  
	
	USE MOD_GLOBAL
      
      REAL AF(N_CTRD_AG,KB),BF(N_CTRD_AG,KB),DF(N_CTRD_AG,KB)
      REAL GF(N_CTRD_AG,KB),TA(N_CTRD_AG,KB),TB(N_CTRD_AG,KB)
      REAL TC(N_CTRD_AG,KB),TD(N_CTRD_AG,KB),TE(N_CTRD_AG,KB)
	REAL TF(N_CTRD_AG,KB)
      REAL DWSBAR(N_CTRD_AG,KB),SSTAR(N_CTRD_AG,KB)
      REAL WSK(N_CTRD_AG,KB),LSK(N_CTRD_AG,KB)
      REAL VHF(N_CTRD_AG,KB),VHPF(N_CTRD_AG,KB),VHPF0(N_CTRD_AG,KB)
      REAL SS(N_CTRD_AG,KB)
      REAL KH0(N_CTRDP1,KB)
      REAL DD(N_CTRDP1)
      REAL FLUX0(N_CTRD_AG),FLUX1(N_CTRD_AG),FLUX2(N_CTRD_AG)
      REAL FENZI(N_CTRD_AG),FENMU(N_CTRD_AG)
      
	REAL S0(N_CTRDP1,KB)

      REAL EPSON
      REAL LL,LKM1,LK
      INTEGER KI
      
      EPSON = 0.00001
      
      AF = 0.
      BF = 0.
      DF = 0.
      GF = 0.
      TA = 0.
      TB = 0.
      TC = 0.
      TD = 0.
      TE = 0.
      TF = 0.
      DWSBAR = 0.
      WSK = 0.
      LSK = 0.
      SSTAR = 0.
      VHF = 0.
      VHPF = 0.
      VHPF0 = 0.
      SS = 0.
     
      DD = D
      KH0 = KH

	DO I = 1,N_CTRD
      IF (D(I)<=5.*DMIN) THEN
          DD(I) = 5.*DMIN
          KH0(I,:) = AMAX1(KH(I,:),0.0001)
      ELSE
	    KH0(I,:) = AMAX1(KH(I,:),0.00001)
      ENDIF
      ENDDO
      
      DO I = 1,N_CTRD_AG
      DO K = 2,KBM1
	    DWSBAR(I,K) = 0.5*(DWS(I,K-1)+DWS(I,K))
      ENDDO
      ENDDO

	
      DO I = 1,N_CTRD_AG
      DO K = 2, KBM1
          LL = DWSBAR(I,K)*DTI
          LKM1 = 0.5*D(I)*DZ(K-1) !LXY
          LK = 0.5*D(I)*DZ(K)
          LL = AMIN1(LL,LKM1)   !IMPORTANT        
          LSK(I,K) = (LKM1-LL)/(LKM1+LK)
      ENDDO
      ENDDO
      WSK = 0.  
	DO K = 2, KBM2
	  DO I = 1,N_CTRD_AG
          IF (FSM(I)*FSM11(I).EQ.1.) THEN
           TA(I,K) = DTI*DWSBAR(I,K)*(1.-LSK(I,K))/(DD(I)*DZ(K))
           TB(I,K) = DTI*DWSBAR(I,K)*LSK(I,K)/(DD(I)*DZ(K))  
           TC(I,K) = -DTI*KH0(I,K)/(DD(I)*DD(I)*DZ(K)*DZZ(K-1))
           TD(I,K) = -DTI*DWSBAR(I,K+1)/(DD(I)*DZ(K))*(1.-LSK(I,K+1))
           TE(I,K) = -DTI*DWSBAR(I,K+1)/(DD(I)*DZ(K))*LSK(I,K+1)
           TF(I,K) = -DTI*KH0(I,K+1)/(DD(I)*DD(I)*DZ(K)*DZZ(K))
           GF(I,K) = -(WSK(I,K)-WSK(I,K+1))/DZ(K) 
           AF(I,K) = -TE(I,K)+TF(I,K)
           BF(I,K) = 1.-TB(I,K)-TC(I,K)-TD(I,K)-TF(I,K)
           DF(I,K) = -TA(I,K)+TC(I,K)
           SSTAR(I,K) = S0(I,K)+DTI/DD(I)*GF(I,K)
		ENDIF
	  ENDDO
	ENDDO	
	DO I = 1,N_CTRD_AG
	  IF (FSM(I)*FSM11(I).EQ.1.) THEN
	    TA(I,1) = -DTI*DWSBAR(I,2)*(1.-LSK(I,2))/(DD(I)*DZ(1))
		TB(I,1) = -DTI*DWSBAR(I,2)*LSK(I,2)/(DD(I)*DZ(1)) !BUG FOUND !LSK(I,J,K)
		TC(I,1) = -DTI*KH0(I,2)/(DD(I)*DD(I)*DZ(1)*DZZ(1))
          AF(I,1) = -TB(I,1)+TC(I,1)
          BF(I,1) = 1.-TA(I,1)-TC(I,1)
          FLUX0(I) = QSI(I)
          FLUX1(I) = WSK(I,2)+FLUX0(I)
          FLUX2(I) = FLUX1(I)*DTI/(DD(I)*DZ(1))
          SSTAR(I,1) = S0(I,1)+FLUX2(I)
	  ENDIF
      ENDDO	
	DO I = 1,N_CTRD_AG
        IF (FSM(I)*FSM11(I).EQ.1.) THEN
          VHF(I,1) = -AF(I,1)/BF(I,1)
          VHPF(I,1) = SSTAR(I,1)/BF(I,1)
        ENDIF
      ENDDO
      
      DO K = 2, KBM2
        DO I = 1,N_CTRD_AG
		IF (FSM(I)*FSM11(I).EQ.1.) THEN
            VHPF0(I,K) = (BF(I,K)+DF(I,K)*VHF(I,K-1))
            VHF(I,K) = - AF(I,K)/VHPF0(I,K)
            VHPF(I,K) = (SSTAR(I,K)-DF(I,K)*VHPF(I,K-1))/VHPF0(I,K)
		ENDIF
	  ENDDO
	ENDDO
	 
      DO K = 1, KBM1
        DO I = 1,N_CTRD_AG
          IF (FSM(I)*FSM11(I).EQ.1.) THEN
            SS(I,K) = S0(I,K)
          ENDIF
	  ENDDO
	ENDDO
      DO I = 1,N_CTRD_AG
	   IF (FSM(I)*FSM11(I).EQ.1.) THEN
		 TA(I,KBM1) = DTI*DWSBAR(I,KBM1)*(1.-LSK(I,KBM1))/(DD(I)*DZ(KBM1))
		 TB(I,KBM1) = DTI*DWSBAR(I,KBM1)*LSK(I,KBM1)/(DD(I)*DZ(KBM1))
		 TC(I,KBM1) = -DTI*KH0(I,KBM1)/(DD(I)*DD(I)*DZ(KBM1)*DZZ(KBM2)) 
           BF(I,KBM1) = 1.-TB(I,KBM1)-TC(I,KBM1)
           DF(I,KBM1) = -TA(I,KBM1)+TC(I,KBM1)
           FLUX0(I) = QDEP(I)-QERO(I)
           FLUX1(I) = WSK(I,KBM1)+FLUX0(I)
           FLUX2(I) = FLUX1(I)*DTI/(DD(I)*DZ(KBM1))
           SSTAR(I,KBM1) = S0(I,KBM1)-FLUX2(I)
	  ENDIF                             
      ENDDO
      DO I = 1,N_CTRD_AG
        IF (FSM(I)*FSM11(I).EQ.1.) THEN
          FENZI(I) = SSTAR(I,KBM1)-DF(I,KBM1)*VHPF(I,KBM2)
          FENMU(I) = BF(I,KBM1)+DF(I,KBM1)*VHF(I,KBM2)
          SS(I,KBM1) = FENZI(I)/FENMU(I)
          SS(I,KBM1) = AMIN1(SS(I,KBM1),12.0)
        ENDIF                             
      ENDDO

      DO K = 2, KBM1
        KI = KB - K
        DO I = 1,N_CTRD_AG
          IF (FSM(I)*FSM11(I).EQ.1.) THEN
            SS(I,KI) = (VHF(I,KI)*SS(I,KI+1)+VHPF(I,KI))
          ENDIF
        ENDDO
      ENDDO
      
	
      DO K = 1,KBM1
        DO I = 1,N_CTRD_AG
          IF (FSM(I)*FSM11(I).EQ.1.) THEN
            S0(I,K) = SS(I,K)
            S0(I,K) = AMAX1(S0(I,K), 0.) 
		  !因通量修正精度有限，泥沙还是会出现很小的负值，强行限制
          ENDIF
        ENDDO
      ENDDO

        
      RETURN
      
      END SUBROUTINE STEP2_SED2