	SUBROUTINE BROUGH
	USE MOD_GLOBAL
      
      DO I = 1, N_CTRD
      IF (DU(I).GE.3.0) THEN  
          Z0 = Z0B  
          CBCMIN = BFRIC 
          CBC_U(I) = AMAX1(CBCMIN,.16/
     *    ALOG((ZZ(KBM1)-Z(KB))*DU(I)/Z0)**2) 
      ELSEIF (DU(I).GE.1..AND.DU(I).LT.3.0) THEN
          Z0 = Z0B  
          CBCMIN = BFRIC    
          CBC_U(I) = AMAX1(CBCMIN,.16/                 
     *        ALOG((ZZ(KBM1)-Z(KB))*3.0/Z0)**2) 
      ELSEIF (HU(I).GE.-10..AND.DU(I).LT.1.0) THEN
          Z0 = Z0B  
          CBCMIN = BFRIC    
          CBC_U(I) = AMAX1(CBCMIN,.16/                 
     *    ALOG((ZZ(KBM1)-Z(KB))*1.0/Z0)**2)      
          C1 = 0.005
          C2 = CBC_U(I)
          CBC_U(I) = (DU(I)-0.)*C2+(1.-DU(I))*C1
      ENDIF
      ENDDO
      
      DO I = 1, N_CTRD
      IF (DV(I).GE.3.0) THEN  
          Z0 = Z0B  
          CBCMIN = BFRIC 
          CBC_V(I) = AMAX1(CBCMIN,.16/
     *    ALOG((ZZ(KBM1)-Z(KB))*DV(I)/Z0)**2) 
      ELSEIF (DV(I).GE.1..AND.DV(I).LT.3.0) THEN
          Z0 = Z0B  
          CBCMIN = BFRIC    
          CBC_V(I) = AMAX1(CBCMIN,.16/                 
     *    ALOG((ZZ(KBM1)-Z(KB))*3.0/Z0)**2) 
      ELSEIF (HV(I).GE.-10..AND.DV(I).LT.1.0) THEN
          Z0 = Z0B  
          CBCMIN = BFRIC    
          CBC_V(I) = AMAX1(CBCMIN,.16/                 
     *    ALOG((ZZ(KBM1)-Z(KB))*1.0/Z0)**2)      
          C1=0.005
          C2=CBC_V(I)
          CBC_V(I)=(DV(I)-0.)*C2+(1.-DV(I))*C1
      ENDIF
      ENDDO
         
      DO I = 1, N_CTRD
      IF (D(I).GE.3.0) THEN  
          Z0 = Z0B  
          CBCMIN = BFRIC 
          CBC_UV(I) = AMAX1(CBCMIN,.16/
     *    ALOG((ZZ(KBM1)-Z(KB))*D(I)/Z0)**2) 
      ELSEIF (D(I).GE.1..AND.D(I).LT.3.0) THEN
          Z0 = Z0B  
          CBCMIN = BFRIC    
          CBC_UV(I) = AMAX1(CBCMIN,.16/                 
     *    ALOG((ZZ(KBM1)-Z(KB))*3.0/Z0)**2) 
      ELSEIF (H(I).GE.-10..AND.D(I).LT.1.0) THEN
          Z0 = Z0B  
          CBCMIN = BFRIC    
          CBC_UV(I) = AMAX1(CBCMIN,.16/                 
     *    ALOG((ZZ(KBM1)-Z(KB))*1.0/Z0)**2)      
          C1=0.005
          C2=CBC_UV(I)
          CBC_UV(I)=(D(I)-0.)*C2+(1.-D(I))*C1
      ENDIF
      ENDDO
   	RETURN
	END SUBROUTINE BROUGH


