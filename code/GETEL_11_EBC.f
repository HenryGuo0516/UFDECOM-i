	SUBROUTINE  GETEL11EBC(KY,KKM,KDAY,TIMES,ADEPTH,NUMEBC,HA,GA,ELT)

      PARAMETER(PI=3.1415926)
      INTEGER NUMEBC
      INTEGER, PARAMETER :: MONTH(12)=
     * (/-1,30,58,89,119,150,180,211,242,272,303,333/)
	INTEGER  KC(11,6)
      
c	DIMENSION MONTH(12),PHASE(11),FREQ(11)
	DIMENSION PHASE(11),FREQ(11)
	DIMENSION A(11),B(11),FF(11),UU(11)
	DIMENSION HA(NUMEBC,11),GA(NUMEBC,11)

      REAL ELT(NUMEBC)
C
C	CHOOSE THE TIDAL CONSTITITE M2 S2 N2 K2 K1 O1 P1 Q1 M4 MS4 M6
C
c	DATA MONTH/-1,30,58,89,119,150,180,211,242,272,303,333/
      DATA ((KC(I,J),J=1,6),I=1,11)
     $       /2,0,0,0,0,0,2,2,-2,0,0,0,2,-1,0,1,0,0,2,2,0,0,0,0,
     $        1,1,0,0,0,1,1,-1,0,0,0,-1,1,1,-2,0,0,-1,1,-2,0,1,0,
     $        -1,4,0,0,0,0,0,4,2,-2,0,0,0,6,0,0,0,0,0/

	IF(ABS(TIMES-86400).LT.1.0) THEN
	KDAY=KDAY+1
	TIMES=0.0
	ENDIF

	NADDI=MONTH(KKM)+KDAY+INT((KY-1901)/4)
	IF(MOD(KY,4).EQ.0.AND.KKM.GT.2) NADDI=NADDI+1

	S=277.02+129.3848*(KY-1900)+13.1764*(NADDI)
	HH=280.19-0.2387*(KY-1900)+0.9857*(NADDI)
	P=334.39+40.6625*(KY-1900)+0.1114*(NADDI)
	XN=100.84+19.3282*(KY-1900)+0.053*(NADDI)
	PP=281.22+0.0172*(KY-1900)+0.00005*0.053*(NADDI)
	TAO=-S+HH
C
C	ATTENTION IN ALL COEFFICENT OF "XN" EQUAL TO 0
C
	DO 10 I=1,11
           PHASE(I)=KC(I,1)*TAO+KC(I,2)*S+KC(I,3)*HH+KC(I,4)*P
     $          +KC(I,5)*PP+KC(I,6)*90
10    CONTINUE
	DO 20,I=1,11
           FREQ(I)=14.49205211*KC(I,1)+0.54901653*KC(I,2)
     $          +0.04106864*KC(I,3)+0.00646183*KC(I,4)  
     $          +0.00000196*KC(I,5)
20    CONTINUE

        C=P*PI/180.0
        D=XN*PI/180.0
        A(1)=0.0005*COS(-2*D)-0.0373*COS(-D)+0.0006*COS(2*C)
     $        +0.0002*COS(2*C+D)+COS(0*D)
        A(4)=-0.0128*COS(-D)+0.2980*COS(D)+0.0324*COS(2*D)+COS(0*D)
        A(5)=0.0002*COS(-2*C-D)+0.0001*COS(-2*D)-0.0198*COS(-D)
     $        +0.1356*COS(D)-0.0029*COS(2*D)+COS(0*D)
        A(6)=-0.0058*COS(-2*D)+0.1885*COS(-D)+0.0002*COS(2*C-D)
     $       -0.0064*COS(2*C)-0.0010*COS(2*C+D)+COS(0*D)     
        A(7)=0.0008*COS(-2*D)-0.0112*COS(-D)-0.0015*COS(2*C)
     $          -0.0003*COS(2*C+D)+COS(0*D)
        B(1)=-0.0005*SIN(-2*D)-0.0373*SIN(-D)+0.0006*SIN(2*C)
     $        +0.0002*SIN(2*C+D)+SIN(0*D)
        B(4)=-0.0128*SIN(-D)+0.2980*SIN(D)+0.0324*SIN(2*D)+SIN(0*D)
        B(5)=0.0002*SIN(-2*C-D)+0.0001*SIN(-2*D)-0.0198*SIN(-D)
     $        +0.1356*SIN(D)-0.0029*SIN(2*D)+SIN(0*D)
        B(6)=-0.0058*SIN(-2*D)+0.1885*SIN(-D)+0.0002*SIN(2*C-D)
     $       -0.0064*SIN(2*C)-0.0010*SIN(2*C+D)+SIN(0*D)     
        B(7)=0.0008*SIN(-2*D)-0.0112*SIN(-D)-0.0015*SIN(2*C)
     $          -0.0003*SIN(2*C+D)+SIN(0*D)
        FF(1)=SQRT(A(1)**2+B(1)**2)
        FF(4)=SQRT(A(4)**2+B(4)**2)
        FF(5)=SQRT(A(5)**2+B(5)**2)
        FF(6)=SQRT(A(6)**2+B(6)**2)
        FF(7)=SQRT(A(7)**2+B(7)**2)
        UU(1)=ATAN(B(1)/A(1))
        UU(4)=ATAN(B(4)/A(4))
        UU(5)=ATAN(B(5)/A(5))
        UU(6)=ATAN(B(6)/A(6))
        UU(7)=ATAN(B(7)/A(7))
        FF(2)=1
        UU(2)=0
        FF(3)=FF(1)
        UU(3)=UU(1)
        FF(8)=FF(6)
        UU(8)=UU(6)
	  FF(9)=FF(1)**2
        UU(9)=UU(1)*2
        FF(10)=FF(1)
        UU(10)=UU(1)
        FF(11)=FF(1)**3
        UU(11)=UU(1)*3
        
	  ELT=ADEPTH
        
      DO N=1,NUMEBC !concurrentable when multicore
      DO 30 I=1,11
C
        ELT(N)=ELT(N)+FF(I)*HA(N,I)*COS((FREQ(I)*TIMES/3600.0+PHASE(I)
     $                -GA(N,I))*PI/180.0+UU(I))
30    ENDDO        
      ENDDO
c        ELT=ELT+ADEPTH
	
	RETURN
	END SUBROUTINE GETEL11EBC
