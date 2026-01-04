#include "DEFS.h"
      SUBROUTINE ADVU_TVD_3RD
      USE MOD_GLOBAL
      
c      INTEGER BD(3)
      INTEGER IC,INUM,IE,K
      REAL X1,X2,X3,Y1,Y2,Y3,AVEL,BVEL,CVEL
      REAL F1,F2,F3,F4,CNUM1,CNUM
	REAL CFLX,CFLY,CFL,ISIGN,ISIGN1,CFLX1,RAT
	REAL DUMT,DUMT1,DVMT,DVMT1,HUP,ELUP,HDOWN,ELDOWN,DTN

c      UF=0.
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO I = 1, N_CTRDP1
c      DO K = 1, KB
      DO concurrent (I=1:N_CTRDP1,K=1:KB)
          UF(I,K)=0.
          CURV4(I,K)=0.
          XFLUX(I,K)=0.
          YFLUX(I,K)=0.
      ENDDO
c      ENDDO
c!$OMP END TARGET
      
c#ifdef OMP
c!$OMP PARALLEL DO PRIVATE(CFLX,CFLY,CFLX1,CFLY1
c     *    ,FUP,FLW,FS,FUP1,FLW1,FS1,R,RATIO,K)
c#endif
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO I = 1,N_CTRD
c      DO K = 1,KBM1   
      DO concurrent (I=1:N_CTRD,K=1:KBM1) 
     * local(CFLX,FUP,CFLX1,FUP1,FLW1,FLW,FS,FS1,R,RATIO,
     *       CFLY,CFLY1) local_init(DTI)
          
          CURV4(I,K) = .125*((V(NJP1(I),K)+V(I,K))
     *    *((DJ(NIP1(I))/H1(NIP1(I))*FSM(NIP1(I))+DJ(I)/H1(I)*FSM(I))
     *    /(FSM(NIP1(I))+FSM(I)+1.E-30)
     *    -(DJ(I)/H1(I)*FSM(I)+DJ(NIM1(I))/H1(NIM1(I))*FSM(NIM1(I)))
     *    /(FSM(I)+FSM(NIM1(I))+1.E-30))
     *    -(U(NIP1(I),K)+U(I,K))
     *    *((DJ(NJP1(I))/H2(NJP1(I))*FSM(NJP1(I))+DJ(I)/H2(I)*FSM(I))
     *    /(FSM(NJP1(I))+FSM(I)+1.E-30)
     *    -(DJ(I)/H2(I)*FSM(I)+DJ(NJM1(I))/H2(NJM1(I))*FSM(NJM1(I)))
     *    /(FSM(I)+FSM(NJM1(I))+1.E-30)))/DJ(I)+.25*COR(I) 
          
          CFLX = 0.5*(U(NIP1(I),K)+U(I,K))*DTI/H1(I)
          IF (CFLX.GE.0.) THEN
              FUP = XMFLUX(I,K)*U(I,K)
              CFLX1 = 0.5*(U(I,K)+U(NIM1(I),K))*DTI/H1(NIM1(I))
              IF (CFLX1.GE.0.) THEN
                  FUP1 = XMFLUX(NIM1(I),K)*U(NIM1(I),K)
              ELSE
                  FUP1 = XMFLUX(I,K)*U(I,K)
              ENDIF
              FLW1 = 0.5*((1+CFLX1)*XMFLUX(NIM1(I),K)*U(NIM1(I),K)
     *        +(1.-CFLX1)*XMFLUX(I,K)*U(I,K))
          ELSE
              FUP = XMFLUX(NIP1(I),K)*U(NIP1(I),K)
              CFLX1 = 0.5*(U(NIP1(I),K)+U(NIP2(I),K))*DTI/H1(NIP1(I))
              IF (CFLX1.GE.0.) THEN
                  FUP1 = XMFLUX(NIP1(I),K)*U(NIP1(I),K)
              ELSE
                  FUP1 = XMFLUX(NIP2(I),K)*U(NIP2(I),K)
              ENDIF
              FLW1 = 0.5*((1+CFLX1)*XMFLUX(NIP1(I),K)*U(NIP1(I),K)
     *        +(1.-CFLX1)*XMFLUX(NIP2(I),K)*U(NIP2(I),K))
          ENDIF
          FLW = 0.5*((1+CFLX)*XMFLUX(I,K)*U(I,K)
     *    +(1.-CFLX)*XMFLUX(NIP1(I),K)*U(NIP1(I),K))
          FS = FLW-FUP
          FS1 = FLW1-FUP1
          IF (ABS(FS).LE.1.E-10) THEN
              R = 0.
          ELSE
              R = FS1/FS
          ENDIF
c          RATIO = RAT(R)
            IF (R.LE.0.) THEN
                RATIO=0.
            ELSE
                RATIO=(R+R)/(1.+R)
            ENDIF

          XFLUX(I,K) = FUP+RATIO*FS 
          
          CFLY = 0.5*(V(I,K)+V(NIM1(I),K))*DTI
     *    /(0.25*(H2(I)+H2(NJM1(I))+H2(NIM1(I))+H2(NIM1JM1(I))))
          IF (CFLY.GE.0.) THEN
              FUP = U(NJM1(I),K)
              CFLY1 = 0.5*(V(NJM1(I),K)+V(NIM1JM1(I),K))*DTI/(0.25*
     *        (H2(NJM1(I))+H2(NJM2(I))+H2(NIM1JM1(I))+H2(NIM1JM2(I))))
              IF (CFLY1.GE.0.) THEN
                  FUP1 = U(NJM2(I),K)
              ELSE
                  FUP1 = U(NJM1(I),K)
              ENDIF
              FLW1=0.5*((1+CFLY1)*U(NJM2(I),K)+(1.-CFLY1)*U(NJM1(I),K))
          ELSE
              FUP = U(I,K)
              CFLY1 = 0.5*(V(NJP1(I),K)+V(NIM1JP1(I),K))*DTI
     *        /(0.25*(H2(I)+H2(NJP1(I))+H2(NIM1(I))+H2(NIM1JP1(I))))
              IF (CFLY1.GE.0.) THEN
                  FUP1 = U(I,K)
              ELSE
                  FUP1 = U(NJP1(I),K)
              ENDIF
              FLW1=0.5*((1+CFLY1)*U(I,K)+(1.-CFLY1)*U(NJP1(I),K))
          ENDIF
          FLW = 0.5*((1+CFLY)*U(NJM1(I),K)+(1.-CFLY)*U(I,K))
          FS = FLW-FUP
          FS1 = FLW1-FUP1     
          IF (ABS(FS).LE.1.E-10) THEN
              R = 0.
          ELSE
              R = FS1/FS
          ENDIF
c          RATIO = RAT(R)
            IF (R.LE.0.) THEN
                RATIO=0.
            ELSE
                RATIO=(R+R)/(1.+R)
            ENDIF
          YFLUX(I,K)=(FUP+RATIO*FS)*0.5*(YMFLUX(I,K)+YMFLUX(NIM1(I),K))
      ENDDO    
c      ENDDO
c!$OMP END TARGET
c#ifdef OMP
c!$OMP END PARALLEL DO
c#endif
      

#if defined TIDE_EL || defined TIDE_FLUX
      IF (NUMEBC.NE.0) THEN
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO N = 1,NUMEBC
c      DO K = 1,KBM1
      DO concurrent (N=1:NUMEBC,K=1:KBM1) local(IE,IC,IA) 
          IE = NETA(N)
          IC = NCON(N)
          IA = NAJA(N)
              XFLUX(IA,K)=.25*(XMFLUX(NIP1(IA),K)+XMFLUX(IA,K))*
     *        (U(NIP1(IA),K)+U(IA,K))
              YFLUX(IA,K)=.25*(YMFLUX(IA,K)+YMFLUX(NIM1(IA),K))*
     *        (U(IA,K)+U(NJM1(IA),K))
              XFLUX(IC,K)=.25*(XMFLUX(NIP1(IC),K)+XMFLUX(IC,K))*
     *        (U(NIP1(IC),K)+U(IC,K))
              YFLUX(IC,K)=.25*(YMFLUX(IC,K)+YMFLUX(NIM1(IC),K))*
     *        (U(IC,K)+U(NJM1(IC),K))
              XFLUX(IE,K)=.25*(XMFLUX(NIP1(IE),K)+XMFLUX(IE,K))*
     *        (U(NIP1(IE),K)+U(IE,K))
              YFLUX(IE,K)=.25*(YMFLUX(IE,K)+YMFLUX(NIM1(IE),K))*
     *        (U(IE,K)+U(NJM1(IE),K))
      ENDDO
c      ENDDO
c!$OMP END TARGET
      ENDIF
#elif defined TIDE_FLATHER
      IF (NUMEBC.NE.0) THEN
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO N = 1,NUMEBC
c      DO K = 1,KBM1
      DO concurrent (N=1:NUMEBC,K=1:KBM1) local(IE,IC,IA) 
          IE = NETA(N)
          IC = NCON(N)
          IA = NAJA(N)
              XFLUX(IA,K)=.25*(XMFLUX(NIP1(IA),K)+XMFLUX(IA,K))*
     *        (U(NIP1(IA),K)+U(IA,K))
              YFLUX(IA,K)=.25*(YMFLUX(IA,K)+YMFLUX(NIM1(IA),K))*
     *        (U(IA,K)+U(NJM1(IA),K))
              XFLUX(IC,K)=.25*(XMFLUX(NIP1(IC),K)+XMFLUX(IC,K))*
     *        (U(NIP1(IC),K)+U(IC,K))
              YFLUX(IC,K)=.25*(YMFLUX(IC,K)+YMFLUX(NIM1(IC),K))*
     *        (U(IC,K)+U(NJM1(IC),K))
              XFLUX(IE,K)=.25*(XMFLUX(NIP1(IE),K)+XMFLUX(IE,K))*
     *        (U(NIP1(IE),K)+U(IE,K))
              YFLUX(IE,K)=.25*(YMFLUX(IE,K)+YMFLUX(NIM1(IE),K))*
     *        (U(IE,K)+U(NJM1(IE),K))
      ENDDO
c      ENDDO
c!$OMP END TARGET
      ENDIF
#endif


      
c#ifdef OMP
c!$OMP PARALLEL DO
c#endif
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(PROD0)
c      DO I = 1,N_CTRD
c      DO K = 1,KBM1     
      DO concurrent (I=1:N_CTRD,K=1:KBM1) local(PROD0) 
          PROD0 = 4.*(D(I)*FSM(I)+D(NIM1(I))*FSM(NIM1(I))
     *    +D(NJM1(I))*FSM(NJM1(I))+D(NIM1JM1(I))*FSM(NIM1JM1(I)))
     *    /(FSM(I)+FSM(NIM1(I))+FSM(NJM1(I))+FSM(NIM1JM1(I))+1.E-5)
     *    *(AAM(I,K)*FSM(I)+AAM(NIM1(I),K)*FSM(NIM1(I))
     *    +AAM(NJM1(I),K)*FSM(NJM1(I))
     *    +AAM(NIM1JM1(I),K)*FSM(NIM1JM1(I)))
     *    /(FSM(I)+FSM(NIM1(I))+FSM(NJM1(I))+FSM(NIM1JM1(I)))
          XFLUX(I,K) = XFLUX(I,K)-
     *    D(I)*AAM(I,K)*FSM(I)*2.*(U(NIP1(I),K)-U(I,K))*H2(I)/H1(I)
          YFLUX(I,K) = YFLUX(I,K)
     *    -PROD0 *((U(I,K)-U(NJM1(I),K))
     *    /(H2(I)+H2(NIM1(I))+H2(NJM1(I))+H2(NIM1JM1(I)))
     *    +(V(I,K)-V(NIM1(I),K))
     *    /(H1(I)+H1(NIM1(I))+H1(NJM1(I))+H1(NIM1JM1(I)))) 
     *    *.25*(H1(I)+H1(NIM1(I))+H1(NJM1(I))+H1(NIM1JM1(I)))
     *    *DUM(I)*DUM(NJM1(I))
      ENDDO    
c      ENDDO
c!$OMP END TARGET
c#ifdef OMP
c!$OMP END PARALLEL DO
c#endif


***********************************************************************
*                             VG BONDOURY                             *

c!$OMP TARGET teams DEFAULTMAP(present: allocatable)
c!$OMP DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K=1,KBM1
c      DO I=1,NNIVGFWALL
      DO concurrent (I=1:NNIVGFWALL,K=1:KBM1) 
     * local(INUM,IC,X1,X3,Y1,Y2,Y3,BVEL,AVEL,CVEL) 
      INUM=NIVGFWALL(I,1)
      IC =NIVGFWALL(I,2)
      IF (U(NIP1(INUM),K).GT.0.) THEN
      IF (XFLUX(IC,K).EQ.0) THEN
          XFLUX(INUM,K)=0.
      ELSE
          X1=(YR(NJM1(IC))-YR(IC))/H1(INUM)
c          X2=0.
          X3=(YR(NJP1(IC))-YR(IC))/H1(INUM)
          Y1=XFLUX(NJM1(IC),K)
          Y2=XFLUX(IC,K)
          Y3=XFLUX(NJP1(IC),K)
c          BVEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/
c     *    (X2-X3)/((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2))
          BVEL=(X3**2*(Y1-Y2)+X1**2*(Y2-Y3))
     */(X1*X3*(X3-X1))
          AVEL=(Y1-Y2-BEL*X1)/(X1**2)
c          CVEL=Y1-AEL*X1**2-BEL*X1
          CVEL=Y2
c          IF (AVEL.EQ.0..AND.BVEL.EQ.0..AND.CVEL.EQ.0.) CD=1 !CBR: Never used
c          F1=(YNODE(1,IC)-YR(IC))/H1(INUM)
c          F2=(YNODE(1,INUM)-YR(IC))/H1(INUM)
c          F3=(YNODE(2,INUM)-YR(IC))/H1(INUM)
c          F4=(YNODE(2,IC)-YR(IC))/H1(INUM)
c          CNUM1=AVEL*(F3**3-F2**3)/3+BVEL*(F3**2-F2**2)/2+CVEL*(F3-F2)
c          CNUM=AVEL*(F4**3-F1**3)/3+BVEL*(F4**2-F1**2)/2+CVEL*(F4-F1)
          XFLUX(INUM,K)=XFLUX(IC,K)*(1./3.)  !*CNUM1/CNUM
      ENDIF
      ENDIF
      ENDDO
c      ENDDO
      
      
c!$OMP DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K=1,KBM1
c      DO I=1,NNIVGFEALL
      DO concurrent (I=1:NNIVGFEALL,K=1:KBM1) 
     * local(INUM,IC,X1,X3,Y1,Y2,Y3,BVEL,AVEL,CVEL) 
      INUM=NIVGFEALL(I,1)
      IC =NIP1(NIVGFEALL(I,2))
      IF (U(NIM1(INUM),K).LT.0.) THEN
      IF (XFLUX(IC,K).EQ.0) THEN
          XFLUX(INUM,K)=0
      ELSE
          X1=(YR(NJM1(IC))-YR(IC))/H1(INUM)
c          X2=0.
          X3=(YR(NJP1(IC))-YR(IC))/H1(INUM)
          Y1=XFLUX(NJM1(NIP1(IC)),K)
          Y2=XFLUX(NIP1(IC),K)
          Y3=XFLUX(NJP1(NIP1(IC)),K)
c          BVEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/
c     *    (X2-X3)/((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2))
          BVEL=(X3**2*(Y1-Y2)+X1**2*(Y2-Y3))
     */(X1*X3*(X3-X1))
c          AVEL=(Y1-Y2-BEL*(X1-X2))/(X1**2-X2**2)
          AVEL=(Y1-Y2-BEL*X1)/(X1**2)
c          CVEL=Y1-AEL*X1**2-BEL*X1
          CVEL=Y2
c          IF (AVEL.EQ.0..AND.BVEL.EQ.0..AND.CVEL.EQ.0.) CD=1 !CBR: Never used
c          F1=(YNODE(1,IC)-YR(IC))/H1(INUM)
c          F2=(YNODE(1,INUM)-YR(IC))/H1(INUM)
c          F3=(YNODE(2,INUM)-YR(IC))/H1(INUM)
c          F4=(YNODE(2,IC)-YR(IC))/H1(INUM)
c          CNUM1=AVEL*(F3**3-F2**3)/3+BVEL*(F3**2-F2**2)/2+CVEL*(F3-F2)
c          CNUM=AVEL*(F4**3-F1**3)/3+BVEL*(F4**2-F1**2)/2+CVEL*(F4-F1)
          XFLUX(INUM,K)=XFLUX(IC,K)*(1./3.)  !*CNUM1/CNUM
      ENDIF
      ENDIF
      ENDDO
c      ENDDO   
c!$OMP END TARGET teams
*                                                                     *
***********************************************************************
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I = 1,N_CTRD
      DO concurrent (I=1:N_CTRD) local(K) local_init(KBM1)
      IF (DUM(I).GT.0.0) THEN
          DO K = 2,KBM1
              UF(I,K) = .25*(W(I,K)+W(NIM1(I),K))*(U(I,K)+U(I,K-1))
          ENDDO
          DO K = 1,KBM1
              UF(I,K) = DZR(K)*(UF(I,K)-UF(I,K+1))*ARU(I)
          ENDDO
      ENDIF
      ENDDO
c!$OMP END TARGET

     
c#ifdef OMP
c!$OMP PARALLEL DO PRIVATE(K)
c#endif
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO I = 1,N_CTRD
c      DO K = 1,KBM1    
      DO concurrent (I=1:N_CTRD,K=1:KBM1)
      IF (DUM(I).GT.0.0) THEN
          UF(I,K)=UF(I,K)
     *    +XFLUX(I,K)-XFLUX(NIM1(I),K)+YFLUX(NJP1(I),K)-YFLUX(I,K)
      ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET
c#ifdef OMP
c!$OMP END PARALLEL DO
c#endif
      
***********************************************************************
*                             VG BONDOURY                             *
c!$OMP TARGET teams DEFAULTMAP(present: allocatable)
c!$OMP DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K=1,KBM1
c      DO II=1,NNIVGFSALL
      DO concurrent (II=1:NNIVGFSALL,K=1:KBM1) local(I) 
          I=NIVGFSALL(II,1)
          IF (DUM(I).GT.0.0) THEN
          UF(I,K)=UF(I,K)-YFLUX(NJP1(I),K)+YFLUX(I,K)
          ENDIF
      ENDDO
c      ENDDO
      
c!$OMP DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K=1,KBM1
c      DO II=1,NNIVGFNALL
      DO concurrent (II=1:NNIVGFNALL,K=1:KBM1) local(I) 
          I=NIVGFNALL(II,1)
          IF (DUM(I).GT.0.0) THEN
          UF(I,K)=UF(I,K)-YFLUX(NJP1(I),K)+YFLUX(I,K)
          ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET teams
*                                                                     *
***********************************************************************
       
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO N = 1, NUMEBC
c      DO K = 1,KBM1
      DO concurrent (N=1:NUMEBC,K=1:KBM1) local(IE,IC) 
          IE = NETA(N)
          IC = NCON(N)
          IF (FSM(NIP1(IE)).EQ.0.0
     *    .AND. (IE==NIM1(IC) .OR. IE==NIP1(IC))) THEN
                  CURV4(IE,K) = CURV4(NIM1(IE),K)
         ELSEIF (FSM(NIM1(IE)).EQ.0.0
     *    .AND. (IE==NIM1(IC) .OR. IE==NIP1(IC))) THEN
                  CURV4(IE,K) = CURV4(NIP1(IE),K)
         ENDIF
       ENDDO
c       ENDDO
c!$OMP END TARGET
        

C*********FORWARD TIME STEP + CURVATURE AND CORIOLIS TERMS
C         + BAROCLINIC PRESSURE TERM***
      
c#ifdef OMP
c!$OMP PARALLEL DO PRIVATE(K)
c#endif
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO I = 1,N_CTRD
c      DO K = 1,KBM1    
      DO concurrent (I=1:N_CTRD,K=1:KBM1)
      IF (DUM(I).GT.0.0) THEN
      UF(I,K) = UF(I,K)
     * -(CURV4(I,K)*D(I)*H1(I)*H2(I)*(V(NJP1(I),K)+V(I,K))
     * +CURV4(NIM1(I),K)*D(NIM1(I))*H1(NIM1(I))*H2(NIM1(I))
     * *(V(NIM1JP1(I),K)+V(NIM1(I),K))
     * -0.5*U(I,K)*(D(I)+D(NIM1(I)))
     * *(H3(I)+H3(NIM1(I)))*(CURV4(I,K)+CURV4(NIM1(I),K)))+DRHOX(I,K)
      UF(I,K) = UF(I,K)
     * -0.25*(H2(I)+H2(NIM1(I)))*(D(I)+D(NIM1(I)))*U(I,K)
     * *(V(NJP1(I),K)*DVM(NJP1(I))+V(I,K)*DVM(I)
     * +V(NIM1JP1(I),K)*DVM(NIM1JP1(I))+V(NIM1(I),K)*DVM(NIM1(I)))
     * *(H3(I)/(H1(I)*H2(I))-H3(NIM1(I))/(H1(NIM1(I))*H2(NIM1(I))))
     * /(DVM(NJP1(I))+DVM(I)+DVM(NIM1JP1(I))+DVM(NIM1(I))+1.E-30)
#ifdef EXP_EL              
      UF(I,K) = UF(I,K)
     * +0.5*(H2(I)+H2(NIM1(I)))*GRAV*DU(I)*(EL(I)-EL(NIM1(I)))
#endif  
      ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET
c#ifdef OMP
c!$OMP END PARALLEL DO
c#endif
     
C******* STEP FORWARD IN TIME *****************************************
      
c#ifdef OMP
c!$OMP PARALLEL DO PRIVATE(DUMT,DUMT1,DVMT,DVMT1,HUP,ELUP,
c     *    HDOWN,ELDOWN,DTN,K)
c#endif
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I = 1,N_CTRD
      DO concurrent (I=1:N_CTRD) local(K,DUMT,DUMT1,DVMT,DVMT1,HUP,
     * ELUP,HDOWN,ELDOWN,DTN) local_init(KBM1,DTI)
      IF (DUM(I) .GT. 0.0) THEN
        
      DUMT = FSMADD(NIM1(I))
      DUMT1 = FSMADD(I)
      DVMT=FSMADD(I)*FSMADD(NJM1(I))*FSMADD(NIM1JM1(I))*FSMADD(NIM1(I))
      DVMT1=FSMADD(I)*FSMADD(NJP1(I))*FSMADD(NIM1JP1(I))*FSMADD(NIM1(I))
           
      HUP = (H(I)*FSMADD(I)+H(NIM1(I))*FSMADD(NIM1(I))
     * +H(NIM1JP1(I))*FSMADD(NIM1JP1(I))
     * +H(NJP1(I))*FSMADD(NJP1(I)))
     * /(FSMADD(I)+FSMADD(NIM1(I))+FSMADD(NIM1JP1(I))+FSMADD(NJP1(I))
     * +1.E-10)
     
      ELUP = (EL(I)*FSM(I)+EL(NIM1(I))*FSM(NIM1(I))
     * +EL(NIM1JP1(I))*FSM(NIM1JP1(I))+EL(NJP1(I))*FSM(NJP1(I)))
     * /(FSM(I)+FSM(NIM1(I))+FSM(NIM1JP1(I))+FSM(NJP1(I))+1.E-10)  
     
      HDOWN = (H(I)*FSMADD(I)+H(NIM1(I))*FSMADD(NIM1(I))
     * +H(NIM1JM1(I))*FSMADD(NIM1JM1(I))
     * +H(NJM1(I))*FSMADD(NJM1(I)))
     * /(FSMADD(I)+FSMADD(NIM1(I))+FSMADD(NIM1JM1(I))+FSMADD(NJM1(I))
     * +1.E-10)
     
      ELDOWN = (EL(I)*FSM(I)+EL(NIM1(I))*FSM(NIM1(I))
     * +EL(NIM1JM1(I))*FSM(NIM1JM1(I))+EL(NJM1(I))*FSM(NJM1(I)))
     * /(FSM(I)+FSM(NIM1(I))+FSM(NIM1JM1(I))+FSM(NJM1(I))+1.E-10)
     
      DTN = 
     * ((HU(I)+2*AMAX1(H(NIM1(I)), -1.*EL(NIM1(I)))*DUMT)/(1.+2.*DUMT)
     * +(HU(I)+2*AMAX1(H(I), -1.*EL(I))*DUMT1)/(1.+2.*DUMT1) 
     * +(HU(I)+2.*AMAX1(HDOWN, -1.*ELDOWN)*DVMT)/(1.+2.*DVMT)
     * +(HU(I)+2.*AMAX1(HUP, -1.*ELUP)*DVMT1)/(1.+2.*DVMT1))/4.
      
      DTN = DTN+(EL(I)*FSM(I)+EL(NIM1(I))*FSM(NIM1(I)))
     * /(FSM(I)+FSM(NIM1(I))+1.E-10)
      DO K = 1,KBM1             
          UF(I,K) = (DTN*ARU(I)*U(I,K)-DTI*UF(I,K))/(DTN*ARU(I))   
      ENDDO  
      ENDIF
      ENDDO
c!$OMP END TARGET
c#ifdef OMP
c!$OMP END PARALLEL DO 
c#endif


      RETURN
        
      END SUBROUTINE ADVU_TVD_3RD
