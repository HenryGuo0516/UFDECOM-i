#include "DEFS.h"
      SUBROUTINE ADVV_TVD_3RD
      USE MOD_GLOBAL
c      DIMENSION AAMX(N_CTRDP1,KB),AAMY(N_CTRDP1,KB) !CBR: Never used.
      INTEGER IC,INUM 
c      INTEGER BD(3)    
      REAL X1,X2,X3,Y1,Y2,Y3,AVEL,BVEL,CVEL
      REAL F1,F2,F3,F4,CNUM1,CNUM
      REAL DVMT,DVMT1,DUMT,DUMT1,HLEFT,ELLEFT,HRIGHT,ELRIGHT,DTN
	REAL CFLX,CFLY,CFL,ISIGN,ISIGN1,CFLX1,RAT
      
c      VF=0.
c      XFLUX=0.
c      YFLUX=0.
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO I = 1, N_CTRDP1
c      DO K = 1, KB
      DO concurrent (I=1:N_CTRDP1,K=1:KB)
          VF(I,K)=0.
          XFLUX(I,K)=0.
          YFLUX(I,K)=0.
      ENDDO
c      ENDDO
c!$OMP END TARGET
         

c#ifdef OMP        
c!$OMP PARALLEL DO PRIVATE(CFLX,FUP,FLW,FS,CFLX1,FUP1,FLW1,FS1
c     * ,R,RATIO,K)
c#endif
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO I = 1,N_CTRD 
c      DO K = 1,KBM1
      DO concurrent (I=1:N_CTRD,K=1:KBM1) 
     * local(CFLX,FUP,CFLX1,FUP1,FLW1,FLW,FS,FS1,R,RATIO,
     *       CFLY,CFLY1) local_init(DTI)
c          PROD(I,K)=
c     * 4.*(D(I)*FSM(I)+D(NIM1(I))*FSM(NIM1(I))+D(NJM1(I))*FSM(NJM1(I))
c     * +D(NIM1JM1(I))*FSM(NIM1JM1(I)))
c     * /(FSM(I)+FSM(NIM1(I))+FSM(NJM1(I))+FSM(NIM1JM1(I))+1.E-30)
c     * *(AAM(I,K)*FSM(I)+AAM(NIM1(I),K)*FSM(NIM1(I))
c     * +AAM(NJM1(I),K)*FSM(NJM1(I))
c     * +AAM(NIM1JM1(I),K)*FSM(NIM1JM1(I)))
c     * /(FSM(I)+FSM(NIM1(I))+FSM(NJM1(I))+FSM(NIM1JM1(I))+1.E-30)
          
          CFLX = 0.5*(U(I,K)+U(NJM1(I),K))*DTI
     *    /(0.25*(H1(I)+H1(NIM1(I))+H1(NJM1(I))+H1(NIM1JM1(I))))
          IF(CFLX.GE.0.) THEN
              FUP = V(NIM1(I),K)
              CFLX1 = 0.5*(U(NIM1(I),K)+U(NIM1JM1(I),K))*DTI/(0.25*
     *        (H1(NIM1(I))+H1(NIM1JM1(I))+H1(NIM2(I))+H1(NIM2JM1(I))))
              IF(CFLX1.GE.0.) THEN
                  FUP1 = V(NIM2(I),K)
              ELSE
                  FUP1 = V(NIM1(I),K)
              ENDIF
              FLW1=0.5*((1+CFLX1)*V(NIM2(I),K)+(1.-CFLX1)*V(NIM1(I),K))
          ELSE
              FUP = V(I,K)
              CFLX1 = 0.5*(U(NIP1(I),K)+U(NIP1JM1(I),K))*DTI
     *        /(0.25*(H1(I)+H1(NIP1(I))+H1(NJM1(I))+H1(NIP1JM1(I))))
              IF(CFLX1.GE.0.) THEN
                  FUP1 = V(I,K)
              ELSE
                  FUP1 = V(NIP1(I),K)
              ENDIF
              FLW1 = 0.5*((1+CFLX1)*V(I,K)+(1.-CFLX1)*V(NIP1(I),K))
          ENDIF
          FLW = 0.5*((1+CFLX)*V(NIM1(I),K)+(1.-CFLX)*V(I,K)) 
          FS = FLW-FUP
          FS1 = FLW1-FUP1
          IF(ABS(FS).LE.1.E-10) THEN
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
          XFLUX(I,K)=(FUP+RATIO*FS)*0.5*(XMFLUX(I,K)+XMFLUX(NJM1(I),K))
          
          CFLY = 0.5*(V(NJP1(I),K)+V(I,K))*DTI/H2(I)
          IF(CFLY.GE.0.) THEN
              FUP = YMFLUX(I,K)*V(I,K)
              CFLY1 = 0.5*(V(I,K)+V(NJM1(I),K))*DTI/H2(NJM1(I))
              IF(CFLY1.GE.0.) THEN
                  FUP1 = YMFLUX(NJM1(I),K)*V(NJM1(I),K)
              ELSE
                  FUP1 = YMFLUX(I,K)*V(I,K)
              ENDIF
              FLW1=0.5*((1+CFLY1)*YMFLUX(NJM1(I),K)*V(NJM1(I),K)
     *        +(1.-CFLY1)*YMFLUX(I,K)*V(I,K))
          ELSE 
              FUP = YMFLUX(NJP1(I),K)*V(NJP1(I),K)
              CFLY1 = 0.5*(V(NJP1(I),K)+V(NJP2(I),K))*DTI/H2(NJP1(I))
              IF(CFLY1.GE.0.) THEN
                  FUP1 = YMFLUX(NJP1(I),K)*V(NJP1(I),K)
              ELSE
                  FUP1 = YMFLUX(NJP2(I),K)*V(NJP2(I),K)
              ENDIF
              FLW1 = 0.5*((1+CFLY1)*YMFLUX(NJP1(I),K)*V(NJP1(I),K)
     *        +(1.-CFLY1)*YMFLUX(NJP2(I),K)*V(NJP2(I),K))
          ENDIF
          FLW=0.5*((1+CFLY)*YMFLUX(I,K)*V(I,K)
     *    +(1.-CFLY)*YMFLUX(NJP1(I),K)*V(NJP1(I),K))
          FS = FLW-FUP
          FS1 = FLW1-FUP1
          IF(ABS(FS).LE.1.E-10) THEN
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
          YFLUX(I,K) = FUP+RATIO*FS
          
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
          XFLUX(IA,K)=.25*(XMFLUX(IA,K)+XMFLUX(NJM1(IA),K))*
     *    (V(IA,K)+V(NIM1(IA),K))
          YFLUX(IA,K)=.25*(YMFLUX(NJP1(IA),K)+YMFLUX(IA,K))*
     *    (V(NJP1(IA),K)+V(IA,K))
          XFLUX(IC,K)=.25*(XMFLUX(IC,K)+XMFLUX(NJM1(IC),K))*
     *    (V(IC,K)+V(NIM1(IC),K))
          YFLUX(IC,K)=.25*(YMFLUX(NJP1(IC),K)+YMFLUX(IC,K))*
     *    (V(NJP1(IC),K)+V(IC,K))
          XFLUX(IE,K)=.25*(XMFLUX(IE,K)+XMFLUX(NJM1(IE),K))*
     *    (V(IE,K)+V(NIM1(IE),K))
          YFLUX(IE,K)=.25*(YMFLUX(NJP1(IE),K)+YMFLUX(IE,K))*
     *    (V(NJP1(IE),K)+V(IE,K))
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
              XFLUX(IA,K)=.25*(XMFLUX(IA,K)+XMFLUX(NJM1(IA),K))*
     *        (V(IA,K)+V(NIM1(IA),K))
              YFLUX(IA,K)=.25*(YMFLUX(NJP1(IA),K)+YMFLUX(IA,K))*
     *        (V(NJP1(IA),K)+V(IA,K))
              XFLUX(IC,K)=.25*(XMFLUX(IC,K)+XMFLUX(NJM1(IC),K))*
     *        (V(IC,K)+V(NIM1(IC),K))
              YFLUX(IC,K)=.25*(YMFLUX(NJP1(IC),K)+YMFLUX(IC,K))*
     *        (V(NJP1(IC),K)+V(IC,K))
              XFLUX(IE,K)=.25*(XMFLUX(IE,K)+XMFLUX(NJM1(IE),K))*
     *        (V(IE,K)+V(NIM1(IE),K))
              YFLUX(IE,K)=.25*(YMFLUX(NJP1(IE),K)+YMFLUX(IE,K))*
     *        (V(NJP1(IE),K)+V(IE,K))
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
          PROD0=
     * 4.*(D(I)*FSM(I)+D(NIM1(I))*FSM(NIM1(I))+D(NJM1(I))*FSM(NJM1(I))
     * +D(NIM1JM1(I))*FSM(NIM1JM1(I)))
     * /(FSM(I)+FSM(NIM1(I))+FSM(NJM1(I))+FSM(NIM1JM1(I))+1.E-30)
     * *(AAM(I,K)*FSM(I)+AAM(NIM1(I),K)*FSM(NIM1(I))
     * +AAM(NJM1(I),K)*FSM(NJM1(I))
     * +AAM(NIM1JM1(I),K)*FSM(NIM1JM1(I)))
     * /(FSM(I)+FSM(NIM1(I))+FSM(NJM1(I))+FSM(NIM1JM1(I))+1.E-30)
          
          XFLUX(I,K) = XFLUX(I,K)-PROD0*((U(I,K)-U(NJM1(I),K))
     *    /(H2(I)+H2(NIM1(I))+H2(NJM1(I))+H2(NIM1JM1(I)))
     *    +(V(I,K)-V(NIM1(I),K))
     *    /(H1(I)+H1(NIM1(I))+H1(NJM1(I))+H1(NIM1JM1(I))))
     *    *.25*(H2(I)+H2(NIM1(I))+H2(NJM1(I))+H2(NIM1JM1(I)))
     *    *DVM(I)*DVM(NIM1(I)) 
          YFLUX(I,K)=YFLUX(I,K)-
     *    D(I)*AAM(I,K)*FSM(I)*2.*(V(NJP1(I),K)-V(I,K))*H1(I)/H2(I)
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
c      DO I=1,NNIVGFSALL
      DO concurrent (I=1:NNIVGFSALL,K=1:KBM1) 
     * local(INUM,IC,X1,X3,Y1,Y2,Y3,BVEL,AVEL,CVEL) 
      INUM=NIVGFSALL(I,1)
      IC =NIVGFSALL(I,2)
      IF (V(NJP1(INUM),K).GT.0.) THEN
      IF (YFLUX(IC,K).EQ.0) THEN
          YFLUX(INUM,K)=0.
      ELSE
          X1=(XR(NIM1(IC))-XR(IC))/H2(INUM)
c          X2=0
          X3=(XR(NIP1(IC))-XR(IC))/H2(INUM)
          Y1=YFLUX(NIM1(IC),K)
          Y2=YFLUX(IC,K)
          Y3=YFLUX(NIP1(IC),K)
c          BVEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/
c     *    (X2-X3)/((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2))
          BVEL=(X3**2*(Y1-Y2)+X1**2*(Y2-Y3))
     */(X1*X3*(X3-X1))
c          AVEL=(Y1-Y2-BEL*(X1-X2))/(X1**2-X2**2)
          AVEL=(Y1-Y2-BEL*X1)/(X1**2)
c          CVEL=Y1-AEL*X1**2-BEL*X1
          CVEL=Y2
c          IF (AVEL.EQ.0..AND.BVEL.EQ.0..AND.CVEL.EQ.0.) CD=1 !CBR: Never used
c          F1=(XNODE(1,IC)-XR(IC))/H2(INUM)
c          F2=(XNODE(1,INUM)-XR(IC))/H2(INUM)
c          F3=(XNODE(4,INUM)-XR(IC))/H2(INUM)
c          F4=(XNODE(4,IC)-XR(IC))/H2(INUM)
c          CNUM1=AVEL*(F3**3-F2**3)/3+BVEL*(F3**2-F2**2)/2+CVEL*(F3-F2)
c          CNUM=AVEL*(F4**3-F1**3)/3+BVEL*(F4**2-F1**2)/2+CVEL*(F4-F1)
          YFLUX(INUM,K)=YFLUX(IC,K)*(1./3.)  !*CNUM1/CNUM
      ENDIF
      ENDIF
      ENDDO
c      ENDDO
      
c!$OMP DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K=1,KBM1
c      DO I=1,NNIVGFNALL
      DO concurrent (I=1:NNIVGFNALL,K=1:KBM1) 
     * local(INUM,IC,X1,X3,Y1,Y2,Y3,BVEL,AVEL,CVEL) 
      INUM=NIVGFNALL(I,1)
      IC =NIVGFNALL(I,2)
      IF (V(NJM1(INUM),K).LT.0) THEN
      IF (YFLUX(IC,K).EQ.0) THEN
          YFLUX(INUM,K)=0.
      ELSE
          X1=(XR(NIM1(IC))-XR(IC))/H2(INUM)
c          X2=0
          X3=(XR(NIP1(IC))-XR(IC))/H2(INUM)
          Y1=YFLUX(NIM1(IC),K)
          Y2=YFLUX(IC,K)
          Y3=YFLUX(NIP1(IC),K)
c          BVEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/
c     *    (X2-X3)/((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2))
          BVEL=(X3**2*(Y1-Y2)+X1**2*(Y2-Y3))
     */(X1*X3*(X3-X1))
          AVEL=(Y1-Y2-BEL*X1)/(X1**2)
c          CVEL=Y1-AEL*X1**2-BEL*X1
          CVEL=Y2
c          IF (AVEL.EQ.0..AND.BVEL.EQ.0..AND.CVEL.EQ.0.) CD=1 !CBR: Never used
c          F1=(XNODE(1,IC)-XR(IC))/H2(INUM)
c          F2=(XNODE(1,INUM)-XR(IC))/H2(INUM)
c          F3=(XNODE(4,INUM)-XR(IC))/H2(INUM)
c          F4=(XNODE(4,IC)-XR(IC))/H2(INUM)
c          CNUM1=AVEL*(F3**3-F2**3)/3+BVEL*(F3**2-F2**2)/2+CVEL*(F3-F2)
c          CNUM=AVEL*(F4**3-F1**3)/3+BVEL*(F4**2-F1**2)/2+CVEL*(F4-F1)
          YFLUX(INUM,K)=YFLUX(IC,K)*(1./3.)  !*CNUM1/CNUM
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
      IF(DVM(I).GT.0.0) THEN
          DO K = 2,KBM1
              VF(I,K) = .25*(W(I,K)+W(NJM1(I),K))*(V(I,K)+V(I,K-1))
          ENDDO
          DO K = 1,KBM1
              VF(I,K) = DZR(K)*(VF(I,K)-VF(I,K+1))*ARV(I)
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
      IF(DVM(I).GT.0.0) THEN
          VF(I,K)=VF(I,K)+
     *    XFLUX(NIP1(I),K)-XFLUX(I,K)+YFLUX(I,K)-YFLUX(NJM1(I),K)
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
c      DO II=1,NNIVGFWALL
      DO concurrent (II=1:NNIVGFWALL,K=1:KBM1) local(I) 
          I=NIVGFWALL(II,1)
          IF (DVM(I).GT.0.0) THEN
          VF(I,K)=VF(I,K)-XFLUX(NIP1(I),K)+XFLUX(I,K)
          ENDIF
      ENDDO
c      ENDDO
      
c!$OMP DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K=1,KBM1
c      DO II=1,NNIVGFEALL
      DO concurrent (II=1:NNIVGFEALL,K=1:KBM1) local(I) 
          I=NIVGFEALL(II,1)
          IF (DVM(I).GT.0.0) THEN
          VF(I,K)=VF(I,K)-XFLUX(NIP1(I),K)+XFLUX(I,K)
          ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET teams
*                                                                     *
***********************************************************************
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO N = 1,NUMEBC
c      DO K = 1,KBM1
      DO concurrent (N=1:NUMEBC,K=1:KBM1) local(IE,IC) 
      IE = NETA(N)
      IC = NCON(N)
      IF(FSM(NJP1(IE)).EQ.0.0.AND.(IE==NJP1(IC).OR.IE==NJM1(IC))) THEN
              CURV4(IE,K) = CURV4(NJM1(IE),K)
      ELSEIF(FSM(NJM1(IE)).EQ.0.0.AND.(IE==NJP1(IC).OR.IE==NJM1(IC)))
     * THEN
              CURV4(IE,K) = CURV4(NJP1(IE),K)
      ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET

      
c#ifdef OMP
c!$OMP PARALLEL DO PRIVATE(K)
c#endif
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO I = 1,N_CTRD
c      DO K = 1,KBM1    
      DO concurrent (I=1:N_CTRD,K=1:KBM1)
      IF(DVM(I).GT.0.0) THEN
      VF(I,K) = VF(I,K)
     * +(CURV4(I,K)*D(I)*H1(I)*H2(I)*(U(NIP1(I),K)+U(I,K))
     * +CURV4(NJM1(I),K)*D(NJM1(I))*H1(NJM1(I))*H2(NJM1(I))
     * *(U(NIP1JM1(I),K)+U(NJM1(I),K))
     * -0.5*(H3(I)+H3(NJM1(I)))*(D(I)+D(NJM1(I)))*V(I,K)
     * *(CURV4(I,K)+CURV4(NJM1(I),K)))+DRHOY(I,K) 
      VF(I,K) = VF(I,K)
     * -0.25*(H1(I)+H1(NJM1(I)))*(D(I)+D(NJM1(I)))*V(I,K)
     * *(U(NIP1(I),K)*DUM(NIP1(I))+U(I,K)*DUM(I)
     * +U(NIP1JM1(I),K)*DUM(NIP1JM1(I))+U(NJM1(I),K)*DUM(NJM1(I)))
     * *(H3(I)/(H1(I)*H2(I))-H3(NJM1(I))/(H1(NJM1(I))*H2(NJM1(I))))
     * /(DUM(NIP1(I))+DUM(I)+DUM(NIP1JM1(I))+DUM(NJM1(I))+1.E-30)
#ifdef EXP_EL
      VF(I,K)=VF(I,K)
     *+0.5*(H1(I)+H1(NJM1(I)))*GRAV*DV(I)*(EL(I)-EL(NJM1(I)))
#endif           
      ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET
c#ifdef OMP
c!$OMP END PARALLEL DO
c#endif
      

      
c#ifdef OMP         
c!$OMP PARALLEL DO PRIVATE(DVMT,DVMT1,DUMT,DUMT1,HLEFT,ELLEFT
c     * ,HRIGHT,ELRIGHT,DTN,K)
c#endif
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I = 1,N_CTRD
      DO concurrent (I=1:N_CTRD) local(K,DUMT,DUMT1,DVMT,DVMT1,HLEFT,
     * ELLEFT,HRIGHT,ELRIGHT,DTN) local_init(KBM1,DTI)
      IF(DVM(I).GT.0.0) THEN
      DVMT = FSMADD(NJM1(I))
      DVMT1 = FSMADD(I)
      DUMT=FSMADD(I)*FSMADD(NJM1(I))*FSMADD(NIM1JM1(I))*FSMADD(NIM1(I))
      DUMT1=FSMADD(I)*FSMADD(NIP1(I))*FSMADD(NIP1JM1(I))*FSMADD(NJM1(I))
          
      HLEFT = (H(I)*FSMADD(I)+H(NIM1(I))*FSMADD(NIM1(I))
     * +H(NIM1JM1(I))*FSMADD(NIM1JM1(I))
     * +H(NJM1(I))*FSMADD(NJM1(I)))
     * /(FSMADD(I)+FSMADD(NIM1(I))+FSMADD(NIM1JM1(I))+FSMADD(NJM1(I))
     * +1.E-10)
      
      ELLEFT = (EL(I)*FSM(I)+EL(NIM1(I))*FSM(NIM1(I))
     * +EL(NIM1JM1(I))*FSM(NIM1JM1(I))+EL(NJM1(I))*FSM(NJM1(I)))
     * /(FSM(I)+FSM(NIM1(I))+FSM(NIM1JM1(I))+FSM(NJM1(I))+1.E-10)
      
      HRIGHT = (H(I)*FSMADD(I)+H(NIP1(I))*FSMADD(NIP1(I))
     * +H(NIP1JM1(I))*FSMADD(NIP1JM1(I))
     * +H(NJM1(I))*FSMADD(NJM1(I)))
     * /(FSMADD(I)+FSMADD(NIP1(I))+FSMADD(NIP1JM1(I))+FSMADD(NJM1(I))
     * +1.E-10)
      
      ELRIGHT = (EL(I)*FSM(I)+EL(NIP1(I))*FSM(NIP1(I))
     * +EL(NIP1JM1(I))*FSM(NIP1JM1(I))+EL(NJM1(I))*FSM(NJM1(I)))
     * /(FSM(I)+FSM(NIP1(I))+FSM(NIP1JM1(I))+FSM(NJM1(I))+1.E-10)
      
      DTN = 
     * ((HV(I)+2*AMAX1(H(NJM1(I)),-1.*EL(NJM1(I)))*DVMT)/(1.+2.*DVMT)
     * +(HV(I)+2*AMAX1(H(I),-1.*EL(I))*DVMT1)/(1.+2.*DVMT1)
     * +(HV(I)+2.*AMAX1(HLEFT,-1.*ELLEFT)*DUMT)/(1.+2.*DUMT)
     * +(HV(I)+2.*AMAX1(HRIGHT,-1.*ELRIGHT)*DUMT1)/(1.+2.*DUMT1))/4.
      
      DTN = DTN+(EL(I)*FSM(I)+EL(NJM1(I))*FSM(NJM1(I)))
     * /(FSM(I)+FSM(NJM1(I))+1.E-10)
      DO K = 1,KBM1
      VF(I,K) = (DTN*ARV(I)*V(I,K)-DTI*VF(I,K))/(DTN*ARV(I))
      ENDDO       
      ENDIF
           
      ENDDO
c!$OMP END TARGET
c#ifdef OMP
c!$OMP END PARALLEL DO 
c#endif



      RETURN
      END SUBROUTINE ADVV_TVD_3RD
