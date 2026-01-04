#include "DEFS.h"
         
      SUBROUTINE ELTEST
      
C     Rewriting by MaRin for special eletion solve way.
C     Any problems can connet 913604577@qq.com  
      
      USE MOD_GLOBAL
      USE MOD_WEIR
	use mod_solver ! == pjy
      USE TGSpeed !CBR

c      INTEGER(KIND=8) ::COUNT1,COUNT2,COUNT_RATE 
      
      integer iter1, iter2   ! == pjy 
      INTEGER SITENUM,NWPOI,K,KK,IC,ID
      INTEGER IPARM(12)
      REAL X1,X2,X3,Y1,Y2,Y3
      REAL F1,F2,F3,F4,DL
      REAL DIAGONALUC1,VC1
      REAL AVEL,BVEL,CVEL
      REAL AEL,BEL,CEL
      REAL A1,A2,A3,A4,A5
      REAL CNUM1,CNUM2,CNUM,DNUM
c      integer NCNUM !CBR
      REAL(KIND=8) RPARM(12)
!      INTEGER,ALLOCATABLE :: IWKSP(:) !never used
!      REAL(KIND=8),ALLOCATABLE :: WKSP(:) !never used
!      REAL(KIND=4),ALLOCATABLE :: UPOI(:) !never used
      CHARACTER*80 FNAME,FN1
!      DATA  LEVEL /0/,IDGTS /2/
!      AAUU=0.
!      AAVV=0. 
!      XMFLUX=0.
!      YMFLUX=0.
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 1, KB
c      DO I = 1, N_CTRDP1
      DO concurrent (I=1:N_CTRDP1,K=1:KB)
          XMFLUX(I,K) = 0.0
          YMFLUX(I,K) = 0.0
      ENDDO
c      ENDDO
c!$OMP END TARGET      
#if defined LBCel_gra0 || defined LBCel_cha || defined LBCel_rad
      CALL BCOND(10)
#else
      CALL BCOND(3)
#endif
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 1,KBM1
c      DO I = 1,N_CTRD
      DO concurrent (I=1:N_CTRD,K=1:KBM1) local(UC1,VC1)
      IF (DUM(I).NE.0.0) THEN
          UC1 = 0.5*UF(I,K)*(H2(I)+H2(NIM1(I)))
     *    -0.25*(H3(I)+H3(NIM1(I)))/(H1(I)+H1(NIM1(I)))
     *    *(V(I,K)+V(NIM1(I),K)+V(NJP1(I),K)+V(NIM1JP1(I),K)) 
          XMFLUX(I,K) = UC1*DU(I)
      ENDIF
      IF (DVM(I).NE.0.0) THEN
          VC1 = 0.5*VF(I,K)*(H1(I)+H1(NJM1(I)))
     *    -0.25*(H3(I)+H3(NJM1(I)))/(H2(I)+H2(NJM1(I)))
     *    *(U(I,K)+U(NIP1(I),K)+U(NJM1(I),K)+U(NIP1JM1(I),K))
          YMFLUX(I,K) = VC1*DV(I)
      ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET  
      
c -------------- Record Speed:
      CALL SYSTEM_CLOCK(COUNT3,COUNT_RATE)
      TG31=TG31+(COUNT3-COUNT2)/REAL(COUNT_RATE)

c!$OMP TARGET teams DEFAULTMAP (present: allocatable)
c!$OMP DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO II=1,NNIAGCW
c      DO K=1,KBM1
      DO concurrent (II=1:NNIAGCW,K=1:KBM1) local(I)
          I=NIAGCW(II)
          IF (DUM(NIP1(I)).GT.0.0) THEN 
          XMFLUX(NIP1(I),K)=XMFLUX(AIJ(I,6),K)+XMFLUX(AIJ(I,7),K)
     *+XMFLUX(AIJ(I,8),K)
          ENDIF
      ENDDO
c      ENDDO
c!$OMP DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO II=1,NNIAGCE
c          DO K=1,KBM1
      DO concurrent (II=1:NNIAGCE,K=1:KBM1) local(I)
          I=NIAGCE(II)
          IF (DUM(I).GT.0.0) THEN 
          XMFLUX(I,K)=XMFLUX(AIJ(I,9),K)+XMFLUX(AIJ(I,10),K)
     *+XMFLUX(AIJ(I,11),K)
          ENDIF
      ENDDO
c      ENDDO
c!$OMP DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO II=1,NNIAGCN
c      DO K=1,KBM1
      DO concurrent (II=1:NNIAGCN,K=1:KBM1) local(I)
          I=NIAGCN(II)
          IF (DVM(I).GT.0.0) THEN 
          YMFLUX(I,K)=YMFLUX(AIJ(I,9),K)+YMFLUX(AIJ(I,10),K)
     *+YMFLUX(AIJ(I,11),K)
          ENDIF
      ENDDO
c      ENDDO
c!$OMP DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO II=1,NNIAGCS
c      DO K=1,KBM1
      DO concurrent (II=1:NNIAGCS,K=1:KBM1) local(I)
          I=NIAGCS(II)
          IF (DVM(NJP1(I)).GT.0.0) THEN 
          YMFLUX(NJP1(I),K)=YMFLUX(AIJ(I,6),K)+YMFLUX(AIJ(I,7),K)
     *+YMFLUX(AIJ(I,8),K)
          ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET teams
      

c -------------- Record Speed:
      CALL SYSTEM_CLOCK(COUNT4,COUNT_RATE)
      TG32=TG32+(COUNT4-COUNT3)/REAL(COUNT_RATE)
      
      HGRAVDTI2=0.5*GRAV*DTI**2
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
!      DO K = 1,KBM1
c      DO I = 1,N_CTRD
      DO concurrent (I=1:N_CTRD) local_init(HGRAVDTI2)
      IF (DUM(I).NE.0.0) THEN
          AAUU(I)=HGRAVDTI2*DU(I)*(H2(I)+H2(NIM1(I)))**2
     *    /(DJ(I)+DJ(NIM1(I)))
      else
          AAUU(I)=0.
      ENDIF
      IF (DVM(I).NE.0.0) THEN
          AAVV(I)=HGRAVDTI2*DV(I)*(H1(I)+H1(NJM1(I)))**2
     *    /(DJ(I)+DJ(NJM1(I)))
      else
          AAVV(I)=0.
      ENDIF
      ENDDO
!      ENDDO
c!$OMP END TARGET  

      CALL BCOND(4)
      
      
#ifdef WEIR
      CALL ADDWEIR(2)
      CALL ADDWEIR(3)
#endif
************************************************************************
************************************************************************
      NPOI=N_CTRD_AG+N_CTRD_EVG+1
      NWPOI=6*NPOI+4*ITMAX
!      ALLOCATE(WKSP(NWPOI))       ;WKSP=0.0 !never used
!      ALLOCATE(UPOI(NPOI))        ;UPOI=0.0 !never used
!      ALLOCATE(IWKSP(3*NPOI))     ;IWKSP=0 !never used
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I=1,ANUM
      DO concurrent (I=1:ANUM)
          APOI(I)=0.
      ENDDO
c!$OMP END TARGET  
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I=1,NPOI
      DO concurrent (I=1:NPOI)
          RHSPOI(I)=0.
      ENDDO
c!$OMP END TARGET  
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I=1,NPOI
      DO concurrent (I=1:NPOI) local(JNUM,K) local_init(KBM1,DTI)
      IF (AIJ(I,1).EQ.5) THEN
          JNUM=IAPOI(I)-1
          APOI(JNUM+AIJ(I,2))=-AAVV(I)
          APOI(JNUM+AIJ(I,3))=-AAUU(I)
          APOI(JNUM+AIJ(I,4))=DJ(I)+AAVV(I)+AAUU(I)
     *+AAUU(NIP1(I))+AAVV(NJP1(I))
          APOI(JNUM+AIJ(I,5))=-AAUU(NIP1(I))
          APOI(JNUM+AIJ(I,6))=-AAVV(NJP1(I))
c          RHSPOI(I)=0.0
          DO K=1,KBM1
              RHSPOI(I)=RHSPOI(I)+DZ(K)*
     *    (XMFLUX(NIP1(I),K)-XMFLUX(I,K)+
     *    YMFLUX(NJP1(I),K)-YMFLUX(I,K))
          ENDDO
          RHSPOI(I)=DJ(I)*EL(I)-RHSPOI(I)*DTI
      ELSEIF(AIJ(I,1).EQ.1) THEN
c          JNUM=IAPOI(I)-1
c          APOI(JNUM+1)=99999.
          APOI(IAPOI(I))=99999.
          RHSPOI(I)=0.
      ENDIF
      ENDDO
c!$OMP END TARGET  


      
#ifdef TIDE_EL      
      DO I=1,NUMEBC
          IC = NETA(I)
          JNUM=IAPOI(IC)-1
          APOI(JNUM+AIJ(IC,2))=0.
          APOI(JNUM+AIJ(IC,3))=0.
          APOI(JNUM+AIJ(IC,5))=0.
          APOI(JNUM+AIJ(IC,6))=0.
          APOI(JNUM+AIJ(IC,4))=99999.
          RHSPOI(IC)=0!99999.*ELF(IC)
      ENDDO
#endif
***********************************************************************
*                                                                     *
*                             AG BOUNDARY                             *
*                                                                     *
***********************************************************************       
c!$OMP TARGET teams DEFAULTMAP(present: allocatable)
c!$OMP DISTRIBUTE PARALLEL DO
c      DO 11 II=1,NNIAGCW
      DO concurrent (II=1:NNIAGCW) local(I,JNUM,K) local_init(KBM1,DTI)
      I=NIAGCW(II)
!      IF (AIJ(I,1).EQ.1) GO TO 11
      IF (AIJ(I,1)/=1) THEN
      JNUM=IAPOI(I)-1
      APOI(JNUM+1)=-AAVV(I)
      APOI(JNUM+2)=-AAUU(I)
      APOI(JNUM+3)=DJ(I)+AAVV(I)+AAUU(I)+AAVV(NJP1(I))
      APOI(JNUM+4)=-AAVV(NJP1(I))
      APOI(JNUM+5)=-AAUU(AIJ(I,6))*DUM(NIP1(I))
      APOI(JNUM+6)=-AAUU(AIJ(I,7))*DUM(NIP1(I))
      APOI(JNUM+7)=-AAUU(AIJ(I,8))*DUM(NIP1(I))
      APOI(JNUM+8)=AAUU(AIJ(I,6))*DUM(NIP1(I))
      APOI(JNUM+9)=AAUU(AIJ(I,7))*DUM(NIP1(I))
      APOI(JNUM+10)=AAUU(AIJ(I,8))*DUM(NIP1(I))
      DO K=1,KBM1
          RHSPOI(I)=RHSPOI(I)+DZ(K)*
     *    (XMFLUX(NIP1(I),K)-XMFLUX(I,K)+
     *    YMFLUX(NJP1(I),K)-YMFLUX(I,K))
      ENDDO
      RHSPOI(I)=DJ(I)*EL(I)-RHSPOI(I)*DTI
      ENDIF
11    ENDDO
c!$OMP DISTRIBUTE PARALLEL DO
c      DO 12 II=1,NNIAGCE
      DO concurrent (II=1:NNIAGCE) local(I,JNUM,K) local_init(KBM1,DTI)
      I=NIAGCE(II)
!      IF (AIJ(I,1).EQ.1) GO TO 12
      IF (AIJ(I,1)/=1) THEN
      JNUM=IAPOI(I)-1
      APOI(JNUM+1)=-AAVV(I)
      APOI(JNUM+2)=DJ(I)+AAVV(I)+AAUU(NIP1(I))+AAVV(NJP1(I))
      APOI(JNUM+3)=-AAUU(NIP1(I))
      APOI(JNUM+4)=-AAVV(NJP1(I))
      APOI(JNUM+5)=-AAUU(AIJ(I,9))*DUM(I)
      APOI(JNUM+6)=-AAUU(AIJ(I,10))*DUM(I)
      APOI(JNUM+7)=-AAUU(AIJ(I,11))*DUM(I)
      APOI(JNUM+8)=AAUU(AIJ(I,9))*DUM(I)
      APOI(JNUM+9)=AAUU(AIJ(I,10))*DUM(I)
      APOI(JNUM+10)=AAUU(AIJ(I,11))*DUM(I)
      DO K=1,KBM1
          RHSPOI(I)=RHSPOI(I)+DZ(K)*
     *    (XMFLUX(NIP1(I),K)-XMFLUX(I,K)+
     *    YMFLUX(NJP1(I),K)-YMFLUX(I,K))
      ENDDO
      RHSPOI(I)=DJ(I)*EL(I)-RHSPOI(I)*DTI
      ENDIF
12    ENDDO
c!$OMP DISTRIBUTE PARALLEL DO
c      DO 13 II=1,NNIAGCN
      DO concurrent (II=1:NNIAGCN) local(I,JNUM,K) local_init(KBM1,DTI)
      I=NIAGCN(II)
!      IF (AIJ(I,1).EQ.1) GO TO 13
      IF (AIJ(I,1)/=1) THEN
      JNUM=IAPOI(I)-1
      APOI(JNUM+1)=-AAUU(I)
      APOI(JNUM+2)=DJ(I)+AAUU(I)+AAUU(NIP1(I))+AAVV(NJP1(I))
      APOI(JNUM+3)=-AAUU(NIP1(I))
      APOI(JNUM+4)=-AAVV(NJP1(I))
      APOI(JNUM+5)=-AAVV(AIJ(I,9))*DVM(I)
      APOI(JNUM+6)=-AAVV(AIJ(I,10))*DVM(I)
      APOI(JNUM+7)=-AAVV(AIJ(I,11))*DVM(I)
      APOI(JNUM+8)=AAVV(AIJ(I,9))*DVM(I)
      APOI(JNUM+9)=AAVV(AIJ(I,10))*DVM(I)
      APOI(JNUM+10)=AAVV(AIJ(I,11))*DVM(I)
      DO K=1,KBM1
          RHSPOI(I)=RHSPOI(I)+DZ(K)*
     *    (XMFLUX(NIP1(I),K)-XMFLUX(I,K)+
     *    YMFLUX(NJP1(I),K)-YMFLUX(I,K))
      ENDDO
      RHSPOI(I)=DJ(I)*EL(I)-RHSPOI(I)*DTI
      ENDIF
13    ENDDO
c!$OMP DISTRIBUTE PARALLEL DO
c      DO 14 II=1,NNIAGCS
      DO concurrent (II=1:NNIAGCS) local(I,JNUM,K) local_init(KBM1,DTI)
      I=NIAGCS(II)
!      IF (AIJ(I,1).EQ.1) GO TO 14
      IF (AIJ(I,1)/=1) THEN
      JNUM=IAPOI(I)-1
      APOI(JNUM+1)=-AAVV(I)
      APOI(JNUM+2)=-AAUU(I)
      APOI(JNUM+3)=DJ(I)+AAUU(I)+AAUU(NIP1(I))+AAVV(I)
      APOI(JNUM+4)=-AAUU(NIP1(I))
      APOI(JNUM+5)=-AAVV(AIJ(I,6))*DVM(NJP1(I))
      APOI(JNUM+6)=-AAVV(AIJ(I,7))*DVM(NJP1(I))
      APOI(JNUM+7)=-AAVV(AIJ(I,8))*DVM(NJP1(I))
      APOI(JNUM+8)=AAVV(AIJ(I,6))*DVM(NJP1(I))
      APOI(JNUM+9)=AAVV(AIJ(I,7))*DVM(NJP1(I))
      APOI(JNUM+10)=AAVV(AIJ(I,8))*DVM(NJP1(I))
      DO K=1,KBM1
          RHSPOI(I)=RHSPOI(I)+DZ(K)*
     *    (XMFLUX(NIP1(I),K)-XMFLUX(I,K)+
     *    YMFLUX(NJP1(I),K)-YMFLUX(I,K))
      ENDDO
      RHSPOI(I)=DJ(I)*EL(I)-RHSPOI(I)*DTI
      ENDIF
14    ENDDO
c!$OMP END TARGET teams
      
c -------------- Record Speed:
      CALL SYSTEM_CLOCK(COUNT3,COUNT_RATE)
      TG33=TG33+(COUNT3-COUNT4)/REAL(COUNT_RATE)
      
      
***********************************************************************
*                                                                     *
*                             VG BOUNDARY                             *
*                                                                     *
*********************************************************************** 
      AVEL=0.
      BVEL=0.
      CVEL=1.
      
      
c!$OMP TARGET teams DEFAULTMAP(present: allocatable)
c!$OMP DISTRIBUTE PARALLEL DO
c      DO 21 II=1,NNIVGFW
      DO concurrent (II=1:NNIVGFW) local(I,JNUM,IC,DL,X1,X3,Y1,Y2,Y3,
     * AEL,BEL,CEL,F1,F2,F3,F4,A1,A2,A3,A4,A5,CNUM1,CNUM2,CNUM,K)
     * local_init(KBM1,DTI)
      I=NIVGFW(II)
!      IF (AIJ(I,1).EQ.1) GO TO 21
      if (AIJ(I,1)/=1) then
      JNUM=IAPOI(I)-1
      IC=AIJ(I,3)
      DL=(YNODE(1,IC)+YNODE(2,IC))/2
      X1=((YNODE(1,NJM1(IC))+YNODE(2,NJM1(IC)))/2-DL)/H1(I)
!      X2=0.0
      X3=((YNODE(1,NJP1(IC))+YNODE(2,NJP1(IC)))/2-DL)/H1(I)
     !! Y1=U(NJM1(IC))
     !! Y2=U(IC)
     !! Y3=U(NJP1(IC))
     !! BVEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/
     !!*(X2-X3)/((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2))
     !! AVEL=(Y1-Y2-BVEL*(X1-X2))/(X1**2-X2**2)
     !! CVEL=Y1-AVEL*X1**2-BVEL*X1
     !! IF (AVEL.EQ.0..AND.BVEL.EQ.0..AND.CVEL.EQ.0.) CVEL=1
      Y1=DU(NJM1(IC))*DUM(NJM1(IC))+DU(IC)*(1-DUM(NJM1(IC)))
      Y2=DU(IC)*DUM(IC)
      Y3=DU(NJP1(IC))*DUM(NJP1(IC))+DU(IC)*(1-DUM(NJP1(IC)))
!      BEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/(X2-X3)/ 
!     *((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2)) 
!      AEL=(Y1-Y2-BEL*(X1-X2))/(X1**2-X2**2)
!      CEL=Y1-AEL*X1**2-BEL*X1
      BEL=(X3**2*(Y1-Y2)+X1**2*(Y2-Y3))/(X1*X3*(X3-X1)) 
      AEL=(Y1-Y2-BEL*X1)/(X1**2)
      CEL=Y2
      IF (AEL.EQ.0..AND.BEL.EQ.0..AND.CEL.EQ.0.) CEL=1
      F1=(YNODE(1,IC)-DL)/H1(I)
      F2=(YNODE(1,I)-DL)/H1(I)
      F3=(YNODE(2,I)-DL)/H1(I)
      F4=(YNODE(2,IC)-DL)/H1(I)
      A1=AVEL*AEL/5
      A2=(AVEL*BEL+AEL*BVEL)/4
      A3=(AVEL*CEL+AEL*CVEL+BVEL*BEL)/3
      A4=(BVEL*CEL+BEL*CVEL)/2
      A5=CVEL*CEL
      CNUM1=A1*(F1**5-F4**5)+A2*(F1**4-F4**4)+A3*(F1**3-F4**3)+
     *A4*(F1**2-F4**2)+A5*(F1-F4)
      CNUM2=A1*(F2**5-F3**5)+A2*(F2**4-F3**4)+A3*(F2**3-F3**3)+
     *A4*(F2**2-F3**2)+A5*(F2-F3)
      CNUM=CNUM2/CNUM1
      APOI(JNUM+1)=-CNUM*AAUU(IC)
      APOI(JNUM+2)=CNUM*AAUU(IC)
      APOI(JNUM+3)=-AAVV(I)
      APOI(JNUM+4)=DJ(I)+AAVV(I)+AAUU(NIP1(I))+AAVV(NJP1(I))
      APOI(JNUM+5)=-AAUU(NIP1(I))
      APOI(JNUM+6)=-AAVV(NJP1(I))
      DO K=1,KBM1
          RHSPOI(I)=RHSPOI(I)+DZ(K)*
     *    (XMFLUX(NIP1(I),K)-CNUM*XMFLUX(IC,K)+
     *    YMFLUX(NJP1(I),K)-YMFLUX(I,K))
      ENDDO
      RHSPOI(I)=DJ(I)*EL(I)-RHSPOI(I)*DTI
      endif
21    ENDDO
c!$OMP DISTRIBUTE PARALLEL DO
c      DO 22 II=1,NNIVGFE
      DO concurrent (II=1:NNIVGFE) local(I,JNUM,IC,DL,X1,X3,Y1,Y2,Y3,
     * AEL,BEL,CEL,F1,F2,F3,F4,A1,A2,A3,A4,A5,CNUM1,CNUM2,CNUM,K)
     * local_init(KBM1,DTI)
      I=NIVGFE(II)
!      IF (AIJ(I,1).EQ.1) GO TO 22
      if (AIJ(I,1)/=1) then
      JNUM=IAPOI(I)-1
      IC=AIJ(I,3)
      DL=(YNODE(1,IC)+YNODE(2,IC))/2
      X1=((YNODE(1,NJM1(IC))+YNODE(2,NJM1(IC)))/2-DL)/H1(I)
!      X2=0.0
      X3=((YNODE(1,NJP1(IC))+YNODE(2,NJP1(IC)))/2-DL)/H1(I)
     !! Y1=U(NJM1(IC))
     !! Y2=U(IC)
     !! Y3=U(NJP1(IC))
     !! BVEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/
     !!*(X2-X3)/((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2))
     !! AVEL=(Y1-Y2-BVEL*(X1-X2))/(X1**2-X2**2)
     !! CVEL=Y1-AVEL*X1**2-BVEL*X1
     !! IF (AVEL.EQ.0..AND.BVEL.EQ.0..AND.CVEL.EQ.0.) CVEL=1
      Y1=DU(NJM1(IC))*DUM(NJM1(IC))+DU(IC)*(1-DUM(NJM1(IC)))
      Y2=DU(IC)*DUM(IC)
      Y3=DU(NJP1(IC))*DUM(NJP1(IC))+DU(IC)*(1-DUM(NJP1(IC)))
!      BEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/(X2-X3)/ 
!     *((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2)) 
!      AEL=(Y1-Y2-BEL*(X1-X2))/(X1**2-X2**2)
!      CEL=Y1-AEL*X1**2-BEL*X1
      BEL=(X3**2*(Y1-Y2)+X1**2*(Y2-Y3))/(X1*X3*(X3-X1)) 
      AEL=(Y1-Y2-BEL*X1)/(X1**2)
      CEL=Y2
      IF (AEL.EQ.0..AND.BEL.EQ.0..AND.CEL.EQ.0.) CEL=1
      F1=(YNODE(1,IC)-DL)/H1(I)
      F2=(YNODE(1,I)-DL)/H1(I)
      F3=(YNODE(2,I)-DL)/H1(I)
      F4=(YNODE(2,IC)-DL)/H1(I)
      A1=AVEL*AEL/5
      A2=(AVEL*BEL+AEL*BVEL)/4
      A3=(AVEL*CEL+AEL*CVEL+BVEL*BEL)/3
      A4=(BVEL*CEL+BEL*CVEL)/2
      A5=CVEL*CEL
      CNUM1=A1*(F1**5-F4**5)+A2*(F1**4-F4**4)+A3*(F1**3-F4**3)+
     *A4*(F1**2-F4**2)+A5*(F1-F4)
      CNUM2=A1*(F2**5-F3**5)+A2*(F2**4-F3**4)+A3*(F2**3-F3**3)+
     *A4*(F2**2-F3**2)+A5*(F2-F3)
      CNUM=CNUM2/CNUM1
      APOI(JNUM+1)=CNUM*AAUU(IC)
      APOI(JNUM+2)=-CNUM*AAUU(IC)
      APOI(JNUM+3)=-AAVV(I)
      APOI(JNUM+4)=-AAUU(I)
      APOI(JNUM+5)=DJ(I)+AAVV(I)+AAUU(I)+AAVV(NJP1(I))
      APOI(JNUM+6)=-AAVV(NJP1(I))
      DO K=1,KBM1
          RHSPOI(I)=RHSPOI(I)+DZ(K)*
     *    (CNUM*XMFLUX(IC,K)-XMFLUX(I,K)+
     *    YMFLUX(NJP1(I),K)-YMFLUX(I,K))
      ENDDO
      RHSPOI(I)=DJ(I)*EL(I)-RHSPOI(I)*DTI
      endif
22    ENDDO
c!$OMP DISTRIBUTE PARALLEL DO
c      DO 23 II=1,NNIVGFN
      DO concurrent (II=1:NNIVGFN) local(I,JNUM,IC,DL,X1,X3,Y1,Y2,Y3,
     * AEL,BEL,CEL,F1,F2,F3,F4,A1,A2,A3,A4,A5,CNUM1,CNUM2,CNUM,K)
     * local_init(KBM1,DTI)
      I=NIVGFN(II)
!      IF (AIJ(I,1).EQ.1) GO TO 23
      if (AIJ(I,1)/=1) then
      JNUM=IAPOI(I)-1
      IC=AIJ(I,3)
      DL=(XNODE(1,IC)+XNODE(4,IC))/2
      X1=((XNODE(1,NIM1(IC))+XNODE(4,NIM1(IC)))/2-DL)/H2(I)
!      X2=0.
      X3=((XNODE(1,NIP1(IC))+XNODE(4,NIP1(IC)))/2-DL)/H2(I)
     !! Y1=V(NIM1(IC))
     !! Y2=V(IC)
     !! Y3=V(NIP1(IC))
     !! BVEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/
     !!*(X2-X3)/((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2))
     !! AVEL=(Y1-Y2-BVEL*(X1-X2))/(X1**2-X2**2)
     !! CVEL=Y1-AVEL*X1**2-BVEL*X1
     !! IF (AVEL.EQ.0..AND.BVEL.EQ.0..AND.CVEL.EQ.0.) CVEL=1
      Y1=DV(NIM1(IC))*DVM(NIM1(IC))+DV(IC)*(1-DVM(NIM1(IC)))
      Y2=DV(IC)*DVM(IC)
      Y3=DV(NIP1(IC))*DVM(NIP1(IC))+DV(IC)*(1-DVM(NIP1(IC)))
!      BEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/(X2-X3)/ 
!     *((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2)) 
!      AEL=(Y1-Y2-BEL*(X1-X2))/(X1**2-X2**2)
!      CEL=Y1-AEL*X1**2-BEL*X1
      BEL=(X3**2*(Y1-Y2)+X1**2*(Y2-Y3))/(X1*X3*(X3-X1)) 
      AEL=(Y1-Y2-BEL*X1)/(X1**2)
      CEL=Y2
      IF (AEL.EQ.0..AND.BEL.EQ.0..AND.CEL.EQ.0.) CEL=1
      F1=(XNODE(1,IC)-DL)/H2(I)
      F2=(XNODE(1,I)-DL)/H2(I)
      F3=(XNODE(4,I)-DL)/H2(I)
      F4=(XNODE(4,IC)-DL)/H2(I)
      A1=AVEL*AEL/5
      A2=(AVEL*BEL+AEL*BVEL)/4
      A3=(AVEL*CEL+AEL*CVEL+BVEL*BEL)/3
      A4=(BVEL*CEL+BEL*CVEL)/2
      A5=CVEL*CEL
      CNUM1=A1*(F1**5-F4**5)+A2*(F1**4-F4**4)+A3*(F1**3-F4**3)+
     *A4*(F1**2-F4**2)+A5*(F1-F4)
      CNUM2=A1*(F2**5-F3**5)+A2*(F2**4-F3**4)+A3*(F2**3-F3**3)+
     *A4*(F2**2-F3**2)+A5*(F2-F3)
      CNUM=CNUM2/CNUM1
      APOI(JNUM+1)=CNUM*AAVV(IC)
      APOI(JNUM+2)=-CNUM*AAVV(IC)
      APOI(JNUM+3)=-AAVV(I)
      APOI(JNUM+4)=-AAUU(I)
      APOI(JNUM+5)=DJ(I)+AAUU(I)+AAUU(NIP1(I))+AAVV(I)
      APOI(JNUM+6)=-AAUU(NIP1(I))
      DO K=1,KBM1
          RHSPOI(I)=RHSPOI(I)+DZ(K)*
     *    (XMFLUX(NIP1(I),K)-XMFLUX(I,K)+
     *    CNUM*YMFLUX(IC,K)-YMFLUX(I,K))
      ENDDO
      RHSPOI(I)=DJ(I)*EL(I)-RHSPOI(I)*DTI
      endif
23    ENDDO
c!$OMP DISTRIBUTE PARALLEL DO
c      DO 24 II=1,NNIVGFS
      DO concurrent (II=1:NNIVGFS) local(I,JNUM,IC,DL,X1,X3,Y1,Y2,Y3,
     * AEL,BEL,CEL,F1,F2,F3,F4,A1,A2,A3,A4,A5,CNUM1,CNUM2,CNUM,K)
     * local_init(KBM1,DTI)
      I=NIVGFS(II)
!      IF (AIJ(I,1).EQ.1) GO TO 24
      if (AIJ(I,1)/=1) then
      JNUM=IAPOI(I)-1
      IC=AIJ(I,3)
      DL=(XNODE(1,IC)+XNODE(4,IC))/2
      X1=((XNODE(1,NIM1(IC))+XNODE(4,NIM1(IC)))/2-DL)/H2(I)
!      X2=0.
      X3=((XNODE(1,NIP1(IC))+XNODE(4,NIP1(IC)))/2-DL)/H2(I)
     !! Y1=V(NIM1(IC))
     !! Y2=V(IC)
     !! Y3=V(NIP1(IC))
     !! BVEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/
     !!*(X2-X3)/((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2))
     !! AVEL=(Y1-Y2-BVEL*(X1-X2))/(X1**2-X2**2)
     !! CVEL=Y1-AVEL*X1**2-BVEL*X1
     !! IF (AVEL.EQ.0..AND.BVEL.EQ.0..AND.CVEL.EQ.0.) CVEL=1
      Y1=DV(NIM1(IC))*DVM(NIM1(IC))+DV(IC)*(1-DVM(NIM1(IC)))
      Y2=DV(IC)*DVM(IC)
      Y3=DV(NIP1(IC))*DVM(NIP1(IC))+DV(IC)*(1-DVM(NIP1(IC)))
!      BEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/(X2-X3)/ 
!     *((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2)) 
!      AEL=(Y1-Y2-BEL*(X1-X2))/(X1**2-X2**2)
!      CEL=Y1-AEL*X1**2-BEL*X1
      BEL=(X3**2*(Y1-Y2)+X1**2*(Y2-Y3))/(X1*X3*(X3-X1)) 
      AEL=(Y1-Y2-BEL*X1)/(X1**2)
      CEL=Y2
      IF (AEL.EQ.0..AND.BEL.EQ.0..AND.CEL.EQ.0.) CEL=1
      F1=(XNODE(1,IC)-DL)/H2(I)
      F2=(XNODE(1,I)-DL)/H2(I)
      F3=(XNODE(4,I)-DL)/H2(I)
      F4=(XNODE(4,IC)-DL)/H2(I)
      A1=AVEL*AEL/5
      A2=(AVEL*BEL+AEL*BVEL)/4
      A3=(AVEL*CEL+AEL*CVEL+BVEL*BEL)/3
      A4=(BVEL*CEL+BEL*CVEL)/2
      A5=CVEL*CEL
      CNUM1=A1*(F1**5-F4**5)+A2*(F1**4-F4**4)+A3*(F1**3-F4**3)+
     *A4*(F1**2-F4**2)+A5*(F1-F4)
      CNUM2=A1*(F2**5-F3**5)+A2*(F2**4-F3**4)+A3*(F2**3-F3**3)+
     *A4*(F2**2-F3**2)+A5*(F2-F3)
      CNUM=CNUM2/CNUM1
      APOI(JNUM+1)=-CNUM*AAVV(IC)
      APOI(JNUM+2)=CNUM*AAVV(IC)
      APOI(JNUM+3)=-AAUU(I)
      APOI(JNUM+4)=DJ(I)+AAUU(I)+AAUU(NIP1(I))+AAVV(NJP1(I))
      APOI(JNUM+5)=-AAUU(NIP1(I))
      APOI(JNUM+6)=-AAVV(NJP1(I))
      DO K=1,KBM1
          RHSPOI(I)=RHSPOI(I)+DZ(K)*
     *    (XMFLUX(NIP1(I),K)-XMFLUX(I,K)+
     *    YMFLUX(NJP1(I),K)-CNUM*YMFLUX(IC,K))
      ENDDO
      RHSPOI(I)=DJ(I)*EL(I)-RHSPOI(I)*DTI
      endif
24    ENDDO
c!$OMP END TARGET teams
      
      
c!$OMP TARGET teams DEFAULTMAP(present: allocatable)
c!$OMP DISTRIBUTE PARALLEL DO
c      DO 25 II=1,NNIVGFWS
      DO concurrent (II=1:NNIVGFWS) local(I,JNUM,IC,DL,X1,X3,Y1,Y2,Y3,
     * AEL,BEL,CEL,F1,F2,F3,F4,A1,A2,A3,A4,A5,CNUM1,CNUM2,CNUM,K)
     * local_init(KBM1,DTI)
      I=NIVGFWS(II)
!      IF (AIJ(I,1).EQ.1) GO TO 25
      if (AIJ(I,1)/=1) then
      JNUM=IAPOI(I)-1
      IC=AIJ(I,4)
      DL=(YNODE(1,IC)+YNODE(2,IC))/2
      X1=((YNODE(1,NJM1(IC))+YNODE(2,NJM1(IC)))/2-DL)/H1(I)
!      X2=0.0
      X3=((YNODE(1,NJP1(IC))+YNODE(2,NJP1(IC)))/2-DL)/H1(I)
     !! Y1=U(NJM1(IC))
     !! Y2=U(IC)
     !! Y3=U(NJP1(IC))
     !! BVEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/
     !!*(X2-X3)/((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2))
     !! AVEL=(Y1-Y2-BVEL*(X1-X2))/(X1**2-X2**2)
     !! CVEL=Y1-AVEL*X1**2-BVEL*X1
     !! IF (AVEL.EQ.0..AND.BVEL.EQ.0..AND.CVEL.EQ.0.) CVEL=1
      Y1=DU(NJM1(IC))*DUM(NJM1(IC))+DU(IC)*(1-DUM(NJM1(IC)))
      Y2=DU(IC)*DUM(IC)
      Y3=DU(NJP1(IC))*DUM(NJP1(IC))+DU(IC)*(1-DUM(NJP1(IC)))
!      BEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/(X2-X3)/ 
!     *((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2)) 
!      AEL=(Y1-Y2-BEL*(X1-X2))/(X1**2-X2**2)
!      CEL=Y1-AEL*X1**2-BEL*X1
      BEL=(X3**2*(Y1-Y2)+X1**2*(Y2-Y3))/(X1*X3*(X3-X1)) 
      AEL=(Y1-Y2-BEL*X1)/(X1**2)
      CEL=Y2
      IF (AEL.EQ.0..AND.BEL.EQ.0..AND.CEL.EQ.0.) CEL=1
      F1=(YNODE(1,IC)-DL)/H1(I)
      F2=(YNODE(1,I)-DL)/H1(I)
      F3=(YNODE(2,I)-DL)/H1(I)
      F4=(YNODE(2,IC)-DL)/H1(I)
      A1=AVEL*AEL/5
      A2=(AVEL*BEL+AEL*BVEL)/4
      A3=(AVEL*CEL+AEL*CVEL+BVEL*BEL)/3
      A4=(BVEL*CEL+BEL*CVEL)/2
      A5=CVEL*CEL
      CNUM=(AVEL*AEL*(F2**5-F3**5)/5.+
     *(AVEL*BEL+AEL*BVEL)*(F2**4-F3**4)/4.+
     *(AVEL*CEL+AEL*CVEL+BVEL*BEL)*(F2**3-F3**3)/3.+
     *(BVEL*CEL+BEL*CVEL)*(F2**2-F3**2)/2.+
     *CVEL*CEL*(F2-F3))/
     *(AVEL*AEL*(F1**5-F4**5)/5.+
     *(AVEL*BEL+AEL*BVEL)*(F1**4-F4**4)/4.+
     *(AVEL*CEL+AEL*CVEL+BVEL*BEL)*(F1**3-F4**3)/3.+
     *(BVEL*CEL+BEL*CVEL)*(F1**2-F4**2)/2.+
     *CVEL*CEL*(F1-F4))
      DL=(XNODE(1,IC)+XNODE(4,IC))/2
      X1=((XNODE(1,NIM1(IC))+XNODE(4,NIM1(IC)))/2-DL)/H2(I)
!      X2=0.
      X3=((XNODE(1,NIP1(IC))+XNODE(4,NIP1(IC)))/2-DL)/H2(I)
     !! Y1=V(NIM1(IC))
     !! Y2=V(IC)
     !! Y3=V(NIP1(IC))
     !! BVEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/
     !!*(X2-X3)/((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2))
     !! AVEL=(Y1-Y2-BVEL*(X1-X2))/(X1**2-X2**2)
     !! CVEL=Y1-AVEL*X1**2-BVEL*X1
     !! IF (AVEL.EQ.0..AND.BVEL.EQ.0..AND.CVEL.EQ.0.) CVEL=1
      Y1=DV(NIM1(IC))*DVM(NIM1(IC))+DV(IC)*(1-DVM(NIM1(IC)))
      Y2=DV(IC)*DVM(IC)
      Y3=DV(NIP1(IC))*DVM(NIP1(IC))+DV(IC)*(1-DVM(NIP1(IC)))
!      BEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/(X2-X3)/ 
!     *((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2)) 
!      AEL=(Y1-Y2-BEL*(X1-X2))/(X1**2-X2**2)
!      CEL=Y1-AEL*X1**2-BEL*X1
      BEL=(X3**2*(Y1-Y2)+X1**2*(Y2-Y3))/(X1*X3*(X3-X1)) 
      AEL=(Y1-Y2-BEL*X1)/(X1**2)
      CEL=Y2
      IF (AEL.EQ.0..AND.BEL.EQ.0..AND.CEL.EQ.0.) CEL=1
      F1=(XNODE(1,IC)-DL)/H2(I)
      F2=(XNODE(1,I)-DL)/H2(I)
      F3=(XNODE(4,I)-DL)/H2(I)
      F4=(XNODE(4,IC)-DL)/H2(I)
      A1=AVEL*AEL/5
      A2=(AVEL*BEL+AEL*BVEL)/4
      A3=(AVEL*CEL+AEL*CVEL+BVEL*BEL)/3
      A4=(BVEL*CEL+BEL*CVEL)/2
      A5=CVEL*CEL
      DNUM=(AVEL*AEL*(F2**5-F3**5)/5.+
     *(AVEL*BEL+AEL*BVEL)*(F2**4-F3**4)/4.+
     *(AVEL*CEL+AEL*CVEL+BVEL*BEL)*(F2**3-F3**3)/3.+
     *(BVEL*CEL+BEL*CVEL)*(F2**2-F3**2)/2.+
     *CVEL*CEL*(F2-F3))/
     *(AVEL*AEL*(F1**5-F4**5)/5.+
     *(AVEL*BEL+AEL*BVEL)*(F1**4-F4**4)/4.+
     *(AVEL*CEL+AEL*CVEL+BVEL*BEL)*(F1**3-F4**3)/3.+
     *(BVEL*CEL+BEL*CVEL)*(F1**2-F4**2)/2.+
     *CVEL*CEL*(F1-F4))
      APOI(JNUM+1)=-DNUM*AAVV(IC)
      APOI(JNUM+2)=-CNUM*AAUU(IC)
      APOI(JNUM+3)=DNUM*AAVV(IC)+CNUM*AAUU(IC)
      APOI(JNUM+4)=DJ(I)+AAUU(NIP1(I))+AAVV(NJP1(I))
      APOI(JNUM+5)=-AAUU(NIP1(I))
      APOI(JNUM+6)=-AAVV(NJP1(I))
      DO K=1,KBM1
          RHSPOI(I)=RHSPOI(I)+DZ(K)*
     *    (XMFLUX(NIP1(I),K)-CNUM*XMFLUX(IC,K)+
     *    YMFLUX(NJP1(I),K)-DNUM*YMFLUX(IC,K))
      ENDDO
      RHSPOI(I)=DJ(I)*EL(I)-RHSPOI(I)*DTI
      endif
25    ENDDO
c!$OMP DISTRIBUTE PARALLEL DO
c      DO 26 II=1,NNIVGFWN
      DO concurrent (II=1:NNIVGFWN) local(I,JNUM,IC,DL,X1,X3,Y1,Y2,Y3,
     * AEL,BEL,CEL,F1,F2,F3,F4,A1,A2,A3,A4,A5,CNUM1,CNUM2,CNUM,K)
     * local_init(KBM1,DTI)
      I=NIVGFWN(II)
!      IF (AIJ(I,1).EQ.1) GO TO 26
      if (AIJ(I,1)/=1) then
      JNUM=IAPOI(I)-1
      IC=AIJ(I,4)
      DL=(XNODE(1,IC)+XNODE(4,IC))/2
      X1=((XNODE(1,NIM1(IC))+XNODE(4,NIM1(IC)))/2-DL)/H2(I)
!      X2=0.
      X3=((XNODE(1,NIP1(IC))+XNODE(4,NIP1(IC)))/2-DL)/H2(I)
     !! Y1=V(NIM1(IC))
     !! Y2=V(IC)
     !! Y3=V(NIP1(IC))
     !! BVEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/
     !!*(X2-X3)/((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2))
     !! AVEL=(Y1-Y2-BVEL*(X1-X2))/(X1**2-X2**2)
     !! CVEL=Y1-AVEL*X1**2-BVEL*X1
     !! IF (AVEL.EQ.0..AND.BVEL.EQ.0..AND.CVEL.EQ.0.) CVEL=1
      Y1=DV(NIM1(IC))*DVM(NIM1(IC))+DV(IC)*(1-DVM(NIM1(IC)))
      Y2=DV(IC)*DVM(IC)
      Y3=DV(NIP1(IC))*DVM(NIP1(IC))+DV(IC)*(1-DVM(NIP1(IC)))
!      BEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/(X2-X3)/ 
!     *((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2)) 
!      AEL=(Y1-Y2-BEL*(X1-X2))/(X1**2-X2**2)
!      CEL=Y1-AEL*X1**2-BEL*X1
      BEL=(X3**2*(Y1-Y2)+X1**2*(Y2-Y3))/(X1*X3*(X3-X1)) 
      AEL=(Y1-Y2-BEL*X1)/(X1**2)
      CEL=Y2
      IF (AEL.EQ.0..AND.BEL.EQ.0..AND.CEL.EQ.0.) CEL=1
      F1=(XNODE(1,IC)-DL)/H2(I)
      F2=(XNODE(1,I)-DL)/H2(I)
      F3=(XNODE(4,I)-DL)/H2(I)
      F4=(XNODE(4,IC)-DL)/H2(I)
      DNUM=(AVEL*AEL*(F2**5-F3**5)/5.+
     *(AVEL*BEL+AEL*BVEL)*(F2**4-F3**4)/4.+
     *(AVEL*CEL+AEL*CVEL+BVEL*BEL)*(F2**3-F3**3)/3.+
     *(BVEL*CEL+BEL*CVEL)*(F2**2-F3**2)/2.+
     *CVEL*CEL*(F2-F3))/
     *(AVEL*AEL*(F1**5-F4**5)/5.+
     *(AVEL*BEL+AEL*BVEL)*(F1**4-F4**4)/4.+
     *(AVEL*CEL+AEL*CVEL+BVEL*BEL)*(F1**3-F4**3)/3.+
     *(BVEL*CEL+BEL*CVEL)*(F1**2-F4**2)/2.+
     *CVEL*CEL*(F1-F4))
      IC=AIJ(I,3)
      DL=(YNODE(1,IC)+YNODE(2,IC))/2
      X1=((YNODE(1,NJM1(IC))+YNODE(2,NJM1(IC)))/2-DL)/H1(I)
!      X2=0.0
      X3=((YNODE(1,NJP1(IC))+YNODE(2,NJP1(IC)))/2-DL)/H1(I)
     !! Y1=U(NJM1(IC))
     !! Y2=U(IC)
     !! Y3=U(NJP1(IC))
     !! BVEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/
     !!*(X2-X3)/((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2))
     !! AVEL=(Y1-Y2-BVEL*(X1-X2))/(X1**2-X2**2)
     !! CVEL=Y1-AVEL*X1**2-BVEL*X1
     !! IF (AVEL.EQ.0..AND.BVEL.EQ.0..AND.CVEL.EQ.0.) CVEL=1
      Y1=DU(NJM1(IC))*DUM(NJM1(IC))+DU(IC)*(1-DUM(NJM1(IC)))
      Y2=DU(IC)*DUM(IC)
      Y3=DU(NJP1(IC))*DUM(NJP1(IC))+DU(IC)*(1-DUM(NJP1(IC)))
!      BEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/(X2-X3)/ 
!     *((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2)) 
!      AEL=(Y1-Y2-BEL*(X1-X2))/(X1**2-X2**2)
!      CEL=Y1-AEL*X1**2-BEL*X1
      BEL=(X3**2*(Y1-Y2)+X1**2*(Y2-Y3))/(X1*X3*(X3-X1)) 
      AEL=(Y1-Y2-BEL*X1)/(X1**2)
      CEL=Y2
      IF (AEL.EQ.0..AND.BEL.EQ.0..AND.CEL.EQ.0.) CEL=1
      F1=(YNODE(1,IC)-DL)/H1(I)
      F2=(YNODE(1,I)-DL)/H1(I)
      F3=(YNODE(2,I)-DL)/H1(I)
      F4=(YNODE(2,IC)-DL)/H1(I)
      CNUM=(AVEL*AEL*(F2**5-F3**5)/5.+
     *(AVEL*BEL+AEL*BVEL)*(F2**4-F3**4)/4.+
     *(AVEL*CEL+AEL*CVEL+BVEL*BEL)*(F2**3-F3**3)/3.+
     *(BVEL*CEL+BEL*CVEL)*(F2**2-F3**2)/2.+
     *CVEL*CEL*(F2-F3))/
     *(AVEL*AEL*(F1**5-F4**5)/5.+
     *(AVEL*BEL+AEL*BVEL)*(F1**4-F4**4)/4.+
     *(AVEL*CEL+AEL*CVEL+BVEL*BEL)*(F1**3-F4**3)/3.+
     *(BVEL*CEL+BEL*CVEL)*(F1**2-F4**2)/2.+
     *CVEL*CEL*(F1-F4))
      APOI(JNUM+1)=-CNUM*AAUU(IC)
      APOI(JNUM+2)=DNUM*AAVV(NJP1(IC))+CNUM*AAUU(IC)
      APOI(JNUM+3)=-DNUM*AAVV(NJP1(IC))
      APOI(JNUM+4)=-AAVV(I)
      APOI(JNUM+5)=DJ(I)+AAVV(I)+AAUU(NIP1(I))
      APOI(JNUM+6)=-AAUU(NIP1(I))
      DO K=1,KBM1
          RHSPOI(I)=RHSPOI(I)+DZ(K)*
     *    (XMFLUX(NIP1(I),K)-CNUM*XMFLUX(IC,K)+
     *    DNUM*YMFLUX(NJP1(IC),K)-YMFLUX(I,K))
      ENDDO
      RHSPOI(I)=DJ(I)*EL(I)-RHSPOI(I)*DTI
      endif
26    ENDDO
c!$OMP DISTRIBUTE PARALLEL DO
c      DO 27 II=1,NNIVGFES
      DO concurrent (II=1:NNIVGFES) local(I,JNUM,IC,DL,X1,X3,Y1,Y2,Y3,
     * AEL,BEL,CEL,F1,F2,F3,F4,A1,A2,A3,A4,A5,CNUM1,CNUM2,CNUM,K)
     * local_init(KBM1,DTI)
      I=NIVGFES(II)
!      IF (AIJ(I,1).EQ.1) GO TO 27
      if (AIJ(I,1)/=1) then
      JNUM=IAPOI(I)-1
      IC=AIJ(I,4)
      DL=(YNODE(1,IC)+YNODE(2,IC))/2
      X1=((YNODE(1,NJM1(IC))+YNODE(2,NJM1(IC)))/2-DL)/H1(I)
!      X2=0.0
      X3=((YNODE(1,NJP1(IC))+YNODE(2,NJP1(IC)))/2-DL)/H1(I)
     !! Y1=U(NJM1(IC))
     !! Y2=U(IC)
     !! Y3=U(NJP1(IC))
     !! BVEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/
     !!*(X2-X3)/((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2))
     !! AVEL=(Y1-Y2-BVEL*(X1-X2))/(X1**2-X2**2)
     !! CVEL=Y1-AVEL*X1**2-BVEL*X1
     !! IF (AVEL.EQ.0..AND.BVEL.EQ.0..AND.CVEL.EQ.0.) CVEL=1
      Y1=DU(NJM1(IC))*DUM(NJM1(IC))+DU(IC)*(1-DUM(NJM1(IC)))
      Y2=DU(IC)*DUM(IC)
      Y3=DU(NJP1(IC))*DUM(NJP1(IC))+DU(IC)*(1-DUM(NJP1(IC)))
!      BEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/(X2-X3)/ 
!     *((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2)) 
!      AEL=(Y1-Y2-BEL*(X1-X2))/(X1**2-X2**2)
!      CEL=Y1-AEL*X1**2-BEL*X1
      BEL=(X3**2*(Y1-Y2)+X1**2*(Y2-Y3))/(X1*X3*(X3-X1)) 
      AEL=(Y1-Y2-BEL*X1)/(X1**2)
      CEL=Y2
      IF (AEL.EQ.0..AND.BEL.EQ.0..AND.CEL.EQ.0.) CEL=1
      F1=(YNODE(1,IC)-DL)/H1(I)
      F2=(YNODE(1,I)-DL)/H1(I)
      F3=(YNODE(2,I)-DL)/H1(I)
      F4=(YNODE(2,IC)-DL)/H1(I)
      A1=AVEL*AEL/5
      A2=(AVEL*BEL+AEL*BVEL)/4
      A3=(AVEL*CEL+AEL*CVEL+BVEL*BEL)/3
      A4=(BVEL*CEL+BEL*CVEL)/2
      A5=CVEL*CEL
      CNUM1=A1*(F1**5-F4**5)+A2*(F1**4-F4**4)+A3*(F1**3-F4**3)+
     *A4*(F1**2-F4**2)+A5*(F1-F4)
      CNUM2=A1*(F2**5-F3**5)+A2*(F2**4-F3**4)+A3*(F2**3-F3**3)+
     *A4*(F2**2-F3**2)+A5*(F2-F3)
      CNUM=CNUM2/CNUM1
      IC=AIJ(I,3)
      DL=(XNODE(1,IC)+XNODE(4,IC))/2
      X1=((XNODE(1,NIM1(IC))+XNODE(4,NIM1(IC)))/2-DL)/H2(I)
!      X2=0.
      X3=((XNODE(1,NIP1(IC))+XNODE(4,NIP1(IC)))/2-DL)/H2(I)
     !! Y1=V(NIM1(IC))
     !! Y2=V(IC)
     !! Y3=V(NIP1(IC))
     !! BVEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/
     !!*(X2-X3)/((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2))
     !! AVEL=(Y1-Y2-BVEL*(X1-X2))/(X1**2-X2**2)
     !! CVEL=Y1-AVEL*X1**2-BVEL*X1
     !! IF (AVEL.EQ.0..AND.BVEL.EQ.0..AND.CVEL.EQ.0.) CVEL=1
      Y1=DV(NIM1(IC))*DVM(NIM1(IC))+DV(IC)*(1-DVM(NIM1(IC)))
      Y2=DV(IC)*DVM(IC)
      Y3=DV(NIP1(IC))*DVM(NIP1(IC))+DV(IC)*(1-DVM(NIP1(IC)))
!      BEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/(X2-X3)/ 
!     *((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2)) 
!      AEL=(Y1-Y2-BEL*(X1-X2))/(X1**2-X2**2)
!      CEL=Y1-AEL*X1**2-BEL*X1
      BEL=(X3**2*(Y1-Y2)+X1**2*(Y2-Y3))/(X1*X3*(X3-X1)) 
      AEL=(Y1-Y2-BEL*X1)/(X1**2)
      CEL=Y2
      IF (AEL.EQ.0..AND.BEL.EQ.0..AND.CEL.EQ.0.) CEL=1
      F1=(XNODE(1,IC)-DL)/H2(I)
      F2=(XNODE(1,I)-DL)/H2(I)
      F3=(XNODE(4,I)-DL)/H2(I)
      F4=(XNODE(4,IC)-DL)/H2(I)
      A1=AVEL*AEL/5
      A2=(AVEL*BEL+AEL*BVEL)/4
      A3=(AVEL*CEL+AEL*CVEL+BVEL*BEL)/3
      A4=(BVEL*CEL+BEL*CVEL)/2
      A5=CVEL*CEL
      CNUM1=A1*(F1**5-F4**5)+A2*(F1**4-F4**4)+A3*(F1**3-F4**3)+
     *A4*(F1**2-F4**2)+A5*(F1-F4)
      CNUM2=A1*(F2**5-F3**5)+A2*(F2**4-F3**4)+A3*(F2**3-F3**3)+
     *A4*(F2**2-F3**2)+A5*(F2-F3)
      DNUM=CNUM2/CNUM1
      APOI(JNUM+1)=-DNUM*AAVV(IC)
      APOI(JNUM+2)=CNUM*AAUU(NIP1(IC))+DNUM*AAVV(IC)
      APOI(JNUM+3)=-CNUM*AAUU(NIP1(IC))
      APOI(JNUM+4)=-AAUU(I)
      APOI(JNUM+5)=DJ(I)+AAUU(I)+AAVV(NJP1(I))
      APOI(JNUM+6)=-AAVV(NJP1(I))
      DO K=1,KBM1
          RHSPOI(I)=RHSPOI(I)+DZ(K)*
     *    (CNUM*XMFLUX(NIP1(IC),K)-XMFLUX(I,K)+
     *    YMFLUX(NJP1(I),K)-DNUM*YMFLUX(IC,K))
      ENDDO
      RHSPOI(I)=DJ(I)*EL(I)-RHSPOI(I)*DTI
      endif
27    ENDDO
c!$OMP DISTRIBUTE PARALLEL DO
c      DO 28 II=1,NNIVGFEN
      DO concurrent (II=1:NNIVGFEN) local(I,JNUM,IC,DL,X1,X3,Y1,Y2,Y3,
     * AEL,BEL,CEL,F1,F2,F3,F4,A1,A2,A3,A4,A5,CNUM1,CNUM2,CNUM,K)
     * local_init(KBM1,DTI)
      I=NIVGFEN(II)
!      IF (AIJ(I,1).EQ.1) GO TO 28
      if (AIJ(I,1)/=1) then
      JNUM=IAPOI(I)-1
      IC=AIJ(I,3)
      DL=(YNODE(1,IC)+YNODE(2,IC))/2
      X1=((YNODE(1,NJM1(IC))+YNODE(2,NJM1(IC)))/2-DL)/H1(I)
!      X2=0.0
      X3=((YNODE(1,NJP1(IC))+YNODE(2,NJP1(IC)))/2-DL)/H1(I)
     !! Y1=U(NJM1(IC))
     !! Y2=U(IC)
     !! Y3=U(NJP1(IC))
     !! BVEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/
     !!*(X2-X3)/((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2))
     !! AVEL=(Y1-Y2-BVEL*(X1-X2))/(X1**2-X2**2)
     !! CVEL=Y1-AVEL*X1**2-BVEL*X1
     !! IF (AVEL.EQ.0..AND.BVEL.EQ.0..AND.CVEL.EQ.0.) CVEL=1
      Y1=DU(NJM1(IC))*DUM(NJM1(IC))+DU(IC)*(1-DUM(NJM1(IC)))
      Y2=DU(IC)*DUM(IC)
      Y3=DU(NJP1(IC))*DUM(NJP1(IC))+DU(IC)*(1-DUM(NJP1(IC)))
!      BEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/(X2-X3)/ 
!     *((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2)) 
!      AEL=(Y1-Y2-BEL*(X1-X2))/(X1**2-X2**2)
!      CEL=Y1-AEL*X1**2-BEL*X1
      BEL=(X3**2*(Y1-Y2)+X1**2*(Y2-Y3))/(X1*X3*(X3-X1)) 
      AEL=(Y1-Y2-BEL*X1)/(X1**2)
      CEL=Y2
      IF (AEL.EQ.0..AND.BEL.EQ.0..AND.CEL.EQ.0.) CEL=1
      F1=(YNODE(1,IC)-DL)/H1(I)
      F2=(YNODE(1,I)-DL)/H1(I)
      F3=(YNODE(2,I)-DL)/H1(I)
      F4=(YNODE(2,IC)-DL)/H1(I)
      A1=AVEL*AEL/5
      A2=(AVEL*BEL+AEL*BVEL)/4
      A3=(AVEL*CEL+AEL*CVEL+BVEL*BEL)/3
      A4=(BVEL*CEL+BEL*CVEL)/2
      A5=CVEL*CEL
      CNUM1=A1*(F1**5-F4**5)+A2*(F1**4-F4**4)+A3*(F1**3-F4**3)+
     *A4*(F1**2-F4**2)+A5*(F1-F4)
      CNUM2=A1*(F2**5-F3**5)+A2*(F2**4-F3**4)+A3*(F2**3-F3**3)+
     *A4*(F2**2-F3**2)+A5*(F2-F3)
      CNUM=CNUM2/CNUM1
      IC=AIJ(I,4)
      DL=(XNODE(1,IC)+XNODE(4,IC))/2
      X1=((XNODE(1,NIM1(IC))+XNODE(4,NIM1(IC)))/2-DL)/H2(I)
!      X2=0.
      X3=((XNODE(1,NIP1(IC))+XNODE(4,NIP1(IC)))/2-DL)/H2(I)
     !! Y1=V(NIM1(IC))
     !! Y2=V(IC)
     !! Y3=V(NIP1(IC))
     !! BVEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/
     !!*(X2-X3)/((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2))
     !! AVEL=(Y1-Y2-BVEL*(X1-X2))/(X1**2-X2**2)
     !! CVEL=Y1-AVEL*X1**2-BVEL*X1
     !! IF (AVEL.EQ.0..AND.BVEL.EQ.0..AND.CVEL.EQ.0.) CVEL=1
      Y1=DV(NIM1(IC))*DVM(NIM1(IC))+DV(IC)*(1-DVM(NIM1(IC)))
      Y2=DV(IC)*DVM(IC)
      Y3=DV(NIP1(IC))*DVM(NIP1(IC))+DV(IC)*(1-DVM(NIP1(IC)))
!      BEL=((X2**2-X3**2)*(Y1-Y2)-(X1**2-X2**2)*(Y2-Y3))/(X1-X2)/(X2-X3)/ 
!     *((X2**2-X3**2)/(X2-X3)-(X1**2-X2**2)/(X1-X2)) 
!      AEL=(Y1-Y2-BEL*(X1-X2))/(X1**2-X2**2)
!      CEL=Y1-AEL*X1**2-BEL*X1
      BEL=(X3**2*(Y1-Y2)+X1**2*(Y2-Y3))/(X1*X3*(X3-X1)) 
      AEL=(Y1-Y2-BEL*X1)/(X1**2)
      CEL=Y2
      IF (AEL.EQ.0..AND.BEL.EQ.0..AND.CEL.EQ.0.) CEL=1
      F1=(XNODE(1,IC)-DL)/H2(I)
      F2=(XNODE(1,I)-DL)/H2(I)
      F3=(XNODE(4,I)-DL)/H2(I)
      F4=(XNODE(4,IC)-DL)/H2(I)
      A1=AVEL*AEL/5
      A2=(AVEL*BEL+AEL*BVEL)/4
      A3=(AVEL*CEL+AEL*CVEL+BVEL*BEL)/3
      A4=(BVEL*CEL+BEL*CVEL)/2
      A5=CVEL*CEL
      CNUM1=A1*(F1**5-F4**5)+A2*(F1**4-F4**4)+A3*(F1**3-F4**3)+
     *A4*(F1**2-F4**2)+A5*(F1-F4)
      CNUM2=A1*(F2**5-F3**5)+A2*(F2**4-F3**4)+A3*(F2**3-F3**3)+
     *A4*(F2**2-F3**2)+A5*(F2-F3)
      DNUM=CNUM2/CNUM1
      IC=AIJ(I,2)
      APOI(JNUM+1)=CNUM*AAUU(NIP1(IC))+DNUM*AAVV(NJP1(IC))
      APOI(JNUM+2)=-CNUM*AAUU(NIP1(IC)) 
      APOI(JNUM+3)=-DNUM*AAVV(NJP1(IC))
      APOI(JNUM+4)=-AAVV(I)
      APOI(JNUM+5)=-AAUU(I)
      APOI(JNUM+6)=DJ(I)+AAUU(I)+AAVV(I)
      DO K=1,KBM1
      RHSPOI(I)=RHSPOI(I)+DZ(K)*
     *(CNUM*XMFLUX(NIP1(IC),K)-XMFLUX(I,K)+
     *DNUM*YMFLUX(NJP1(IC),K)-YMFLUX(I,K))
      ENDDO
      RHSPOI(I)=DJ(I)*EL(I)-RHSPOI(I)*DTI
      endif
28    ENDDO
c!$OMP END TARGET  teams 

      
c --- Merged:      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I=1,NPOI
c          XPOI(I)=0.
c      ENDDO
c!$OMP END TARGET 
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
cc      DO I=1,N_CTRD_AG+N_CTRD_EVG
c      DO I=1,NPOI
      DO concurrent (I=1:NPOI)
          IF (FSMADD(I).GT.0.) THEN
              XPOI(I)=ELF(I)
          ELSE
              XPOI(I)=0.
          ENDIF
      ENDDO
c!$OMP END TARGET 

c -------------- Record Speed:
      CALL SYSTEM_CLOCK(COUNT4,COUNT_RATE)
      TG34=TG34+(COUNT4-COUNT3)/REAL(COUNT_RATE)
      
c ===== pjy  
#ifdef TIDE_EL
c --- V2412:
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I=1,NUMEBC
      do concurrent (I=1:NUMEBC) local(IC,IE)
      IC = NETA(I)
      IE = NCON(I)
      IF (IC==NIM1(IE)) THEN
          RHSPOI(IE)=RHSPOI(IE)+AAUU(IE)*ELF(IC)
      ELSEIF (IC==NIP1(IE)) THEN
          RHSPOI(IE)=RHSPOI(IE)+AAUU(NIP1(IE))*ELF(IC)
      ELSEIF (IC==NJM1(IE)) THEN
          RHSPOI(IE)=RHSPOI(IE)+AAVV(IE)*ELF(IC)
      ELSE
          RHSPOI(IE)=RHSPOI(IE)+AAVV(NJP1(IE))*ELF(IC)
      ENDIF
      ENDDO
c!$OMP END TARGET 
c      APOI(ANUM)=1. !moved to after update
c      RHSPOI(N_CTRD_AG+1)=0 !moved to after update
#endif

	!call sor(APOI,IAPOI,JAPOI,RHSPOI,UPOI,ITMAX,iter1,omega)
      CALL SOLVER_SOR_TEST

c -------------- Record Speed:
      CALL SYSTEM_CLOCK(COUNT3,COUNT_RATE)
      TG35=TG35+(COUNT3-COUNT4)/REAL(COUNT_RATE)
      
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I=1,N_CTRD_AG+N_CTRD_EVG
      DO concurrent (I=1:N_CTRD_AG+N_CTRD_EVG)
          IF (FSMADD(I).GT.0.) THEN
          ELF(I)=XPOI(I)
          ENDIF
      ENDDO
c!$OMP END TARGET 
      

c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO II=1,N_CTRD_IVG
      DO concurrent (II=1:N_CTRD_IVG) local(I,ELF0,KK,K)
      I=NAG4VGC(II,1)
      IF (FSM(I).GT.0) THEN
          ELF0=0.
          DO K=1,9
              KK=NAG4VGC(II,K+1)
              ELF0=ELF0+DJ(KK)*ELF(KK)*FSM(KK)
          ENDDO
          ELF(I)=ELF0/DJ(I)
      ENDIF
      ENDDO
c!$OMP END TARGET 
***********************************************************************
***********************************************************************
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K=1,KBM1
c      DO I=1,N_CTRD
      DO concurrent (I=1:N_CTRD,K=1:KBM1) local_init(DTI)
      IF (DUM(I).GT.0.0) THEN 
      UF(I,K)=(DU(I)*UF(I,K)-GRAV*DTI*DU(I)*(ELF(I)-ELF(NIM1(I)))*
     *(H2(I)+H2(NIM1(I)))/(DJ(I)+DJ(NIM1(I))))/
     *(HU(I)+.5*(ELF(I)+ELF(NIM1(I))))
#ifdef AIRPRESSURE
      UF(I,K)=UF(I,K)-DTI*(APR(I)-APR(NIM1(I)))/(0.5*
     *(H1(I)+H1(NIM1(I)))*1.028E3)
#endif
          IF (DU(I).LT.3..AND.UF(I,K).GT.2.) THEN
              UF(I,K)=2.
          ELSEIF (DU(I).LT.3..AND.UF(I,K).LT.-2.) THEN
              UF(I,K)=-2.   !CBR: -2?
          ENDIF      
      ELSE
          UF(I,K)=0.
      ENDIF
      IF (DVM(I).GT.0.0) THEN 
      VF(I,K)=(DV(I)*VF(I,K)-GRAV*DTI*DV(I)*(ELF(I)-ELF(NJM1(I)))*
     *(H1(I)+H1(NJM1(I)))/(DJ(I)+DJ(NJM1(I))))/
     *(HV(I)+.5*(ELF(I)+ELF(NJM1(I))))  
#ifdef AIRPRESSURE
      VF(I,K)=VF(I,K)-DTI*(APR(I)-APR(NJM1(I)))/(0.5*
     *(H2(I)+H2(NJM1(I)))*1.028E3)
#endif 
          IF (DV(I).LT.3..AND.VF(I,K).GT.2.) THEN
              VF(I,K)=2.
          ELSEIF (DV(I).LT.3..AND.VF(I,K).LT.-2.) THEN
              VF(I,K)=-2.
          ENDIF 
      ELSE
          VF(I,K)=0.
      ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET 
      
c   ------ CBR Editted (Exchange Line 2 and Line 3):    
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO II=1,NNVGF4AG
c      DO K=1,9
      DO concurrent (II=1:NNVGF4AG,K=1:9) local(IC,I)
          IC=NVGF4AG(II,1)
          IF(FSM(IC).GT.0.0) THEN
              I=NVGF4AG(II,K+1)
              ELF(I)=ELF(IC)
      ENDIF 
      ENDDO
c      ENDDO
c!$OMP END TARGET 
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 1,KBM1
c      DO I = 1,N_CTRD
      DO concurrent (I=1:N_CTRD,K=1:KBM1) local(UC1,VC1)
      IF (DUM(I).NE.0.0) THEN
          UC1 = 0.5*UF(I,K)*(H2(I)+H2(NIM1(I)))
     *    -0.25*(H3(I)+H3(NIM1(I)))/(H1(I)+H1(NIM1(I)))
     *    *(VF(I,K)+VF(NIM1(I),K)+VF(NJP1(I),K)+VF(NIM1JP1(I),K)) 
          XMFLUX(I,K) = UC1*DU(I)
      ENDIF
      IF (DVM(I).NE.0.0) THEN
          VC1 = 0.5*VF(I,K)*(H1(I)+H1(NJM1(I)))
     *    -0.25*(H3(I)+H3(NJM1(I)))/(H2(I)+H2(NJM1(I)))
     *    *(UF(I,K)+UF(NIP1(I),K)+UF(NJM1(I),K)+UF(NIP1JM1(I),K))
          YMFLUX(I,K) = VC1*DV(I)
      ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET 
      
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO II=1,N_CTRD_IVG
c      DO K=1,KBM1
      DO concurrent (II=1:N_CTRD_IVG,K=1:KBM1) local(I)
      I=NAG4VGC(II,1)
      IF (DVM(I).GT.0) THEN
      YMFLUX(I,K)=YMFLUX(NAG4VGC(II,2),K)+YMFLUX(NAG4VGC(II,5),K)
     *+YMFLUX(NAG4VGC(II,8),K)
      VF(I,K)=YMFLUX(I,K)/DV(I)*2/(H1(I)+H1(NJM1(I)))
      ENDIF
      IF (DUM(I).GT.0) THEN
      XMFLUX(I,K)=XMFLUX(NAG4VGC(II,2),K)+XMFLUX(NAG4VGC(II,3),K)
     *+XMFLUX(NAG4VGC(II,4),K)
      UF(I,K)=XMFLUX(I,K)/DU(I)*2/(H2(I)+H2(NIM1(I)))
      ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET 
      
      
c!$OMP TARGET teams DEFAULTMAP(present: allocatable)
c!$OMP DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO II=1,NNIAGCW
c      DO K=1,KBM1
      DO concurrent (II=1:NNIAGCW,K=1:KBM1) local(I,IC)
      I=NIAGCW(II)
      IC=NIP1(I)
!      IF (DUM(NIP1(I)).GT.0.0) THEN 
      IF (DUM(IC).GT.0.0) THEN 
      XMFLUX(IC,K)=XMFLUX(AIJ(I,6),K)+XMFLUX(AIJ(I,7),K)
     *+XMFLUX(AIJ(I,8),K)
      UF(IC,K)=(XMFLUX(IC,K)/DU(IC)
     *+0.25*(H3(IC)+H3(NIM1(IC)))/(H1(IC)+H1(NIM1(IC)))
     **(VF(IC,K)+VF(NIM1(IC),K)+VF(NJP1(IC),K)+VF(NIM1JP1(IC),K))) 
     **2/(H2(IC)+H2(NIM1(IC)))
      ENDIF
      ENDDO
c      ENDDO
c!$OMP DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO II=1,NNIAGCE
c      DO K=1,KBM1
      DO concurrent (II=1:NNIAGCE,K=1:KBM1) local(I)
      I=NIAGCE(II)
      IF (DUM(I).GT.0.0) THEN 
      XMFLUX(I,K)=XMFLUX(AIJ(I,9),K)+XMFLUX(AIJ(I,10),K)
     *+XMFLUX(AIJ(I,11),K)
      UF(I,K)=(XMFLUX(I,K)/DU(I)
     *+0.25*(H3(I)+H3(NIM1(I)))/(H1(I)+H1(NIM1(I)))
     **(VF(I,K)+VF(NIM1(I),K)+VF(NJP1(I),K)+VF(NIM1JP1(I),K))) 
     **2/(H2(I)+H2(NIM1(I)))
      ENDIF
      ENDDO
c      ENDDO
c!$OMP DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO II=1,NNIAGCN
c      DO K=1,KBM1
      DO concurrent (II=1:NNIAGCN,K=1:KBM1) local(I)
      I=NIAGCN(II)
      IF (DVM(I).GT.0.0) THEN 
      YMFLUX(I,K)=YMFLUX(AIJ(I,9),K)+YMFLUX(AIJ(I,10),K)
     *+YMFLUX(AIJ(I,11),K)
      VF(I,K)=(YMFLUX(I,K)/DV(I)
     *+0.25*(H3(I)+H3(NJM1(I)))/(H2(I)+H2(NJM1(I)))
     **(UF(I,K)+UF(NIP1(I),K)+UF(NJM1(I),K)+UF(NIP1JM1(I),K))) 
     **2/(H1(I)+H1(NJM1(I)))
      ENDIF
      ENDDO
c      ENDDO
c!$OMP DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO II=1,NNIAGCS
c      DO K=1,KBM1
      DO concurrent (II=1:NNIAGCS,K=1:KBM1) local(I,IC)
      I=NIAGCS(II)
      IC=NJP1(I)
!      IF (DVM(NJP1(I)).GT.0.0) THEN 
      IF (DVM(IC).GT.0.0) THEN 
      YMFLUX(IC,K)=YMFLUX(AIJ(I,6),K)+YMFLUX(AIJ(I,7),K)
     *+YMFLUX(AIJ(I,8),K)
      VF(IC,K)=(YMFLUX(IC,K)/DV(IC)
     *+0.25*(H3(IC)+H3(NJM1(IC)))/(H2(IC)+H2(NJM1(IC)))
     **(UF(IC,K)+UF(NIP1(IC),K)+UF(NJM1(IC),K)+UF(NIP1JM1(IC),K))) 
     **2/(H1(IC)+H1(NJM1(IC)))
      ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET teams

c -------------- Record Speed:
      CALL SYSTEM_CLOCK(COUNT4,COUNT_RATE)
      TG36=TG36+(COUNT4-COUNT3)/REAL(COUNT_RATE)
      
      
#ifdef WEIR
      CALL ADDWEIR(2)
      CALL ADDWEIR(3)
#endif
      RETURN
      END SUBROUTINE ELTEST