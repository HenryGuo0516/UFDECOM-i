#INCLUDE "DEFS.h"
      
      SUBROUTINE DISTINGUISH
      USE MOD_GLOBAL
      IMPLICIT NONE
      
      INTEGER DIU,INUM,CNUM
      INTEGER I,II,IND,KK,JJ,K,J
      INTEGER ICE,ICW,ICS,ICN,POS,IE,IC
      REAL XB(4),YB(4)
      
      NPOI=N_CTRD_AG+N_CTRD_EVG+1
      
      ALLOCATE (AIJ(NPOI,11))         ;AIJ=0
      ALLOCATE (IAPOI(NPOI+1))        ;IAPOI=0.
      ALLOCATE (RHSPOI(NPOI))         ;RHSPOI=0.0
      ALLOCATE (XPOI(NPOI))           ;XPOI=0.0
c      ALLOCATE (DPOI(NPOI))           ;DPOI=0.0  !CBR: scalarized
      
      DO 10 I=1,NPOI-1
      IF (FSMADD(I).GT.0.0)THEN   
          AIJ(I,1)=5
          AIJ(I,2)=1
          AIJ(I,3)=2
          AIJ(I,4)=3
          AIJ(I,5)=4
          AIJ(I,6)=5
      ELSE
          AIJ(I,1)=1
          AIJ(I,2)=I
      ENDIF
10    CONTINUE
      AIJ(NPOI,1)=1
      AIJ(NPOI,2)=NPOI
***********************************************************************
*                                                                     *
*                             AG BOUNDARY                             *
*                                                                     *
***********************************************************************
      DO 11 I=1,NNIAGCW
      INUM=NIAGCW(I)
      IF (FSMADD(INUM).EQ.0) GO TO 11
      AIJ(INUM,1)=10
      AIJ(INUM,2)=NJM1(INUM)
      AIJ(INUM,3)=NIM1(INUM)
      AIJ(INUM,4)=INUM
      AIJ(INUM,5)=NJP1(INUM)
      DIU=0
      DO II=1,4
          XB(II)=XNODE(II,NIP1(INUM))
          YB(II)=YNODE(II,NIP1(INUM))
      ENDDO
      DO II=1,NNIAGFW
          CALL INSIDE(XR(NIAGFW(II)),YR(NIAGFW(II)),XB,YB,NB,IND)
          IF(IND==1) THEN
              DIU=DIU+1
              AIJ(INUM,5+DIU)=NIAGFW(II)
              AIJ(INUM,8+DIU)=NIM1(NIAGFW(II))
          ENDIF
      ENDDO
11    CONTINUE
      DO 12 I=1,NNIAGCE
      INUM=NIAGCE(I)
      IF (FSMADD(INUM).EQ.0) GO TO 12
      AIJ(INUM,1)=10
      AIJ(INUM,2)=NJM1(INUM)
      AIJ(INUM,3)=INUM
      AIJ(INUM,4)=NIP1(INUM)
      AIJ(INUM,5)=NJP1(INUM)
      DIU=0
      DO II=1,4
          XB(II)=XNODE(II,NIM1(INUM))
          YB(II)=YNODE(II,NIM1(INUM))
      ENDDO
      DO II=1,NNIAGFE
          CALL INSIDE(XR(NIAGFE(II)),YR(NIAGFE(II)),XB,YB,NB,IND)
          IF(IND==1) THEN
              DIU=DIU+1
              AIJ(INUM,5+DIU)=NIAGFE(II)
              AIJ(INUM,8+DIU)=NIP1(NIAGFE(II))
          ENDIF
      ENDDO
12    CONTINUE     
      DO 13 I=1,NNIAGCN
      INUM=NIAGCN(I)
      IF (FSMADD(INUM).EQ.0) GO TO 13
      AIJ(INUM,1)=10
      AIJ(INUM,2)=NIM1(INUM)
      AIJ(INUM,3)=INUM
      AIJ(INUM,4)=NIP1(INUM)
      AIJ(INUM,5)=NJP1(INUM)
      DIU=0
      DO II=1,4
          XB(II)=XNODE(II,NJM1(INUM))
          YB(II)=YNODE(II,NJM1(INUM))
      ENDDO
      DO II=1,NNIAGFN
          CALL INSIDE(XR(NIAGFN(II)),YR(NIAGFN(II)),XB,YB,NB,IND)
          IF(IND==1) THEN
              DIU=DIU+1
              AIJ(INUM,5+DIU)=NIAGFN(II)
              AIJ(INUM,8+DIU)=NJP1(NIAGFN(II))
          ENDIF
      ENDDO
13    CONTINUE 
      DO 14 I=1,NNIAGCS
      INUM=NIAGCS(I)
      IF (FSMADD(INUM).EQ.0) GO TO 14
      AIJ(INUM,1)=10
      AIJ(INUM,2)=NJM1(INUM)
      AIJ(INUM,3)=NIM1(INUM)
      AIJ(INUM,4)=INUM
      AIJ(INUM,5)=NIP1(INUM)
      DIU=0
      DO II=1,4
          XB(II)=XNODE(II,NJP1(INUM))
          YB(II)=YNODE(II,NJP1(INUM))
      ENDDO
      DO II=1,NNIAGFS
          CALL INSIDE(XR(NIAGFS(II)),YR(NIAGFS(II)),XB,YB,NB,IND)
          IF(IND==1) THEN
              DIU=DIU+1
              AIJ(INUM,5+DIU)=NIAGFS(II)
              AIJ(INUM,8+DIU)=NJM1(NIAGFS(II))
          ENDIF
      ENDDO
14    CONTINUE
      DO II=1,NNIAGFS
      I=NIAGFS(II)
      IF(FSMADD(I).GT.0.) THEN
      AIJ(I,2)=5
      AIJ(I,3)=AIJ(I,3)-1
      AIJ(I,4)=AIJ(I,4)-1
      AIJ(I,5)=AIJ(I,5)-1
      AIJ(I,6)=AIJ(I,6)-1
      ENDIF
      I=NJM1(I)
      IF(FSMADD(I).GT.0.) THEN
      AIJ(I,2)=2
      AIJ(I,3)=3
      AIJ(I,4)=4
      AIJ(I,5)=5
      AIJ(I,6)=1
      ENDIF
      ENDDO
      DO II=1,NNIAGFW
      I=NIAGFW(II)
      IF(FSMADD(I).GT.0.) THEN
      IF(AIJ(I,2).EQ.5) AIJ(I,2)=AIJ(I,2)-1
      AIJ(I,3)=5
      AIJ(I,4)=AIJ(I,4)-1
      AIJ(I,5)=AIJ(I,5)-1
      AIJ(I,6)=AIJ(I,6)-1
      ENDIF
      I=NIM1(I)
      IF(FSMADD(I).GT.0.) THEN
      AIJ(I,2)=2
      AIJ(I,3)=3
      AIJ(I,4)=4
      AIJ(I,5)=1
      AIJ(I,6)=5
      ENDIF
      ENDDO
      DO II=1,NNIAGFE
      I=NIAGFE(II)
      IF(FSMADD(I).GT.0.) THEN
      IF(AIJ(I,2).EQ.5) AIJ(I,2)=AIJ(I,2)-1
      IF(AIJ(I,3).EQ.5) AIJ(I,3)=AIJ(I,3)-1
      IF(AIJ(I,4).EQ.5) AIJ(I,4)=AIJ(I,4)-1
      AIJ(I,5)=5
      AIJ(I,6)=AIJ(I,6)-1
      ENDIF
      I=NIP1(I)
      IF(FSMADD(I).GT.0.) THEN
      AIJ(I,2)=2
      AIJ(I,3)=1
      AIJ(I,4)=3
      AIJ(I,5)=4
      AIJ(I,6)=5
      ENDIF
      ENDDO
      DO II=1,NNIAGFN
      I=NIAGFN(II)
      IF(FSMADD(I).GT.0.) THEN
      IF(AIJ(I,2).EQ.5) AIJ(I,2)=AIJ(I,2)-1
      IF(AIJ(I,3).EQ.5) AIJ(I,3)=AIJ(I,3)-1
      IF(AIJ(I,4).EQ.5) AIJ(I,4)=AIJ(I,4)-1
      IF(AIJ(I,5).EQ.5) AIJ(I,5)=AIJ(I,5)-1
      AIJ(I,6)=5
      ENDIF
      I=NJP1(I)
      IF(FSMADD(I).GT.0.) THEN
      AIJ(I,2)=1
      AIJ(I,3)=2
      AIJ(I,4)=3
      AIJ(I,5)=4
      AIJ(I,6)=5
      ENDIF
      ENDDO
***********************************************************************
*                                                                     *
*                             VG BOUNDARY                             *
*                                                                     *
***********************************************************************
      KK=1
      DO 21 I=1,NNIVGFW
      INUM=NIVGFW(I)
      IF (FSMADD(INUM).EQ.0) GO TO 21
      AIJ(INUM,1)=6
      AIJ(INUM,4)=NJM1(INUM)
      AIJ(INUM,5)=INUM
      AIJ(INUM,6)=NIP1(INUM)
      AIJ(INUM,7)=NJP1(INUM)
      DO II=1,4
          XB(II)=XNODE(II,KK)
          YB(II)=YNODE(II,KK)
      ENDDO
      CALL INSIDE(XR(INUM),YR(INUM),XB,YB,NB,IND)
      IF(IND==1) THEN
          AIJ(INUM,2)=NIM1(KK)
          AIJ(INUM,3)=KK
          DO JJ=1,NNIVGFWALL
              IF (NIVGFWALL(JJ,1).EQ.INUM) NIVGFWALL(JJ,2)=KK
          ENDDO
          GO TO 21
      ENDIF
      KK=NJP1(KK)
      DO II=1,4
          XB(II)=XNODE(II,KK)
          YB(II)=YNODE(II,KK)
      ENDDO
      CALL  INSIDE(XR(INUM),YR(INUM),XB,YB,NB,IND)
      IF(IND==1) THEN
          AIJ(INUM,2)=NIM1(KK)
          AIJ(INUM,3)=KK
          DO JJ=1,NNIVGFWALL
              IF (NIVGFWALL(JJ,1).EQ.INUM) NIVGFWALL(JJ,2)=KK
          ENDDO
          GO TO 21
      ENDIF
      DO 211  KK=1,N_CTRD_AG
          DO II=1,4
              XB(II)=XNODE(II,KK)
              YB(II)=YNODE(II,KK)
          ENDDO
          CALL  INSIDE(XR(INUM),YR(INUM),XB,YB,NB,IND)
          IF(IND==1) THEN
              AIJ(INUM,2)=NIM1(KK)
              AIJ(INUM,3)=KK
              DO JJ=1,NNIVGFWALL
                  IF (NIVGFWALL(JJ,1).EQ.INUM) NIVGFWALL(JJ,2)=KK
              ENDDO
              GO TO 21
          ENDIF
211   CONTINUE
      PRINT*,'FUCCK',INUM
      PAUSE
21    CONTINUE
      KK=1
      DO 22 I=1,NNIVGFE
      INUM=NIVGFE(I)
      IF (FSMADD(INUM).EQ.0) GO TO 22
      AIJ(INUM,1)=6
      AIJ(INUM,4)=NJM1(INUM)
      AIJ(INUM,5)=NIM1(INUM)
      AIJ(INUM,6)=INUM
      AIJ(INUM,7)=NJP1(INUM)
      DO II=1,4
          XB(II)=XNODE(II,KK)
          YB(II)=YNODE(II,KK)
      ENDDO
      CALL INSIDE(XR(INUM),YR(INUM),XB,YB,NB,IND)
      IF(IND==1) THEN
          AIJ(INUM,2)=KK
          AIJ(INUM,3)=NIP1(KK)
          DO JJ=1,NNIVGFEALL
              IF (NIVGFEALL(JJ,1).EQ.INUM) NIVGFEALL(JJ,2)=KK
          ENDDO
          GO TO 22
      ENDIF
      KK=NJP1(KK)
      DO II=1,4
          XB(II)=XNODE(II,KK)
          YB(II)=YNODE(II,KK)
      ENDDO
      CALL  INSIDE(XR(INUM),YR(INUM),XB,YB,NB,IND)
      IF(IND==1) THEN
          AIJ(INUM,2)=KK
          AIJ(INUM,3)=NIP1(KK)
          DO JJ=1,NNIVGFEALL
              IF (NIVGFEALL(JJ,1).EQ.INUM) NIVGFEALL(JJ,2)=KK
          ENDDO
          GO TO 22
      ENDIF
      DO 221  KK=1,N_CTRD_AG
          DO II=1,4
              XB(II)=XNODE(II,KK)
              YB(II)=YNODE(II,KK)
          ENDDO
          CALL  INSIDE(XR(INUM),YR(INUM),XB,YB,NB,IND)
          IF(IND==1) THEN
              AIJ(INUM,2)=KK
              AIJ(INUM,3)=NIP1(KK)
              DO JJ=1,NNIVGFEALL
                  IF (NIVGFEALL(JJ,1).EQ.INUM) NIVGFEALL(JJ,2)=KK
              ENDDO
              GO TO 22
          ENDIF
221   CONTINUE
      PRINT*,'FUCCK',INUM
      PAUSE
22    CONTINUE
      KK=1
      DO 23 I=1,NNIVGFS
      INUM=NIVGFS(I)
      IF (FSMADD(INUM).EQ.0) GO TO 23
      AIJ(INUM,1)=6
      AIJ(INUM,4)=NIM1(INUM)
      AIJ(INUM,5)=INUM
      AIJ(INUM,6)=NIP1(INUM)
      AIJ(INUM,7)=NJP1(INUM)
      DO II=1,4
          XB(II)=XNODE(II,KK)
          YB(II)=YNODE(II,KK)
      ENDDO
      CALL INSIDE(XR(INUM),YR(INUM),XB,YB,NB,IND)
      IF(IND==1) THEN
          AIJ(INUM,2)=NJM1(KK)
          AIJ(INUM,3)=KK
          DO JJ=1,NNIVGFSALL
              IF (NIVGFSALL(JJ,1).EQ.INUM) NIVGFSALL(JJ,2)=KK
          ENDDO
          GO TO 23
      ENDIF
      KK=NIP1(KK)
      DO II=1,4
          XB(II)=XNODE(II,KK)
          YB(II)=YNODE(II,KK)
      ENDDO
      CALL  INSIDE(XR(INUM),YR(INUM),XB,YB,NB,IND)
      IF(IND==1) THEN
          AIJ(INUM,2)=NJM1(KK)
          AIJ(INUM,3)=KK
          DO JJ=1,NNIVGFSALL
              IF (NIVGFSALL(JJ,1).EQ.INUM) NIVGFSALL(JJ,2)=KK
          ENDDO
          GO TO 23
      ENDIF
      DO 231  KK=1,N_CTRD_AG
          DO II=1,4
              XB(II)=XNODE(II,KK)
              YB(II)=YNODE(II,KK)
          ENDDO
          CALL  INSIDE(XR(INUM),YR(INUM),XB,YB,NB,IND)
          IF(IND==1) THEN
              AIJ(INUM,2)=NJM1(KK)
              AIJ(INUM,3)=KK
              DO JJ=1,NNIVGFSALL
                  IF (NIVGFSALL(JJ,1).EQ.INUM) NIVGFSALL(JJ,2)=KK
              ENDDO
              GO TO 23
          ENDIF
231   CONTINUE
      PRINT*,'FUCCK',INUM
      PAUSE
23    CONTINUE
      KK=1
      DO 24 I=1,NNIVGFN
      INUM=NIVGFN(I)
      IF (FSMADD(INUM).EQ.0) GO TO 24
      AIJ(INUM,1)=6
      AIJ(INUM,4)=NJM1(INUM)
      AIJ(INUM,5)=NIM1(INUM)
      AIJ(INUM,6)=INUM
      AIJ(INUM,7)=NIP1(INUM)
      DO II=1,4
          XB(II)=XNODE(II,KK)
          YB(II)=YNODE(II,KK)
      ENDDO
      CALL INSIDE(XR(INUM),YR(INUM),XB,YB,NB,IND)
      IF(IND==1) THEN
          AIJ(INUM,2)=KK
          AIJ(INUM,3)=NJP1(KK)
          DO JJ=1,NNIVGFNALL
              IF (NIVGFNALL(JJ,1).EQ.INUM) NIVGFNALL(JJ,2)=KK
          ENDDO
          GO TO 24
      ENDIF
      KK=NIP1(KK)
      DO II=1,4
          XB(II)=XNODE(II,KK)
          YB(II)=YNODE(II,KK)
      ENDDO
      CALL  INSIDE(XR(INUM),YR(INUM),XB,YB,NB,IND)
      IF(IND==1) THEN
          AIJ(INUM,2)=KK
          AIJ(INUM,3)=NJP1(KK)
          DO JJ=1,NNIVGFNALL
              IF (NIVGFNALL(JJ,1).EQ.INUM) NIVGFNALL(JJ,2)=KK
          ENDDO
          GO TO 24
      ENDIF
      DO 241  KK=1,N_CTRD_AG
          DO II=1,4
              XB(II)=XNODE(II,KK)
              YB(II)=YNODE(II,KK)
          ENDDO
          CALL  INSIDE(XR(INUM),YR(INUM),XB,YB,NB,IND)
          IF(IND==1) THEN
              AIJ(INUM,2)=KK
              AIJ(INUM,3)=NJP1(KK)
              DO JJ=1,NNIVGFNALL
                  IF (NIVGFNALL(JJ,1).EQ.INUM) NIVGFNALL(JJ,2)=KK
              ENDDO
              GO TO 24
          ENDIF
241   CONTINUE
      PRINT*,'FUCCK',INUM
      PAUSE
24    CONTINUE
      DO 25 I=1,NNIVGFWS
      INUM=NIVGFWS(I)
      IF (FSMADD(INUM).EQ.0) GO TO 25
      AIJ(INUM,1)=6
      AIJ(INUM,5)=INUM
      AIJ(INUM,6)=NIP1(INUM)
      AIJ(INUM,7)=NJP1(INUM)
      DO 251  KK=1,N_CTRD_AG
          DO II=1,4
              XB(II)=XNODE(II,KK)
              YB(II)=YNODE(II,KK)
          ENDDO
          CALL  INSIDE(XR(INUM),YR(INUM),XB,YB,NB,IND)
          IF(IND==1) THEN
              AIJ(INUM,2)=NJM1(KK)
              AIJ(INUM,3)=NIM1(KK)
              AIJ(INUM,4)=KK
              DO JJ=1,NNIVGFWALL
                  IF (NIVGFWALL(JJ,1).EQ.INUM) NIVGFWALL(JJ,2)=KK
              ENDDO
              DO JJ=1,NNIVGFSALL
                  IF (NIVGFSALL(JJ,1).EQ.INUM) NIVGFSALL(JJ,2)=KK
              ENDDO
              GO TO 25
          ENDIF
251   CONTINUE
      PRINT*,'FUCCK',INUM
      PAUSE
25    CONTINUE
      DO 26 I=1,NNIVGFWN
      INUM=NIVGFWN(I)
      IF (FSMADD(INUM).EQ.0) GO TO 26
      AIJ(INUM,1)=6
      AIJ(INUM,5)=NJM1(INUM)
      AIJ(INUM,6)=INUM
      AIJ(INUM,7)=NIP1(INUM)
      DO 261  KK=1,N_CTRD_AG
          DO II=1,4
              XB(II)=XNODE(II,KK)
              YB(II)=YNODE(II,KK)
          ENDDO
          CALL  INSIDE(XR(INUM),YR(INUM),XB,YB,NB,IND)
          IF(IND==1) THEN
              AIJ(INUM,2)=NIM1(KK)
              AIJ(INUM,3)=KK
              AIJ(INUM,4)=NJP1(KK)
              DO JJ=1,NNIVGFWALL
                  IF (NIVGFWALL(JJ,1).EQ.INUM) NIVGFWALL(JJ,2)=KK
              ENDDO
              DO JJ=1,NNIVGFNALL
                  IF (NIVGFNALL(JJ,1).EQ.INUM) NIVGFNALL(JJ,2)=KK
              ENDDO
              GO TO 26
          ENDIF
261   CONTINUE
      PRINT*,'FUCCK',INUM
      PAUSE
26    CONTINUE
      DO 27 I=1,NNIVGFES
      INUM=NIVGFES(I)
      IF (FSMADD(INUM).EQ.0) GO TO 27
      AIJ(INUM,1)=6
      AIJ(INUM,5)=NIM1(INUM)
      AIJ(INUM,6)=INUM
      AIJ(INUM,7)=NJP1(INUM)
      DO 271  KK=1,N_CTRD_AG
          DO II=1,4
              XB(II)=XNODE(II,KK)
              YB(II)=YNODE(II,KK)
          ENDDO
          CALL  INSIDE(XR(INUM),YR(INUM),XB,YB,NB,IND)
          IF(IND==1) THEN
              AIJ(INUM,2)=NJM1(KK)
              AIJ(INUM,3)=KK
              AIJ(INUM,4)=NIP1(KK)
              DO JJ=1,NNIVGFEALL
                  IF (NIVGFEALL(JJ,1).EQ.INUM) NIVGFEALL(JJ,2)=KK
              ENDDO
              DO JJ=1,NNIVGFSALL
                  IF (NIVGFSALL(JJ,1).EQ.INUM) NIVGFSALL(JJ,2)=KK
              ENDDO
              GO TO 27
          ENDIF
271   CONTINUE
      PRINT*,'FUCCK',INUM
      PAUSE
27    CONTINUE
      DO 28 I=1,NNIVGFEN
      INUM=NIVGFEN(I)
      IF (FSMADD(INUM).EQ.0) GO TO 28
      AIJ(INUM,1)=6
      AIJ(INUM,5)=NJM1(INUM)
      AIJ(INUM,6)=NIM1(INUM)
      AIJ(INUM,7)=INUM
      DO 281  KK=1,N_CTRD_AG
          DO II=1,4
              XB(II)=XNODE(II,KK)
              YB(II)=YNODE(II,KK)
          ENDDO
          CALL  INSIDE(XR(INUM),YR(INUM),XB,YB,NB,IND)
          IF(IND==1) THEN
              AIJ(INUM,2)=KK
              AIJ(INUM,3)=NIP1(KK)
              AIJ(INUM,4)=NJP1(KK)
              DO JJ=1,NNIVGFEALL
                  IF (NIVGFEALL(JJ,1).EQ.INUM) NIVGFEALL(JJ,2)=KK
              ENDDO
              DO JJ=1,NNIVGFNALL
                  IF (NIVGFNALL(JJ,1).EQ.INUM) NIVGFNALL(JJ,2)=KK
              ENDDO
              GO TO 28
          ENDIF
281   CONTINUE
      PRINT*,'FUCCK',INUM
      PAUSE
28    CONTINUE
      ANUM=0
      IAPOI(1)=1
      DO I=1,NPOI
          ANUM=ANUM+AIJ(I,1)
          IAPOI(I+1)=ANUM+1
      ENDDO
      ALLOCATE(JAPOI(ANUM))   ;JAPOI=0
      ALLOCATE(JAPOI2(ANUM))  ;JAPOI2=0
      ALLOCATE(APOI(ANUM))    ;APOI=0. 
      ANUM=0
      DO I=1,NPOI
      IF (AIJ(I,1).NE.5) THEN
          DO II=1,AIJ(I,1)
              JAPOI(ANUM+II)=AIJ(I,II+1)
              IF(JAPOI(ANUM+II).GE.NPOI) JAPOI(ANUM+II)=NPOI
          ENDDO
          ANUM=ANUM+AIJ(I,1)
      ELSE
          JAPOI(ANUM+AIJ(I,2))=NJM1(I)
          JAPOI(ANUM+AIJ(I,3))=NIM1(I)
          JAPOI(ANUM+AIJ(I,4))=I
          JAPOI(ANUM+AIJ(I,5))=NIP1(I)
          JAPOI(ANUM+AIJ(I,6))=NJP1(I)
          DO II=1,5
              IF (JAPOI(ANUM+II).GE.NPOI) JAPOI(ANUM+II)=NPOI
          ENDDO
          ANUM=ANUM+5
      ENDIF
      ENDDO
c --- CBR: The only application for JAPOI(K) is via KAPOI(I)
c --- No need to calculate KAPOI(I) every step:
      ALLOCATE(KAPOI(NPOI))    ;KAPOI=0. 
      DO I=1,NPOI
c          IDX_S=IAPOI(I)
c          IDX_E=IAPOI(I+1)-1
          DO K=IAPOI(I),IAPOI(I+1)-1
              IF (JAPOI(K).EQ.I) THEN !don't have to for every step
                  KAPOI(I)=K
c                  DPOI(I)=1./APOI(K)
                  EXIT
              ENDIF
          ENDDO
      ENDDO
      
      DO II=1,N_CTRD_IVG
      I=N_CTRD_AG+N_CTRD_EVG+II
      IF (FSM(I).GT.0) THEN
      H(I)=0.
      DO K=1,9
          CNUM=NAG4VGC(II,K)
          H(I)=H(I)+DJ(CNUM)*H(CNUM)*FSM(CNUM)
      ENDDO
      H(I)=H(I)/DJ(I)
      ENDIF
      ENDDO
      DO II=1,NNVGF4AG
      IC=NVGF4AG(II,1)
      IF(FSM(IC).GT.0.0) THEN
          DO K=1,9
              I=NVGF4AG(II,K+1)
              H(I)=H(IC)
          ENDDO
      ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DISTINGUISH