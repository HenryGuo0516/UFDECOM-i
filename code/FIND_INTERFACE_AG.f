#INCLUDE "DEFS.h"
      
      SUBROUTINE  FIND_INTERFACE_AG 
!=======================================================================
! THIS SUBROUTINE FINDS THE ACTUAL GRIDS AT THE NESTED INTERFACE.
!=======================================================================
      USE MOD_GLOBAL
      IMPLICIT NONE
      
      INTEGER I, J, K, II, JJ ,KK,IND,IC
      INTEGER ICW, ICE, ICN, ICS
      INTEGER IFW, IFE, IFN, IFS
      INTEGER VFW, VFE, VFN, VFS
      INTEGER VFWALL, VFEALL, VFNALL, VFSALL
      INTEGER VFWS, VFES, VFWN, VFEN 
      REAL XB(4),YB(4)
      
      NNVGF4AG=N_CTRD_EVG/9
***********************************************************************
*                                                                     *
*                                                                     *
*                                                                     *
***********************************************************************      
      ALLOCATE (NVGF4AG(NNVGF4AG,10))
      ALLOCATE (NAG4VGC(N_CTRD_IVG,10))
      NVGF4AG=0.
      NAG4VGC=0.
      JJ=0
      DO 101 II=1,N_CTRD_EVG
      I=N_CTRD_AG+II
      DO KK=1,N_CTRD_AG
          DO K=1,4
              XB(K)=XNODE(K,KK)
              YB(K)=YNODE(K,KK)
          ENDDO
          CALL INSIDE(XR(I),YR(I),XB,YB,NB,IND)
          IF (IND.EQ.1) THEN
              DO J=1,NNVGF4AG
                  IF (KK.EQ.NVGF4AG(J,1)) GOTO 101
              ENDDO
              IF (FSM(KK).EQ.0) GOTO 101
              JJ=JJ+1
              NVGF4AG(JJ,1)=KK
              
              IC=I
              NVGF4AG(JJ,2)=IC
              NVGF4AG(JJ,3)=NIP1(IC)
              NVGF4AG(JJ,4)=NIP2(IC)
              IC=NJP1(IC)
              NVGF4AG(JJ,5)=IC
              NVGF4AG(JJ,6)=NIP1(IC)
              NVGF4AG(JJ,7)=NIP2(IC)
              IC=NJP1(IC)
              NVGF4AG(JJ,8)=IC
              NVGF4AG(JJ,9)=NIP1(IC)
              NVGF4AG(JJ,10)=NIP2(IC)
          ENDIF
      ENDDO
101   CONTINUE 
      NNVGF4AG=JJ
      DO 102 II=1,N_CTRD_IVG
          I=N_CTRD_AG+N_CTRD_EVG+II
          DO K=1,4
              XB(K)=XNODE(K,I)
              YB(K)=YNODE(K,I)
          ENDDO
          DO KK=1,N_CTRD_AG
              CALL INSIDE(XR(KK),YR(KK),XB,YB,NB,IND)
              IF (IND.EQ.1) THEN
                  NAG4VGC(II,1)=I
                  IC=KK
                  NAG4VGC(II,2)=IC
                  NAG4VGC(II,3)=NIP1(IC)
                  NAG4VGC(II,4)=NIP2(IC)
                  IC=NJP1(IC)
                  NAG4VGC(II,5)=IC
                  NAG4VGC(II,6)=NIP1(IC)
                  NAG4VGC(II,7)=NIP2(IC)
                  IC=NJP1(IC)
                  NAG4VGC(II,8)=IC
                  NAG4VGC(II,9)=NIP1(IC)
                  NAG4VGC(II,10)=NIP2(IC)
                  GO TO 102
              ENDIF
          ENDDO
102   CONTINUE
      NNIVGW=0.
      NNIVGS=0
      DO  I = N_CTRD_AG+N_CTRD_EVG+1,N_CTRD
          IF (NIP1(I).GT.N_CTRD) THEN
              NNIVGW=NNIVGW+1
          ELSEIF (NJP1(I).GT.N_CTRD) THEN
              NNIVGS=NNIVGS+1
          ENDIF
      ENDDO
      ALLOCATE (NIVGW(NNIVGW,4))
      ALLOCATE (NIVGS(NNIVGS,4))
      ICW=0
      ICS=0
      DO 103 I = N_CTRD_AG+N_CTRD_EVG+1,N_CTRD
          IF (NIP1(I).GT.N_CTRD) THEN
              ICW=ICW+1
              NIVGW(ICW,1)=I
              DO K=1,4
                  XB(K)=XNODE(K,I)
                  YB(K)=YNODE(K,I)
              ENDDO
              DO KK=1,N_CTRD_AG
                  CALL INSIDE(XR(KK),YR(KK),XB,YB,NB,IND)
                  IF (IND.EQ.1) THEN
                      NIVGW(ICW,2)=NIP3(KK);
                      NIVGW(ICW,3)=NJP1(NIP3(KK));
                      NIVGW(ICW,4)=NJP2(NIP3(KK));
                      GO TO 103
                  ENDIF
              ENDDO
          ELSEIF (NJP1(I).GT.N_CTRD) THEN
              ICS=ICS+1
              NIVGS(ICS,1)=I
              DO K=1,4
                  XB(K)=XNODE(K,I)
                  YB(K)=YNODE(K,I)
              ENDDO
              DO KK=1,N_CTRD_AG
                  CALL INSIDE(XR(KK),YR(KK),XB,YB,NB,IND)
                  IF (IND.EQ.1) THEN
                      NIVGS(ICS,2)=NJP3(KK);
                      NIVGS(ICS,3)=NIP1(NJP3(KK));
                      NIVGS(ICS,4)=NIP2(NJP3(KK));
                      GO TO 103
                  ENDIF
              ENDDO
          ENDIF
103   CONTINUE     
               
***********************************************************************
*                                                                     *
*                             AG BOUNDARY                             *
*                                                                     *
***********************************************************************
      NNIAGCW = 0
      NNIAGCE = 0
      NNIAGCN = 0
      NNIAGCS = 0
      NNIAGFW = 0
      NNIAGFE = 0
      NNIAGFN = 0
      NNIAGFS = 0
      DO 10 I = 1,N_CTRD_AG
      IF (NIP1(I).GT.N_CTRD_AG+N_CTRD_EVG .AND. NIP1(I).LE.N_CTRD)  
     * NNIAGCW = NNIAGCW+1
      IF (NIM1(I).GT.N_CTRD_AG .AND. NIM1(I).LE.N_CTRD_AG+N_CTRD_EVG) 
     * NNIAGFW = NNIAGFW+1
      IF(NIM1(I).GT.N_CTRD_AG+N_CTRD_EVG .AND. NIM1(I).LE.N_CTRD)
     * NNIAGCE = NNIAGCE+1
      IF(NIP1(I).GT.N_CTRD_AG .AND. NIP1(I).LE.N_CTRD_AG+N_CTRD_EVG)
     * NNIAGFE = NNIAGFE+1   
      IF (NJM1(I).GT.N_CTRD_AG+N_CTRD_EVG .AND. NJM1(I).LE.N_CTRD) 
     *NNIAGCN = NNIAGCN+1  
      IF(NJP1(I).GT.N_CTRD_AG .AND. NJP1(I).LE.N_CTRD_AG+N_CTRD_EVG)
     *NNIAGFN = NNIAGFN+1
      IF (NJP1(I).GT.N_CTRD_AG+N_CTRD_EVG .AND. NJP1(I).LE.N_CTRD)
     *NNIAGCS = NNIAGCS+1    
      IF (NJM1(I).GT.N_CTRD_AG .AND. NJM1(I).LE.N_CTRD_AG+N_CTRD_EVG)
     *NNIAGFS = NNIAGFS+1 
10    CONTINUE
      ALLOCATE (NIAGCW(NNIAGCW),NIAGCE(NNIAGCE))
      ALLOCATE (NIAGCN(NNIAGCN),NIAGCS(NNIAGCS))
      ALLOCATE (NIAGFW(NNIAGFW),NIAGFE(NNIAGFE))
      ALLOCATE (NIAGFN(NNIAGFN),NIAGFS(NNIAGFS))
      ICW = 0
      ICE = 0
      ICN = 0
      ICS = 0
      IFW = 0
      IFE = 0
      IFN = 0
      IFS = 0
      DO 20 I = 1,N_CTRD_AG
      IF (NJM1(I).GT.N_CTRD_AG+N_CTRD_EVG .AND. NJM1(I).LE.N_CTRD)
     *THEN
          ICN = ICN+1
          NIAGCN(ICN) = I 
      ENDIF
      IF(NJP1(I).GT.N_CTRD_AG .AND. NJP1(I).LE.N_CTRD_AG+N_CTRD_EVG)
     *THEN! NORTH INTERFACE FINE GRID
          IFN = IFN+1
          NIAGFN(IFN) = I
      ENDIF
      IF(NIM1(I).GT.N_CTRD_AG+N_CTRD_EVG .AND. NIM1(I).LE.N_CTRD)
     *THEN
          ICE = ICE+1
          NIAGCE(ICE) = I
      ENDIF
      IF(NIP1(I).GT.N_CTRD_AG .AND. NIP1(I).LE.N_CTRD_AG+N_CTRD_EVG)
     *THEN
          IFE = IFE+1
          NIAGFE(IFE) = I
      ENDIF
      IF (NIP1(I).GT.N_CTRD_AG+N_CTRD_EVG .AND. NIP1(I).LE.N_CTRD) 
     *THEN
          ICW = ICW+1
          NIAGCW(ICW) = I
      ENDIF
      IF (NIM1(I).GT.N_CTRD_AG .AND. NIM1(I).LE.N_CTRD_AG+N_CTRD_EVG)
     *THEN 
          IFW = IFW+1
          NIAGFW(IFW)=I
      ENDIF
      IF (NJP1(I).GT.N_CTRD_AG+N_CTRD_EVG .AND. NJP1(I).LE.N_CTRD)
     *THEN
          ICS = ICS+1
          NIAGCS(ICS) = I  
      ENDIF
      IF(NJM1(I).GT.N_CTRD_AG .AND. NJM1(I).LE.N_CTRD_AG+N_CTRD_EVG) 
     *THEN
          IFS = IFS+1
          NIAGFS(IFS) = I
      ENDIF 
20    CONTINUE
***********************************************************************
*                                                                     *
*                             VG BOUNDARY                             *
*                                                                     *
***********************************************************************      
      DO 30 I=N_CTRD_AG,N_CTRD_AG+N_CTRD_EVG
      IF (NIM1(I).LT.N_CTRDP1.AND.NIP1(I).LT.N_CTRDP1.AND.
     *NJM1(I).LT.N_CTRDP1.AND.NJP1(I).LT.N_CTRDP1) THEN
          GO TO 30
      ELSEIF (NIM1(I).EQ.N_CTRDP1.AND.NIP1(I).LT.N_CTRDP1.AND.
     *NJM1(I).LT.N_CTRDP1.AND.NJP1(I).LT.N_CTRDP1) THEN
          NNIVGFW=NNIVGFW+1
      ELSEIF (NIM1(I).LT.N_CTRDP1.AND.NIP1(I).EQ.N_CTRDP1.AND.
     *NJM1(I).LT.N_CTRDP1.AND.NJP1(I).LT.N_CTRDP1) THEN
          NNIVGFE=NNIVGFE+1
      ELSEIF (NIM1(I).LT.N_CTRDP1.AND.NIP1(I).LT.N_CTRDP1.AND.
     *NJM1(I).EQ.N_CTRDP1.AND.NJP1(I).LT.N_CTRDP1) THEN
          NNIVGFS=NNIVGFS+1
      ELSEIF (NIM1(I).LT.N_CTRDP1.AND.NIP1(I).LT.N_CTRDP1.AND.
     *NJM1(I).LT.N_CTRDP1.AND.NJP1(I).EQ.N_CTRDP1) THEN
          NNIVGFN=NNIVGFN+1
      ELSEIF (NIM1(I).EQ.N_CTRDP1.AND.NIP1(I).LT.N_CTRDP1.AND.
     *NJM1(I).EQ.N_CTRDP1.AND.NJP1(I).LT.N_CTRDP1) THEN
          NNIVGFWS=NNIVGFWS+1
      ELSEIF (NIM1(I).EQ.N_CTRDP1.AND.NIP1(I).LT.N_CTRDP1.AND.
     *NJM1(I).LT.N_CTRDP1.AND.NJP1(I).EQ.N_CTRDP1) THEN
          NNIVGFWN=NNIVGFWN+1
      ELSEIF (NIM1(I).LT.N_CTRDP1.AND.NIP1(I).EQ.N_CTRDP1.AND.
     *NJM1(I).EQ.N_CTRDP1.AND.NJP1(I).LT.N_CTRDP1) THEN
          NNIVGFES=NNIVGFES+1
      ELSEIF (NIM1(I).LT.N_CTRDP1.AND.NIP1(I).EQ.N_CTRDP1.AND.
     *NJM1(I).LT.N_CTRDP1.AND.NJP1(I).EQ.N_CTRDP1) THEN
          NNIVGFEN=NNIVGFEN+1
      ELSE
          PRINT*,'CHECK GRID ',I
          PAUSE
      ENDIF
30    CONTINUE
      NNIVGCW=NNIVGFW+NNIVGFWS+NNIVGFWN
      NNIVGCE=NNIVGFE+NNIVGFES+NNIVGFEN
      NNIVGCN=NNIVGFN+NNIVGFWN+NNIVGFEN
      NNIVGCS=NNIVGFS+NNIVGFWS+NNIVGFES
      IF (MOD(NNIVGCW,3).NE.0) THEN
          PRINT*,'WRONG IN WEST VG BOUNDARY'
          PAUSE
          STOP
      ELSE
          NNIVGFWALL=NNIVGCW
          NNIVGCW=NNIVGCW/3
      ENDIF
      IF (MOD(NNIVGCE,3).NE.0) THEN
          PRINT*,'WRONG IN EASST VG BOUNDARY'
          PAUSE
          STOP
      ELSE
          NNIVGFEALL=NNIVGCE
          NNIVGCE=NNIVGCE/3
      ENDIF
      IF (MOD(NNIVGCN,3).NE.0) THEN
          PRINT*,'WRONG IN NORTH VG BOUNDARY'
          PAUSE
          STOP
      ELSE
          NNIVGFNALL=NNIVGCN
          NNIVGCN=NNIVGCN/3
      ENDIF
      IF (MOD(NNIVGCS,3).NE.0) THEN
          PRINT*,'WRONG IN SOUTH VG BOUNDARY'
          PAUSE
          STOP
      ELSE
          NNIVGFSALL=NNIVGCS
          NNIVGCS=NNIVGCS/3
      ENDIF
      ALLOCATE (NIVGFW(NNIVGFW),NIVGFE(NNIVGFE))
      ALLOCATE (NIVGFN(NNIVGFN),NIVGFS(NNIVGFS))
      ALLOCATE (NIVGFWALL(NNIVGFWALL,2),NIVGFEALL(NNIVGFEALL,2))
      ALLOCATE (NIVGFNALL(NNIVGFNALL,2),NIVGFSALL(NNIVGFSALL,2))
      ALLOCATE (NIVGFWS(NNIVGFWS),NIVGFES(NNIVGFES))
      ALLOCATE (NIVGFWN(NNIVGFWN),NIVGFEN(NNIVGFEN))
      NIVGFWALL=N_CTRDP1
      NIVGFEALL=N_CTRDP1
      NIVGFNALL=N_CTRDP1
      NIVGFSALL=N_CTRDP1
      VFW=0
      VFE=0
      VFS=0
      VFN=0
      VFWALL=0
      VFEALL=0
      VFSALL=0
      VFNALL=0
      VFWS=0
      VFES=0
      VFWN=0
      VFEN=0
      DO 40 I=N_CTRD_AG,N_CTRD_AG+N_CTRD_EVG
      IF (NIM1(I).LT.N_CTRDP1.AND.NIP1(I).LT.N_CTRDP1.AND.
     *NJM1(I).LT.N_CTRDP1.AND.NJP1(I).LT.N_CTRDP1) THEN
          GO TO 40 
      ELSEIF (NIM1(I).EQ.N_CTRDP1.AND.NIP1(I).LT.N_CTRDP1.AND.
     *NJM1(I).LT.N_CTRDP1.AND.NJP1(I).LT.N_CTRDP1) THEN
          VFW=VFW+1
          NIVGFW(VFW)=I
          VFWALL=VFWALL+1
          NIVGFWALL(VFWALL,1)=I
      ELSEIF (NIM1(I).LT.N_CTRDP1.AND.NIP1(I).EQ.N_CTRDP1.AND.
     *NJM1(I).LT.N_CTRDP1.AND.NJP1(I).LT.N_CTRDP1) THEN
          VFE=VFE+1
          NIVGFE(VFE)=I
          VFEALL=VFEALL+1
          NIVGFEALL(VFEALL,1)=I
      ELSEIF (NIM1(I).LT.N_CTRDP1.AND.NIP1(I).LT.N_CTRDP1.AND.
     *NJM1(I).EQ.N_CTRDP1.AND.NJP1(I).LT.N_CTRDP1) THEN
          VFS=VFS+1
          NIVGFS(VFS)=I
          VFSALL=VFSALL+1
          NIVGFSALL(VFSALL,1)=I
      ELSEIF (NIM1(I).LT.N_CTRDP1.AND.NIP1(I).LT.N_CTRDP1.AND.
     *NJM1(I).LT.N_CTRDP1.AND.NJP1(I).EQ.N_CTRDP1) THEN     
          VFN=VFN+1
          NIVGFN(VFN)=I
          VFNALL=VFNALL+1
          NIVGFNALL(VFNALL,1)=I
      ELSEIF (NIM1(I).EQ.N_CTRDP1.AND.NIP1(I).LT.N_CTRDP1.AND.
     *NJM1(I).EQ.N_CTRDP1.AND.NJP1(I).LT.N_CTRDP1) THEN
          VFWS=VFWS+1
          NIVGFWS(VFWS)=I
          VFWALL=VFWALL+1
          NIVGFWALL(VFWALL,1)=I
          VFSALL=VFSALL+1
          NIVGFSALL(VFSALL,1)=I
      ELSEIF (NIM1(I).EQ.N_CTRDP1.AND.NIP1(I).LT.N_CTRDP1.AND.
     *NJM1(I).LT.N_CTRDP1.AND.NJP1(I).EQ.N_CTRDP1) THEN
          VFWN=VFWN+1
          NIVGFWN(VFWN)=I
          VFWALL=VFWALL+1
          NIVGFWALL(VFWALL,1)=I
          VFNALL=VFNALL+1
          NIVGFNALL(VFNALL,1)=I
      ELSEIF (NIM1(I).LT.N_CTRDP1.AND.NIP1(I).EQ.N_CTRDP1.AND.
     *NJM1(I).EQ.N_CTRDP1.AND.NJP1(I).LT.N_CTRDP1) THEN
          VFES=VFES+1
          NIVGFES(VFES)=I
          VFEALL=VFEALL+1
          NIVGFEALL(VFEALL,1)=I
          VFSALL=VFSALL+1
          NIVGFSALL(VFSALL,1)=I
      ELSEIF (NIM1(I).LT.N_CTRDP1.AND.NIP1(I).EQ.N_CTRDP1.AND.
     *NJM1(I).LT.N_CTRDP1.AND.NJP1(I).EQ.N_CTRDP1) THEN
          VFEN=VFEN+1
          NIVGFEN(VFEN)=I
          VFEALL=VFEALL+1
          NIVGFEALL(VFEALL,1)=I
          VFNALL=VFNALL+1
          NIVGFNALL(VFNALL,1)=I
      ENDIF
40    CONTINUE
      RETURN
      END SUBROUTINE FIND_INTERFACE_AG