#include "DEFS.h"

	MODULE MOD_LAG
      SAVE
      CONTAINS
***********************************************************************
*                                                                     *      
*                        Transplanted by MaRin                        *      
*                                                                     *
***********************************************************************       
      SUBROUTINE PARTICLE_TRACK
      USE MOD_GLOBAL
      INTEGER NN,PIEND
      REAL XEND,YEND
      REAL PIXEND,PIYEND,PSEND
      NIT=7488+N_TSR+N_SEC_XY+20
      DO NN=1,N_LAG
      IF (THOUR.GE.PT_LAG(NN)) THEN
          IF (PB_LAG(NN).EQ.0) THEN
              PS_LAG(NN)=S(PI_LAG(NN),PF_LAG(NN))
              IF (PS_LAG(NN).LT.0.) THEN
                  PS_LAG(NN)=0
              ELSEIF(PS_LAG(NN).GT.35) THEN
                  PS_LAG(NN)=35
              ENDIF    
              PB_LAG(NN)=1
          ENDIF
          CALL PART_2D(PI_LAG(NN),PIEND,
     *    PIX_LAG(NN),XEND,
     *    PIY_LAG(NN),YEND,
     *    PENDX_LAG(NN),PIXEND,
     *    PENDY_LAG(NN),PIYEND,
     *    PS_LAG(NN),PSEND,PF_LAG(NN))
          
          PI_LAG(NN)=PIEND
          PIX_LAG(NN)=XEND
          PIY_LAG(NN)=YEND
          PENDX_LAG(NN)=PIXEND
          PENDY_LAG(NN)=PIYEND
          PS_LAG(NN)=PSEND
      
          SECOND = THOUR*3600.
      ENDIF
      ENDDO
      
      IF (INT(SECOND/N_OPT) .GT.NN_OPT2) THEN
      DO NN=1,N_LAG
          NIT=NIT+1
          IF (THOUR.GT.PT_LAG(NN)) THEN
          WRITE(NIT,2000) THOUR,PI_LAG(NN),PF_LAG(NN),
     *    PENDX_LAG(NN),PENDY_LAG(NN),PS_LAG(NN),S(PI_LAG(NN),1)
          ENDIF
      ENDDO    
      NN_OPT2=NN_OPT2+1  
      ENDIF
2000  FORMAT (F12.4,I10,I5,4F14.2)
      
      END SUBROUTINE PARTICLE_TRACK
      
      
      SUBROUTINE PART_2D(I0,I1,PX0,PX1,PY0,PY1,X0,X1,Y0,Y1,S0,S1,ZZZ)
      USE MOD_GLOBAL
      IMPLICIT NONE
      INTEGER I0,I1,II,KK,ZZZ,I
      REAL PX0,PX1,PY0,PY1,X0,X1,Y0,Y1,S1,S0,R1,R2,D1,D2,S2
      REAL UP,ULEFT,URIGHT,VP,VUP,VDOWN
      REAL DD1,DD2,DHMIN,KHPUP,KHPDOWN,AAMP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INTERPOLATE U
      UP=(U(I0,ZZZ)*DUM(I0)+U(NIP1(I0),ZZZ)*DUM(NIP1(I0)))/
     *(DUM(I0)+DUM(NIP1(I0))+1.E-5)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INTERPOLATE V
      VP=(V(I0,ZZZ)*DVM(I0)+V(NJP1(I0),ZZZ)*DVM(NJP1(I0)))/
     *(DVM(I0)+DVM(NJP1(I0))+1.E-5)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! NEW POSITION  
      IF (PX0.GE.0) THEN
          DD1=SQRT((XR(NIP1(I0))-XR(I0))**2+(YR(NIP1(I0))-YR(I0))**2)
          X1=X0+(XR(NIP1(I0))-XR(I0))*UP*DTI/DD1
          Y1=Y0+(YR(NIP1(I0))-YR(I0))*UP*DTI/DD1   
      ELSE
          DD1=-SQRT((XR(NIM1(I0))-XR(I0))**2+(YR(NIM1(I0))-YR(I0))**2)
          X1=X0+(XR(NIM1(I0))-XR(I0))*UP*DTI/DD1
          Y1=Y0+(YR(NIM1(I0))-YR(I0))*UP*DTI/DD1    
      ENDIF
      IF (PY0.GE.0) THEN
          DD2=SQRT((XR(NJP1(I0))-XR(I0))**2+(YR(NJP1(I0))-YR(I0))**2)
          X1=X1+(XR(NJP1(I0))-XR(I0))*VP*DTI/DD2
          Y1=Y1+(YR(NJP1(I0))-YR(I0))*VP*DTI/DD2    
      ELSE
          DD2=-SQRT((XR(NJM1(I0))-XR(I0))**2+(YR(NJM1(I0))-YR(I0))**2)
          X1=X1+(XR(NJM1(I0))-XR(I0))*VP*DTI/DD2
          Y1=Y1+(YR(NJM1(I0))-YR(I0))*VP*DTI/DD2
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! NEW I NUM
      PX1=PX0+UP*DTI/DD1
      PY1=PY0+VP*DTI/DD2
      I1=I0
      IF (PX1.GT.0.5) THEN
          I1=NIP1(I1)
      ELSEIF (PX1.LT.-0.5) THEN
          I1=NIM1(I1)
      ENDIF 
      IF (PY1.GT.0.5) THEN
          I1=NJP1(I1)
      ELSEIF(PY1.LT.-0.5) THEN
          I1=NJM1(I1)
      ENDIF
74881 CONTINUE
      DHMIN=1000000
      IF (I1.GT.N_CTRD_AG) THEN
      DO I=1,NNIAGCW
      II=NIAGCW(I)
      DD1=SQRT((X1-XR(II))**2+(Y1-YR(II))**2)
      IF (DD1.LE.DHMIN.AND.FSM(II).GT.0) THEN
          DHMIN=DD1
          I1=II
      ENDIF
      DO KK=6,8
          II=AIJ(NIAGCW(I),KK)
          DD1=SQRT((X1-XR(II))**2+(Y1-YR(II))**2)
          IF (DD1.LE.DHMIN.AND.FSM(II).GT.0) THEN
              DHMIN=DD1
              I1=II
          ENDIF
      ENDDO
      ENDDO
      DO I=1,NNIAGCE
      II=NIAGCE(I)
      DD1=SQRT((X1-XR(II))**2+(Y1-YR(II))**2)
      IF (DD1.LE.DHMIN.AND.FSM(II).GT.0) THEN
          DHMIN=DD1
          I1=II
      ENDIF
      DO KK=6,8
          II=AIJ(NIAGCE(I),KK)
          DD1=SQRT((X1-XR(II))**2+(Y1-YR(II))**2)
          IF (DD1.LE.DHMIN.AND.FSM(II).GT.0) THEN
              DHMIN=DD1
              I1=II
          ENDIF
      ENDDO
      ENDDO
      DO I=1,NNIAGCS
      II=NIAGCS(I)
      DD1=SQRT((X1-XR(II))**2+(Y1-YR(II))**2)
      IF (DD1.LE.DHMIN.AND.FSM(II).GT.0) THEN
          DHMIN=DD1
          I1=II
      ENDIF
      DO KK=6,8
          II=AIJ(NIAGCS(I),KK)
          DD1=SQRT((X1-XR(II))**2+(Y1-YR(II))**2)
          IF (DD1.LE.DHMIN.AND.FSM(II).GT.0) THEN
              DHMIN=DD1
              I1=II
          ENDIF
      ENDDO
      ENDDO
      DO I=1,NNIAGCN
      II=NIAGCN(I)
      DD1=SQRT((X1-XR(II))**2+(Y1-YR(II))**2)
      IF (DD1.LE.DHMIN.AND.FSM(II).GT.0) THEN
          DHMIN=DD1
          I1=II
      ENDIF
      DO KK=6,8
          II=AIJ(NIAGCN(I),KK)
          DD1=SQRT((X1-XR(II))**2+(Y1-YR(II))**2)
          IF (DD1.LE.DHMIN.AND.FSM(II).GT.0) THEN
              DHMIN=DD1
              I1=II
          ENDIF
      ENDDO
      ENDDO
      GO TO 74882                      
      ENDIF
      
      IF (FSM(I1).GT.0) THEN
          GO TO 74882
      ELSE
      IF (I1.EQ.I0) THEN
          IF (PX0.GE.0.AND.PY0.GE.0.AND.FSM(NIP1JP1(I0)).GT.0) THEN
              I1=NIP1JP1(I0)
              X1=XR(I1)
              Y1=YR(I1)
              GO TO 74882
          ELSEIF (PX0.GE.0.AND.PY0.LT.0.AND.FSM(NIP1JM1(I0)).GT.0) THEN
              I1=NIP1JM1(I0)
              X1=XR(I1)
              Y1=YR(I1)
              GO TO 74882
          ELSEIF (PX0.LT.0.AND.PY0.GE.0.AND.FSM(NIM1JP1(I0)).GT.0) THEN
              I1=NIM1JP1(I0)
              X1=XR(I1)
              Y1=YR(I1)
              GO TO 74882
          ELSEIF (PX0.LT.0.AND.PY0.LT.0.AND.FSM(NIM1JP1(I0)).GT.0) THEN
              I1=NIM1JM1(I0)
              X1=XR(I1)
              Y1=YR(I1)
              GO TO 74882
          ENDIF
          IF (PX0.GE.0.AND.FSM(NIP1(I0)).GT.0) THEN
              I1=NIP1(I0)
              X1=XR(I1)
              Y1=YR(I1)
              GO TO 74882
          ELSEIF (PX0.LT.0.AND.FSM(NIM1(I0)).GT.0) THEN
              I1=NIM1(I0)
              X1=XR(I1)
              Y1=YR(I1)
              GO TO 74882
          ELSEIF (PY0.GE.0.AND.FSM(NJP1(I0)).GT.0) THEN
              I1=NJP1(I0)
              X1=XR(I1)
              Y1=YR(I1)
              GO TO 74882
          ELSEIF (PY0.LT.0.AND.FSM(NJM1(I0)).GT.0) THEN
              I1=NJM1(I0)
              X1=XR(I1)
              Y1=YR(I1)
              GO TO 74882
          ENDIF
      ENDIF
      ENDIF
      DHMIN=1000000
      DO II=1,N_CTRD_AG
          DD1=SQRT((X1-XR(II))**2+(Y1-YR(II))**2)
          IF (DD1.LE.DHMIN.AND.FSM(II).GT.0) THEN
              DHMIN=DD1
              I1=II
          ENDIF
      ENDDO
      X1=XR(I1)
      Y1=YR(I1)
74882 CONTINUE  
      R1=X1-XR(I1)
      R2=Y1-YR(I1)
      D1=XR(NIP1(I1))-XR(I1)
      D2=YR(NIP1(I1))-YR(I1)
      DD1=(R1*D1+R2*D2)/SQRT(D1**2+D2**2)
      IF (DD1.GE.0) THEN
          PX1=DD1/SQRT(D1**2+D2**2)
      ELSE
          D1=XR(NIM1(I1))-XR(I1)
          D2=YR(NIM1(I1))-YR(I1)
          DD1=(R1*D1+R2*D2)/SQRT(D1**2+D2**2)
          PX1=-DD1/SQRT(D1**2+D2**2)
      ENDIF
      D1=XR(NJP1(I1))-XR(I1)
      D2=YR(NJP1(I1))-YR(I1)
      DD1=(R1*D1+R2*D2)/SQRT(D1**2+D2**2)
      IF (DD1.GT.0) THEN
          PY1=DD1/SQRT(D1**2+D2**2)
      ELSE
          D1=XR(NJM1(I1))-XR(I1)
          D2=YR(NJM1(I1))-YR(I1)
          DD1=(R1*D1+R2*D2)/SQRT(D1**2+D2**2)
          PY1=-DD1/SQRT(D1**2+D2**2)
      ENDIF
      KHPUP=KH(I1,ZZZ)
      KHPDOWN=KH(I1,ZZZ+1)
      IF (KHPUP.GT.10/DTI)  KHPUP=10/DTI !KH 过大时容易爆
      IF (KHPDOWN.GT.10/DTI)  KHPDOWN=10/DTI !KH 过大时容易爆
      IF (KHPUP.LT.1/DTI)  KHPUP=1/DTI !KH 过小时扩散不开
      IF (KHPDOWN.LT.1/DTI)  KHPDOWN=1/DTI !KH 过小时扩散不开
      S1=S0
      DO I=1,10
      S1=S1+DTI/10*
     *(KHPUP*(S(I1,ZZZ)-S1)-KHPDOWN*(S1-S(I1,ZZZ)))!理论上，是ds/dz，这里不合理,MaRin
      IF (ABS(S1-S(I1,ZZZ)).LT.0.01) EXIT
      ENDDO  
      
      END SUBROUTINE PART_2D 
	END MODULE MOD_LAG
	  
	  