#include "DEFS.h" 
      SUBROUTINE ADVQ(Q,QF)
      USE MOD_GLOBAL      
      DIMENSION Q(N_CTRDP1,KB), QF(N_CTRDP1,KB)
c      DIMENSION AAMAX(N_CTRD,KB), AAMAY(N_CTRD,KB)

      
c      QF= 0.0
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 1, KB
c      DO I = 1, N_CTRDP1
      DO concurrent (I=1:N_CTRDP1,K=1:KB)
          QF(I,K) = 0.0
      ENDDO
c      ENDDO
c!$OMP END TARGET      
          
c      DO K = 2,KBM1
c#ifdef OMP
c!$OMP PARALLEL DO
c#endif
c      DO I = 1,N_CTRD
c          XFLUX(I,K) = 
c     *    .25*(Q(I,K)+Q(NIM1(I),K))*(XMFLUX(I,K)+XMFLUX(I,K-1))
c          YFLUX(I,K) = 
c     *    .25*(Q(I,K)+Q(NJM1(I),K))*(YMFLUX(I,K)+YMFLUX(I,K-1))
c          AAMAX(I,K)=.5*(AAM(I,K)+AAM(NIM1(I),K))
c          AAMAY(I,K)=.5*(AAM(I,K)+AAM(NJM1(I),K))
c      ENDDO
c#ifdef OMP
c!$OMP END PARALLEL DO
c#endif
c      ENDDO
        
c      DO K = 2,KBM1
c#ifdef OMP
c!$OMP PARALLEL DO
c#endif
c      DO I = 1,N_CTRD
c          XFLUX(I,K)=XFLUX(I,K)-AAMAX(I,K)
c     *    *(H(I)+H(NIM1(I)))*(Q(I,K)-Q(NIM1(I),K))
c     *    *DUM(I)/(H1(I)+H1(NIM1(I)))
c     *    *.5*(H2(I)+H2(NIM1(I)))
c          YFLUX(I,K)=YFLUX(I,K)-AAMAY(I,K)
c     *    *(H(I)+H(NJM1(I)))*(Q(I,K)-Q(NJM1(I),K))
c     *    *DVM(I)/(H2(I)+H2(NJM1(I)))
c     *    *.5*(H1(I)+H1(NJM1(I)))
c      ENDDO
c#ifdef OMP
c!$OMP END PARALLEL DO
c#endif
c      ENDDO
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 2,KBM1
c      DO I = 1,N_CTRD
      DO concurrent (I=1:N_CTRD,K=2:KBM1)
          XFLUX(I,K) = 
     *    .25*(Q(I,K)+Q(NIM1(I),K))*(XMFLUX(I,K)+XMFLUX(I,K-1))
     *    -.5*(AAM(I,K)+AAM(NIM1(I),K))
     *    *(H(I)+H(NIM1(I)))*(Q(I,K)-Q(NIM1(I),K))
     *    *DUM(I)/(H1(I)+H1(NIM1(I)))
     *    *.5*(H2(I)+H2(NIM1(I)))
          YFLUX(I,K) = 
     *    .25*(Q(I,K)+Q(NJM1(I),K))*(YMFLUX(I,K)+YMFLUX(I,K-1))
     *    -.5*(AAM(I,K)+AAM(NJM1(I),K))
     *    *(H(I)+H(NJM1(I)))*(Q(I,K)-Q(NJM1(I),K))
     *    *DVM(I)/(H2(I)+H2(NJM1(I)))
     *    *.5*(H1(I)+H1(NJM1(I)))
      ENDDO
c      ENDDO
c!$OMP END TARGET      
      
      IF (NUMEBC.NE.0) THEN
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO N = 1,NUMEBC
c      DO K = 2,KBM1
      DO concurrent (N=1:NUMEBC,K=2:KBM1) local(IE,IC)
#if defined TIDE_EL || defined TIDE_FLUX  || defined TIDE_FLATHER
          IE = NETA(N)
          IC = NCON(N)
          IF (IE==NIP1(IC)) THEN  !EAST BOUNDARY
              XFLUX(IE,K) = 0.5*Q(IC,K)*(XMFLUX(IE,K)+XMFLUX(IE,K-1))
          ELSEIF (IE==NIM1(IC)) THEN  !WEST BOUNDARY
              XFLUX(IC,K) = 0.5*Q(IC,K)*(XMFLUX(IC,K)+XMFLUX(IC,K-1))
          ELSEIF (IE==NJP1(IC)) THEN  !NORTH BOUNDARY
              YFLUX(IE,K) = 0.5*Q(IC,K)*(YMFLUX(IE,K)+YMFLUX(IE,K-1))
          ELSEIF (IE==NJM1(IC)) THEN  !SOUTH BOUNDARY
              YFLUX(IC,K) = 0.5*Q(IC,K)*(XMFLUX(IC,K)+XMFLUX(IC,K-1))
          ENDIF

#endif
      ENDDO
c      ENDDO
c!$OMP END TARGET      
      ENDIF
      
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 2,KBM1
c      DO I = 1,N_CTRD_AG
      DO concurrent (I=1:N_CTRD_AG,K=2:KBM1) local(QF0)
      IF(FSMADD(I).GT.0.0) THEN
          QF0=(W(I,K-1)*Q(I,K-1)-W(I,K+1)*Q(I,K+1))/(DZ(K)+DZ(K-1))
     * *DJ(I)+XFLUX(NIP1(I),K)-XFLUX(I,K)+YFLUX(NJP1(I),K)-YFLUX(I,K) 
          QF(I,K)=((H(I)+EL(I))*DJ(I)*Q(I,K)-DTI*QF0)/(D(I)*DJ(I))
      ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET 
      
      RETURN
      END SUBROUTINE ADVQ