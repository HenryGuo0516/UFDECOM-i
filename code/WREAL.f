#include "DEFS.h"	      
      SUBROUTINE WREAL
      USE MOD_GLOBAL
   
c      DIMENSION UC(N_CTRDP1,KB), VC(N_CTRDP1,KB)  

c      WR=0.0
c      UC=0.0
c      VC=0.0

c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO I = 1, N_CTRD
c      DO K = 1, KB
      DO concurrent (I=1:N_CTRD,K=1:KB)
      IF(DUM(I).GT.0.0) THEN      
          UC(I,K) = 0.5*U(I,K)*(H2(I)+H2(NIM1(I)))
     *    -(H3(I)+H3(NIM1(I)))/(H1(I)+H1(NIM1(I)))
     *    *(V(I,K)*DVM(I)+V(NIM1(I),K)*DVM(NIM1(I))
     *    +V(NJP1(I),K)*DVM(NJP1(I))+V(NIM1JP1(I),K)*DVM(NIM1JP1(I)))
     *    /(DVM(I)+DVM(NIM1(I))+DVM(NJP1(I))+DVM(NIM1JP1(I))+1.E-30)
          UC(I,K) = 2.*UC(I,K)/(DJ(I)+DJ(NIM1(I)))
      ELSE
          UC(I,K) = 0.0
      ENDIF
      IF(DVM(I).GT.0.0) THEN  
          VC(I,K) = 0.5*(H1(I)+H1(NJM1(I)))*V(I,K)
     *    -(H3(I)+H3(NJM1(I)))/(H2(I)+H2(NJM1(I)))
     *    *(U(I,K)*DUM(I)+U(NIP1(I),K)*DUM(NIP1(I))
     *    +U(NJM1(I),K)*DUM(NJM1(I))+U(NIP1JM1(I),K)*DUM(NIP1JM1(I)))
     *    /(DUM(I)+DUM(NIP1(I))+DUM(NJM1(I))+DUM(NIP1JM1(I))+1.E-30)
          VC(I,K) = 2.*VC(I,K)/(DJ(I)+DJ(NJM1(I)))
      ELSE
          VC(I,K) = 0.0
      ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET      


      DO K = 1, KBM1
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I = 1, N_CTRD
      DO concurrent (I=1:N_CTRD)
            TPS(I) = ZZ(K) * D(I) + ELF(I)
      ENDDO   
c!$OMP END TARGET      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO
c      DO I = 1, N_CTRD
      DO concurrent (I=1:N_CTRD) local_init(DTI)
      IF (FSM(I).GT.0.0) THEN 
      WR(I,K)=0.5*(W(I,K)+W(I,K+1))
     * +0.5*(UC(NIP1(I),K)*(TPS(NIP1(I))-TPS(I))
     * +UC(I,K)*(TPS(I)-TPS(NIM1(I)))
     * +VC(NJP1(I),K)*(TPS(NJP1(I))-TPS(I))
     * +VC(I,K)*(TPS(I)-TPS(NJM1(I))))
     * +(1.+ZZ(K))*(ELF(I)-EL(I))/DTI
      ENDIF
      ENDDO   
c!$OMP END TARGET      
      ENDDO
      
c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO N = 1,NUMEBC
c        DO K = 1,KBM1
      DO concurrent (N=1:NUMEBC,K=1:KBM1) local(IE,IC)
        IE = NETA(N)
        IC = NCON(N)
          IF (IE.EQ.NIP1(IC)) THEN
            WR(NIP1(IE),K) = WR(IE,K)
          ELSEIF (IE.EQ.NIM1(IC)) THEN
            WR(NIM1(IE),K) = WR(IE,K)
          ELSEIF (IE.EQ.NJP1(IC)) THEN
            WR(NJP1(IE),K) = WR(IE,K)
          ELSEIF (IE.EQ.NJM1(IC)) THEN
            WR(NJM1(IE),K) = WR(IE,K)
          ENDIF           
        ENDDO
c      ENDDO
c!$OMP END TARGET      
      
      
       RETURN
       END SUBROUTINE WREAL
