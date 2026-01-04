#include "DEFS.h"
      
      SUBROUTINE UVCon
      
!======================================================================
! This subroutine applies update-mix_low scheme (Debreu et al., 2012)
! to the northern and eastern nested interface.
!
! Sylvers Ding
! 2020-10-20
!======================================================================
      
      USE MOD_GLOBAL
      
      implicit none

      integer i, j, k, ii, jj ,kk
      integer ind, ntemp, rfn
      integer nref(5)
      real xb(4),yb(4)
      real sumall, weiall
  
      
!======================================================================
!                     Western nested interface
!======================================================================
      ! Recalculate UF at the western interface
      DO i = 1,NniAGcW
        ii = NIP1(niAGcW(i))
        
        DO k = 1,4
	    xb(k) = xnode(k,ii)
	    yb(k) = ynode(k,ii)
        ENDDO
        
        DO j = 1,N_CTRD_AG
          IF (COORD.eq.'XY') THEN
	  	  CALL INSIDE(xr(j),yr(j),xb,yb,nb,ind)
	    ELSE
	  	  CALL INSIDE(lon(j),lat(j),xb,yb,nb,ind)
          ENDIF
          
          IF (ind.eq.1) THEN
            nref(1) = j
            exit
          ENDIF
        ENDDO
        
        ntemp = nref(1)
        DO j = 2,5
          ntemp = njp1(ntemp)
          IF (COORD.eq.'XY') THEN
	  	  CALL INSIDE(xr(ntemp),yr(ntemp),xb,yb,nb,ind)
	    ELSE
	  	  CALL INSIDE(lon(ntemp),lat(ntemp),xb,yb,nb,ind)
          ENDIF
          
          IF (ind.eq.1) THEN
            nref(j) = ntemp
          ELSE
            ! refining multiple
            rfn = j-1
            exit
          ENDIF
        ENDDO
        
        DO k = 1,KBM1
          weiall = 0.
          sumall = 0.
          DO kk = 1,rfn
            weiall = weiall+H2(nref(kk))
            sumall = sumall+UF(nref(kk),k)*H2(nref(kk))
          ENDDO
          
          IF (weiall.eq.0.) THEN
            UF(ii,k) = 0.
          ELSE
            UF(ii,k) = sumall/weiall
          ENDIF
        ENDDO
        
      ENDDO
      
!======================================================================
!                   End: Western nested interface
!======================================================================
      
      
!======================================================================
!                     Eastern nested interface
!======================================================================
      
      ! Recalculate UF at the eastern interface
      DO i = 1,NniAGcE
        ii = niAGcE(i)
        
        DO k = 1,4
	    xb(k) = xnode(k,ii)
	    yb(k) = ynode(k,ii)
        ENDDO
        
        DO j = 1,N_CTRD_EVG
          IF (COORD.eq.'XY') THEN
	  	  CALL INSIDE(xr(N_CTRD_AG+j),yr(N_CTRD_AG+j),xb,yb,nb,ind)
	    ELSE
	  	  CALL INSIDE(lon(N_CTRD_AG+j),lat(N_CTRD_AG+j),xb,yb,nb,ind)
          ENDIF
          
          IF (ind.eq.1) THEN
            nref(1) = N_CTRD_AG+j
            exit
          ENDIF
        ENDDO
        
        ntemp = nref(1)
        DO j = 2,5
          ntemp = njp1(ntemp)
          IF (COORD.eq.'XY') THEN
	  	  CALL INSIDE(xr(ntemp),yr(ntemp),xb,yb,nb,ind)
	    ELSE
	  	  CALL INSIDE(lon(ntemp),lat(ntemp),xb,yb,nb,ind)
          ENDIF
          
          IF (ind.eq.1) THEN
            nref(j) = ntemp
          ELSE
            ! refining multiple
            rfn = j-1
            exit
          ENDIF
        ENDDO
        
        
        DO k = 1,KBM1
          weiall = 0.
          sumall = 0.
          DO kk = 1,rfn
            weiall = weiall+H2(nref(kk))
            sumall = sumall+UF(nref(kk),k)*H2(nref(kk))
          ENDDO
          
          IF (weiall.eq.0.) THEN
            UF(ii,k) = 0.
          ELSE
            UF(ii,k) = sumall/weiall
          ENDIF
        ENDDO
          
      ENDDO
      
      
      
      
!======================================================================
!                    End: Eastern nested interface
!======================================================================
      
      
!======================================================================
!                     Northern nested interface
!======================================================================
      ! Recalculate VF at the Northern interface
      DO i = 1,NniAGcN
        ii = niAGcN(i)
        
        DO k = 1,4
	    xb(k) = xnode(k,ii)
	    yb(k) = ynode(k,ii)
        ENDDO
        
        DO j = 1,N_CTRD_EVG
          IF (COORD.eq.'XY') THEN
	  	  CALL INSIDE(xr(N_CTRD_AG+j),yr(N_CTRD_AG+j),xb,yb,nb,ind)
	    ELSE
	  	  CALL INSIDE(lon(N_CTRD_AG+j),lat(N_CTRD_AG+j),xb,yb,nb,ind)
          ENDIF
          
          IF (ind.eq.1) THEN
            nref(1) = N_CTRD_AG+j
            exit
          ENDIF
        ENDDO
        
        ntemp = nref(1)
        DO j = 2,5
          ntemp = nip1(ntemp)
          IF (COORD.eq.'XY') THEN
	  	  CALL INSIDE(xr(ntemp),yr(ntemp),xb,yb,nb,ind)
	    ELSE
	  	  CALL INSIDE(lon(ntemp),lat(ntemp),xb,yb,nb,ind)
          ENDIF
          
          IF (ind.eq.1) THEN
            nref(j) = ntemp
          ELSE
            ! refining multiple
            rfn = j-1
            exit
          ENDIF
        ENDDO
        
        DO k = 1,KBM1
          weiall = 0.
          sumall = 0.
          DO kk = 1,rfn
            weiall = weiall+H1(nref(kk))
            sumall = sumall+VF(nref(kk),k)*H1(nref(kk))
          ENDDO
          
          IF (weiall.eq.0.) THEN
            VF(ii,k) = 0.
          ELSE
            VF(ii,k) = sumall/weiall
          ENDIF
        ENDDO
        
      ENDDO
      
!======================================================================
!                    End: Northern nested interface
!======================================================================
      
      
!======================================================================
!                     Southern nested interface
!======================================================================
      ! Recalculate VF at the Southern interface
      DO i = 1,NniAGcS
        ii = NJP1(niAGcN(i))
        
        DO k = 1,4
	    xb(k) = xnode(k,ii)
	    yb(k) = ynode(k,ii)
        ENDDO
        
        DO j = 1,N_CTRD_AG
          IF (COORD.eq.'XY') THEN
	  	  CALL INSIDE(xr(j),yr(j),xb,yb,nb,ind)
	    ELSE
	  	  CALL INSIDE(lon(j),lat(j),xb,yb,nb,ind)
          ENDIF
          
          IF (ind.eq.1) THEN
            nref(1) = j
            exit
          ENDIF
        ENDDO
        
        ntemp = nref(1)
        DO j = 2,5
          ntemp = nip1(ntemp)
          IF (COORD.eq.'XY') THEN
	  	  CALL INSIDE(xr(ntemp),yr(ntemp),xb,yb,nb,ind)
	    ELSE
	  	  CALL INSIDE(lon(ntemp),lat(ntemp),xb,yb,nb,ind)
          ENDIF
          
          IF (ind.eq.1) THEN
            nref(j) = ntemp
          ELSE
            ! refining multiple
            rfn = j-1
            exit
          ENDIF
        ENDDO
        
        DO k = 1,KBM1
          weiall = 0.
          sumall = 0.
          DO kk = 1,rfn
            weiall = weiall+H1(nref(kk))
            sumall = sumall+VF(nref(kk),k)*H1(nref(kk))
          ENDDO
          
          IF (weiall.eq.0.) THEN
            VF(ii,k) = 0.
          ELSE
            VF(ii,k) = sumall/weiall
          ENDIF
        ENDDO
        
      ENDDO
      
!======================================================================
!                    End: Southern nested interface
!======================================================================
      
      
      RETURN
      
      END SUBROUTINE UVCon