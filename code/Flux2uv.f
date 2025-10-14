      SUBROUTINE FLUX2UV
      
!=======================================================================
! THIS SUBROUTINE RECALCULATE VELOCITY FROM FLUX AT THE NESTED INTERFACES
! REVERSELY.
!
! SYLVERS DING
! 2020-11-10
!=======================================================================
      USE MOD_GLOBAL

      IMPLICIT NONE
      INTEGER I, J, K, II, JJ, KK
!=======================================================================
!                     WESTERN NESTED INTERFACE
!=======================================================================   
      DO I = 1,NNIAGCW
      II = NIP1(NIAGCW(I))
      IF (DUM(II).NE.0.0) THEN
      DO K = 1,KBM1
          UF(II,K) = XMFLUX(II,K)/DU(II)
          UF(II,K) = UF(II,K)
     *    +0.25*(H3(II)+H3(NIM1(II)))/(H1(II)+H1(NIM1(II)))
     *    *(VF(II,K)+VF(NIM1(II),K)+VF(NJP1(II),K)+VF(NIM1JP1(II),K)) 
          UF(II,K) = UF(II,K)/(0.5*(H2(II)+H2(NIM1(II))))
      ENDDO
      ENDIF
      ENDDO
!=======================================================================
!                   END: WESTERN NESTED INTERFACE
!=======================================================================

      
!=======================================================================
!                     EASTERN NESTED INTERFACE
!=======================================================================
      DO I = 1,NNIAGCE
      II = NIAGCE(I)
      IF (DUM(II).NE.0.0) THEN
      DO K = 1,KBM1
          UF(II,K) = XMFLUX(II,K)/DU(II)
          UF(II,K) = UF(II,K)
     *    +0.25*(H3(II)+H3(NIM1(II)))/(H1(II)+H1(NIM1(II)))
     *    *(VF(II,K)+VF(NIM1(II),K)+VF(NJP1(II),K)+VF(NIM1JP1(II),K)) 
          UF(II,K) = UF(II,K)/(0.5*(H2(II)+H2(NIM1(II))))
      ENDDO
      ENDIF
      ENDDO
!=======================================================================
!                    END: EASTERN NESTED INTERFACE
!=======================================================================
      
      
!=======================================================================
!                     NORTHERN NESTED INTERFACE
!=======================================================================
      DO I = 1,NNIAGCN
      II = NIAGCN(I)
      IF (DVM(II).NE.0.0) THEN
      DO K = 1,KBM1
          VF(II,K) = YMFLUX(II,K)/DV(II)
          VF(II,K) = VF(II,K)
     *    +0.25*(H3(II)+H3(NJM1(II)))/(H2(II)+H2(NJM1(II)))
     *    *(UF(II,K)+UF(NIP1(II),K)+UF(NJM1(II),K)+UF(NIP1JM1(II),K))
          VF(II,K) = VF(II,K)/(0.5*(H1(II)+H1(NJM1(II))))
      ENDDO
      ENDIF
      ENDDO
!=======================================================================
!                    END: NORTHERN NESTED INTERFACE
!=======================================================================
      
      
!=======================================================================
!                     SOUTHERN NESTED INTERFACE
!=======================================================================      
      DO I = 1,NNIAGCS
      II = NJP1(NIAGCN(I))
      IF (DVM(II).NE.0.0) THEN
      DO K = 1,KBM1
          VF(II,K) = YMFLUX(II,K)/DV(II)
          VF(II,K) = VF(II,K)
     *    +0.25*(H3(II)+H3(NJM1(II)))/(H2(II)+H2(NJM1(II)))
     *    *(UF(II,K)+UF(NIP1(II),K)+UF(NJM1(II),K)+UF(NIP1JM1(II),K))
          VF(II,K) = VF(II,K)/(0.5*(H1(II)+H1(NJM1(II))))
      ENDDO
      ENDIF
      ENDDO
!=======================================================================
!                    END: SOUTHERN NESTED INTERFACE
!=======================================================================
      RETURN
      END SUBROUTINE FLUX2UV