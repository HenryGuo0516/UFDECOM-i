      SUBROUTINE CAL_EL
!----------------------------------------------------------------------
! COMPUTING EL
! 
! BY SYLVERS DING, OCT. 19, 2018.
!----------------------------------------------------------------------
      USE MOD_GLOBAL
!----------------------------------------------------------------------
!                  SOLVE CONTINUITY EQUATION FOR EL
!----------------------------------------------------------------------
      XMFLUX=0.
      YMFLUX=0.
      
      DO K = 1,KBM1
#ifdef OMP
!$OMP PARALLEL DO PRIVATE(UC1,VC1)
#endif
      DO I = 1,N_CTRD
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
#ifdef OMP
!$OMP END PARALLEL DO
#endif
      ENDDO
      
#ifdef NIFluxCon
      CALL FluxCon(XMFLUX,YMFLUX)
      CALL Flux2uv
#endif
      
#ifdef WEIR
      CALL ADDWEIR(2)
      CALL ADDWEIR(3)
#endif
#if defined LBCel_gra0 || defined LBCel_cha || defined LBCel_rad
      CALL BCOND(10)
#else
      CALL BCOND(3)
#endif
      CALL BCOND(4) !FOR FLUX BOUNDARY
      BBBB= 0.0
      DO K = 1,KBM1
#ifdef OMP
!$OMP PARALLEL DO
#endif
      DO I = 1,N_CTRD_AG 
      IF (FSMADD(I).NE. 0.0) THEN     
          BBBB(I) = BBBB(I)+DZ(K)*(XMFLUX(NIP1(I),K)-XMFLUX(I,K)
     *    +YMFLUX(NJP1(I),K)-YMFLUX(I,K)) 
      ENDIF
      ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif
      ENDDO
#ifdef OMP
!$OMP PARALLEL DO
#endif
      DO I = 1,N_CTRD_AG
      IF (FSMADD(I) .NE. 0.0) THEN     
          ELF(I) = EL(I) - DTI * BBBB(I) / DJ(I)
      ENDIF
      ENDDO
#ifdef OMP
!$OMP END PARALLEL DO
#endif

#if defined LBCel_gra0 || defined LBCel_cha || defined LBCel_rad
      CALL BCOND(10)
#else
      CALL BCOND(3)
#endif
      RETURN
      END