#include "DEFS.h"
      SUBROUTINE RESTART
      USE MOD_GLOBAL
#ifdef MODULE_SED
      USE MOD_SED
#endif
!======================================================================
!						     COLD/HOT START
!======================================================================
!----------------------------------------------------------------------
      IF (RESTAR .EQ. 'cold') THEN   
#if defined MODULE_SAL ||  defined MODULE_TMP
          CALL TANDS3C
#endif

#ifdef MODULE_SED
		CALL INITSED
#  ifdef TAUTYP_USERDEF
		CALL INIT_TAUE
#  endif
#endif
	    CALL DENS  
	    TDAY = 0.
	    CALL INIT_COF
#ifdef Z0BTYP_USERDEF
          CALL INIT_CBCCOF
#endif
!----------------------------------------------------------------------
      ELSE
          
#ifdef MODULE_SED
#  ifdef TAUTYP_USERDEF
          CALL INIT_TAUE
#  endif
#endif         
      OPEN (IURRS,FORM='UNFORMATTED',FILE='RESTART',STATUS='UNKNOWN')
      READ (IURRS) TDAY
      READ (IURRS) NSTEP, DZR, Z, ZZ, DZ, DZZ, H, H1, H2, H3, DJ,
     * D, EL, ARU, ARV,DUM, DVM, FSM, COR, WUBOT, WVBOT,
     * KM, KH, KQ, Q2, Q2L, L, U, V, W, T, S, RHO
     * FSMADD, ELF
      READ (IURRS) UF,VF,HU,HV,TPS,AAUU,AAVV,BBBB,AAZZ
      READ (IURRS) XFLUX,YFLUX,XMFLUX,YMFLUX,AAM,DU,DV
      READ (IURRS) XXI,XETA,YXI,YETA
      READ (IURRS) CBC_COF,KM_COF,KH_COF
#ifdef MODULE_SED
      READ (IURRS) SED,D50,TAU,TAUE,TAUD,TAUE0,TAUD0,
     *QDEP,QERO,ZBED,ZBEDD,SUMH !LXY
#endif
#ifdef MODULE_WAVE
      READ (IURRS) HSIG,TSIG,WAVEDIR
#endif     
#ifdef MODULE_MATERIAL
      READ (IURRS) CONCENTRATION
#endif 
          CLOSE (IURRS)
 	    DO I = 1,JHM
          IF (TDAY.GE.HIST(I,2)/24.) THEN
              LOG_ARC(I) = .FALSE.
          ENDIF
          ENDDO
#if defined MODULE_SAL ||  defined MODULE_TMP
		CALL DENS 
          CALL INIT_COF
          CALL INIT_CBCCOF
#endif
      ENDIF 
      RETURN
      END SUBROUTINE RESTART