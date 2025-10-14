
      SUBROUTINE DENSD(TD,SD,RHOD)
C     VERSION(03/02/90)
C
      USE MOD_GLOBAL
C
C         THIS FUNCTION COMPUTES DENSITY-1.0
C
            RHOD = SD * SD * SD *
     *          6.76786136E-6 - SD * SD * 4.8249614E-4 +
     *          SD * 8.14876577E-1 - 0.22584586E0
C
            RHOD = RHOD * (TD*TD*
     *          TD*1.667E-8-TD*TD*8.164E-7+
     *          TD*1.803E-5)
C
            RHOD = RHOD + 1. - TD * TD *
     *          TD * 1.0843E-6 + TD * TD *
     *          9.8185E-5 - TD * 4.786E-3
C
            RHOD = RHOD * (SD*SD*
     *          SD*6.76786136E-6-SD*SD*
     *          4.8249614E-4+SD*8.14876577E-1+3.895414E-2)
C
            RHOD = RHOD - (TD-3.98) ** 2 * (
     *          TD+283.) / (503.57*(TD+67.26))
C
              RHOD = RHOD * 1.E-3
      RETURN
      END SUBROUTINE DENSD
