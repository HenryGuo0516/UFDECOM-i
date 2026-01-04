      SUBROUTINE MAXMIN(F,FMASK,NN,FMAX,FMIN)
!    VERSION(03/02/90)

      DIMENSION F(NN), FMASK(NN)

      DATA BIG /1.E20/

      FMAX = -BIG
      FMIN = BIG

      DO I=1,NN
        IF (FMASK(I).NE.0.) THEN
          FMAX = AMAX1(FMAX,F(I))
          FMIN = AMIN1(FMIN,F(I))
        ENDIF
      ENDDO
      
      RETURN
      END