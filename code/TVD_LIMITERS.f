        REAL FUNCTION RAT(R)
        REAL R,MIN1,MIN2
        
C        MIN1=AMIN1(2.*R,1.)
C        MIN2=AMIN1(R,2.)
C        RAT=AMAX1(0.,MIN1,MIN2)
C         RAT=AMAX1(0.,AMIN1(1.,R))
        IF (R.LE.0.) THEN
        RAT=0.
        ELSE
        RAT=(R+ABS(R))/(1.+R)
        ENDIF
C         RAT=MAX(0.,AMIN1(2.,2.*R,(1.+R)/2.))
        RETURN
        END
        
        REAL FUNCTION COF(R)
        REAL R
        COF=AMAX1(0.,AMIN1(R,1.))   !MINMOD
        RETURN
        END
        
        REAL FUNCTION COF3(R,BETA)
        REAL R
        COF3=AMAX1(0.,AMIN1(2.,2*R,BETA))   !MINMOD
        RETURN
        END
        
        
        REAL FUNCTION COF5(R,BETA)
        REAL R,BETA
        COF5=AMAX1(0.,AMIN1(2.,2*R,BETA))   !MINMOD
        RETURN
        END
        
        REAL FUNCTION COFMPL(A,B,C)
        REAL A,B,C
        
        COFMPL=AMAX1(0.,AMIN1(B*A,B,C))
        RETURN
        END