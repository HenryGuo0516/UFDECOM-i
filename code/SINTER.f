      PURE SUBROUTINE SINTER(X,A,Y,B,M,N)
c      DIMENSION X(M), A(M), Y(N), B(N)
!$omp declare target
      INTEGER, INTENT(IN) ::  M, N
      REAL, INTENT(IN) ::  X(M), A(M), Y(N)
      REAL, INTENT(OUT) ::  B(N)
c      REAL, PRIVATE :: I1
      LOGICAL RUNNEXT
      
c ----- rivision 1:   
c      I=1
c      J=1
c      DO WHILE (I<=N)
c          IF (Y(I)>X(1)) THEN 
c              B(I)=A(1)+((A(1)-A(2))/(X(1)-X(2)))
c     *      *(Y(I)-X(1))
c              I=I+1
c          ELSEIF (J>M) THEN
c              B(I)=A(M)
c              I=I+1
c          ELSEIF (Y(I)>X(J)) THEN
c              B(I) = A(J-1) - (A(J-1)-
c     *        A(J)) * (X(J-1)-Y(I)) / (X(J-1)-X(J))
c              I=I+1
c          ELSE
c              J=J+1
c          ENDIF
c      ENDDO
      
c ----- original code 0:
      
c      DO 10 I = 1, N
c        IF (Y(I).GT.X(1)) B(I)=A(1)+((A(1)-A(2))/(X(1)-X(2)))
c     *      *(Y(I)-X(1))
c        IF (Y(I).LT.X(M)) B(I)=A(M)
c10    CONTINUE

c      NM = M - 1
c      DO 30 I = 1, N
c        DO 20 J = 1, NM
c          IF (Y(I).LE.X(J).AND.Y(I).GE.X(J+1)) B(I) = A(J) - (A(J)-
c     *        A(J+1)) * (X(J)-Y(I)) / (X(J)-X(J+1))
c20      CONTINUE
c30    CONTINUE

c ----- rivision 2:
      
      I=1
      RUNNEXT=(Y(1)>X(1))
      DO WHILE (RUNNEXT)
          B(I)=A(1)+((A(1)-A(2))/(X(1)-X(2)))*(Y(I)-X(1))
          I=I+1
          if (I<=N) then
                RUNNEXT=(Y(I)>X(1))
          else
                RUNNEXT=.false.
          endif
      ENDDO
      DO J=1,M-1
          if (I<=N) then
              RUNNEXT=(Y(I)>=X(J+1))
          else
              RUNNEXT=.false.
          endif
          
          DO WHILE (RUNNEXT)
              B(I) = A(J) - (A(J)-A(J+1))
     *         * (X(J)-Y(I)) / (X(J)-X(J+1))
              I=I+1
              if (I<=N) then
                  RUNNEXT=(Y(I)>=X(J+1))
              else
                  RUNNEXT=.false.
              endif
          Enddo
      Enddo
      I1=I
      DO I=I1,N
          B(I)=A(M)
      ENDDO

      RETURN
      END SUBROUTINE SINTER