      SUBROUTINE NO_INVERSION !LTY
      USE MOD_GLOBAL 
c      USE DFPORT
      
      
      IMPLICIT NONE
      INTEGER,PARAMETER::IMNMC=337,JMNMC=225,KBM1NMC=10,KBNMC=11
      INTEGER COUNTNMC,NUMNMC
      
	CHARACTER*80 FNNMC,FNSNMC,F1NMC,F2NMC
      CHARACTER*80 CLASSICNMC
	REAL, ALLOCATABLE :: DIFNMC(:,:,:)
	REAL, ALLOCATABLE :: DIF1NMC(:,:)
c	REAL, ALLOCATABLE :: IINMC(:),JJNMC(:)
	INTEGER, ALLOCATABLE :: IINMC(:),JJNMC(:)
	REAL, ALLOCATABLE :: ANMC(:,:),BNMC(:),A1NMC(:),A2NMC(:),A3NMC(:,:)
     *,B1NMC(:)
      REAL, ALLOCATABLE :: AINMC(:,:)
	REAL, ALLOCATABLE :: XNNMC(:),XFNMC(:)
      REAL, ALLOCATABLE :: SSNMC(:,:),SFNMC(:,:)
      REAL, ALLOCATABLE :: HNMC(:)
	REAL, ALLOCATABLE :: RINMC(:)
	REAL, ALLOCATABLE :: Y0NMC(:,:)
	REAL:: F22NMC
	REAL:: THOURNMC,XXXXNMC
      INTEGER INMC,JNMC,KNMC,MNMC,NNMC,I1NMC,J1NMC,RELNMC
      REAL A11NMC,A22NMC,B11NMC,SSSNMC,DIIFNMC,DIIF1NMC
      LOGICAL RE_INMC

      CLASSICNMC=''  ! CLASSICNMC\ 
      
      COUNTNMC=0
      
	F22NMC=0.0000001
	
      OPEN(113,FILE='input/RI.DAT')
      
c      DO WHILE(.NOT.EOF(113))
      DO WHILE(1)
          READ(113,*,end=301)
          COUNTNMC=COUNTNMC+1
      ENDDO
301   continue
      REWIND(113)
      NUMNMC=COUNTNMC  

      ALLOCATE(RINMC(NUMNMC))
      ALLOCATE(Y0NMC(NUMNMC,KBM1NMC))
       
	DO INMC=1,NUMNMC
	    READ(113,*) RINMC(INMC)
      END DO
      CLOSE(113)

	OPEN(10,FILE='output/NMC/'//TRIM(CLASSICNMC)//'NN1.DAT')
      OPEN(20,FILE='output/NMC/'//TRIM(CLASSICNMC)//'H.DAT')
      OPEN(30,FILE='output/NMC/29MODEL_RESULT.DAT')
	DO INMC=1,NUMNMC
	    READ(30,*)XXXXNMC,XXXXNMC,(Y0NMC(INMC,KNMC),KNMC=1,KBM1NMC)
	END DO

	OPEN (199,FILE='input/CH_HZBC_SI_BEF.DAT')    
    	OPEN (200,FILE='input/CH_HZBC_SI.DAT')   
      
      DO I1NMC=1,N_CTRD
          WRITE(199,'(100G)')I1NMC,(S(I1NMC,KNMC),
     *KNMC=1,KBM1NMC)
      ENDDO
      DO INMC=1,NUMNMC !观测点个数
	    READ(10,*)NNMC
	    ALLOCATE (HNMC(NNMC))
	    ALLOCATE (XNNMC(NNMC))
          ALLOCATE (XFNMC(NNMC))
	    ALLOCATE (SSNMC(NNMC,KBM1NMC))
	    ALLOCATE (SFNMC(NNMC,KBM1NMC))
          ALLOCATE (DIF1NMC(NNMC,NNMC))
	    ALLOCATE (ANMC(NNMC,NNMC))
	    ALLOCATE (A1NMC(NNMC))
	    ALLOCATE (A2NMC(NNMC))
	    ALLOCATE (A3NMC(NNMC,NNMC))
	    ALLOCATE (AINMC(NNMC,NNMC))
	    ALLOCATE (BNMC(NNMC))
	    ALLOCATE (B1NMC(NNMC))
	    ALLOCATE (IINMC(NNMC))
	    ALLOCATE (JJNMC(NNMC))

	    READ(20,*)
	    DO KNMC=1,NNMC
	        READ(20,*)HNMC(KNMC)
	    END DO
      
	    DO MNMC=1,NNMC
	        READ(10,*)IINMC(MNMC)
	        DO KNMC=1,KBM1NMC
	            SSNMC(MNMC,KNMC)=S(IINMC(MNMC),KNMC)
	        END DO
          END DO
      
	    WRITE(F1NMC,1000)INMC
          WRITE(F2NMC,1000)INMC*100

          FNNMC='output/NMC/DIFF/'//TRIM(CLASSICNMC)//
     *TRIM(F1NMC)//'.DAT'
	    FNSNMC='output/NMC/B/'//TRIM(F2NMC)//'.DAT'

          OPEN(60,FILE=FNNMC,FORM='UNFORMATTED',CONVERT='BIG_ENDIAN')
	    OPEN(111,FILE=FNSNMC)
          
c          RE_INMC=SYSTEM("MKDIR "//"output/NMC/B/"//TRIM(F2NMC)//"/")
          CALL SYSTEM("mkdir output/NMC/B/"//TRIM(F2NMC)//"/")
          
          OPEN(55,FILE='output/NMC/B/'//TRIM(F2NMC)//'/A.DAT')
          OPEN(551,FILE='output/NMC/B/'//TRIM(F2NMC)//'/A1.DAT')
          OPEN(552,FILE='output/NMC/B/'//TRIM(F2NMC)//'/A2.DAT')
          OPEN(553,FILE='output/NMC/B/'//TRIM(F2NMC)//'/A3.DAT')
          OPEN(555,FILE='output/NMC/B/'//TRIM(F2NMC)//'/BB.DAT')
          OPEN(66,FILE='output/NMC/B/'//TRIM(F2NMC)//'/B.DAT')
          OPEN(666,FILE='output/NMC/B/'//TRIM(F2NMC)//'/XB.DAT')

	    DO KNMC=1,KBM1NMC
	        READ(60)DIF1NMC

              DO I1NMC=1,NNMC
	            A11NMC=0.
                  DO J1NMC=1,NNMC 
                      A11NMC=DIF1NMC(I1NMC,J1NMC)*HNMC(J1NMC)+A11NMC
	            END DO
	            A1NMC(I1NMC)=A11NMC
	        ENDDO

	        DO I1NMC=1,NNMC
	            A22NMC=0.
                  DO J1NMC=1,NNMC 
                      A22NMC=DIF1NMC(I1NMC,J1NMC)*HNMC(J1NMC)+A22NMC
	            END DO
	            A2NMC(I1NMC)=A22NMC
	        ENDDO

	        DO I1NMC=1,NNMC
              DO J1NMC=1,NNMC
	            A3NMC(I1NMC,J1NMC)=A1NMC(I1NMC)*A2NMC(J1NMC)*RINMC(INMC)
	        ENDDO
	        ENDDO

	        DO I1NMC=1,NNMC
              DO J1NMC=1,NNMC
	            ANMC(I1NMC,J1NMC)=A3NMC(I1NMC,J1NMC)+DIF1NMC(I1NMC,J1NMC)
	        ENDDO
	        ENDDO

	
              B11NMC=0.
              DO I1NMC=1,NNMC		
	            B11NMC=SSNMC(I1NMC,KNMC)*HNMC(I1NMC)+B11NMC
	        ENDDO
	
	        DO I1NMC=1,NNMC
	            BNMC(I1NMC)=(Y0NMC(INMC,KNMC)-B11NMC)*RINMC(INMC)*A1NMC
     *(I1NMC)
              END DO
      
	        DO I1NMC=1,NNMC
	            XNNMC(I1NMC)=0.
	            XFNMC(I1NMC)=0.
	        END DO

	        DO I1NMC=1,NNMC
	            WRITE(55,1001)(ANMC(I1NMC,J1NMC),J1NMC=1,NNMC)
              END DO
                  
              DO I1NMC=1,NNMC
	            WRITE(551,1001)(A1NMC(I1NMC))
              END DO
                  
              DO I1NMC=1,NNMC
	            WRITE(552,1001)(A2NMC(I1NMC))
              END DO
                  
              DO I1NMC=1,NNMC
	            WRITE(553,1001)(A3NMC(I1NMC,J1NMC),J1NMC=1,NNMC)
	        END DO
                  
	        DO I1NMC=1,NNMC
	            WRITE(555,1001)(DIF1NMC(I1NMC,J1NMC),J1NMC=1,NNMC)
              END DO
                  
	        DO I1NMC=1,NNMC
	            WRITE(66,*)BNMC(I1NMC)
              END DO
                  
	        DO I1NMC=1,NNMC
	            WRITE(666,*)SSNMC(I1NMC,KNMC)
              END DO
	        !END IF

	        CALL PCG(ANMC,XNNMC,BNMC,XFNMC,NNMC,F22NMC)

              DO I1NMC=1,NNMC
	            SSSNMC=0.
	            DO J1NMC=1,NNMC
	                SSSNMC=DIF1NMC(I1NMC,J1NMC)*XFNMC(J1NMC)+SSSNMC
	            ENDDO
                  SFNMC(I1NMC,KNMC)=SSNMC(I1NMC,KNMC)+SSSNMC
	        END DO

              !IF(INMC.EQ.3)THEN
	        DO I1NMC=1,NNMC
                  WRITE(111,*) SFNMC(I1NMC,KNMC)
	        END DO
              !END IF
          
	    END DO

          DIIFNMC=0.
	    DO I1NMC=1,NNMC
	    DO KNMC=1,KBM1NMC
	        DIIF1NMC=ABS(SSNMC(I1NMC,KNMC)-SFNMC(I1NMC,KNMC))
	        DIIFNMC=MAX(DIIFNMC,DIIF1NMC)
              S(IINMC(I1NMC),KNMC)=SFNMC(I1NMC,KNMC)
	    END DO
          ENDDO
      
C          PRINT*,DIIFNMC
	
	    DEALLOCATE (HNMC)
	    DEALLOCATE (XNNMC)
          DEALLOCATE (XFNMC)
	    DEALLOCATE (SSNMC)
	    DEALLOCATE (SFNMC)
          DEALLOCATE (DIF1NMC)
	    DEALLOCATE (ANMC)
	    DEALLOCATE (A1NMC)
	    DEALLOCATE (A2NMC)
	    DEALLOCATE (A3NMC)
	    DEALLOCATE (AINMC)
	    DEALLOCATE (BNMC)
	    DEALLOCATE (B1NMC)
	    DEALLOCATE (IINMC)
	    DEALLOCATE (JJNMC)

          CLOSE(55)
          CLOSE(551)
          CLOSE(552)
          CLOSE(553)
          CLOSE(555)
          CLOSE(66)
          CLOSE(666)

	END DO

      DO I1NMC=1,N_CTRD
      DO KNMC=1,KBM1NMC
          IF(S(I1NMC,KNMC).LT.0)THEN
              S(I1NMC,KNMC)=0.
          ENDIF
          IF(S(I1NMC,KNMC).GE.35)THEN
              S(I1NMC,KNMC)=35.
          ENDIF
      ENDDO
      WRITE(200,'(100G)')I1NMC,(S(I1NMC,KNMC),KNMC=1,KBM1NMC)
      ENDDO

      CLOSE(200)
1000  FORMAT (I6.6)
1001  FORMAT (15555G)
      END SUBROUTINE NO_INVERSION
      


   
      SUBROUTINE PCG(ANMC,XNNMC,BNMC,XFNMC,NNMC,F22NMC)
C	IMPLICIT NONE
	INTEGER NNMC,LOOP,INMC,MM,JNMC
	DIMENSION ANMC(NNMC,NNMC)
	DIMENSION XNNMC(NNMC)
	DIMENSION BNMC(NNMC)
	DIMENSION XFNMC(NNMC)
      REAL F2NMC
	DIMENSION R(NNMC)
	DIMENSION RF(NNMC)
	DIMENSION P1(NNMC)
	DIMENSION P(NNMC)
	DIMENSION PF(NNMC)
	REAL E
	REAL U,F22NMC,E1,E2,E3,RF1,AA,U1,U2

      R=0.
	RF=0.
	P1=0.
	P=0.
      PF=0.

      !!INITIAL
	DO INMC=1,NNMC
	    AA=0.
	    DO JNMC=1,NNMC
	        AA=ANMC(INMC,JNMC)*XNNMC(JNMC)+AA
	    END DO
	    R(INMC)=BNMC(INMC)-AA
	ENDDO

      DO INMC=1,NNMC
	    P(INMC)=R(INMC)
	END DO

	F2NMC=0.

	DO INMC=1,NNMC
	    F2NMC=MAX(ABS(R(INMC)),F2NMC)
	END DO


	!!!!!!!!!!
	DO MM=1,200
	    F2NMC=0.
	    DO INMC=1,NNMC
	        F2NMC=MAX(ABS(R(INMC)),F2NMC)
	    END DO

	    IF(F2NMC.GE.F22NMC)THEN
              E1=0.
              DO INMC=1,NNMC
	            E1=R(INMC)*R(INMC)+E1
	        ENDDO

	        DO INMC=1,NNMC
	            E2=0.
	            DO JNMC=1,NNMC
	                E2=P(JNMC)*ANMC(INMC,JNMC)+E2
	            END DO
	            P1(INMC)=E2
	        END DO

	        E3=0.
              DO INMC=1,NNMC
	            E3=P1(INMC)*P(INMC)+E3
	        ENDDO
	        E=E1/E3
              
	        DO INMC=1,NNMC
	            XFNMC(INMC)=XNNMC(INMC)+E*P(INMC)
	        END DO

	        DO INMC=1,NNMC
                  RF1=0.
	            DO JNMC=1,NNMC
	                RF1=E*ANMC(INMC,JNMC)*P(JNMC)+RF1
	            END DO
	            RF(INMC)=R(INMC)-RF1
	        END DO

	        U1=0.
	        U2=0.
	        DO INMC=1,NNMC
	            U1=RF(INMC)*RF(INMC)+U1
	            U2=R(INMC)*R(INMC)+U2
	        END DO
	        U=U1/U2

	        DO INMC=1,NNMC
	            PF(INMC)=RF(INMC)+U*P(INMC)
	        END DO
	
	        DO INMC=1,NNMC
	            R(INMC)=RF(INMC)
	            P(INMC)=PF(INMC)
	            XNNMC(INMC)=XFNMC(INMC)
	        END DO
	    END IF 
	END DO
      
      END SUBROUTINE
