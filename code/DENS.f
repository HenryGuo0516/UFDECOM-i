
      MODULE DENSFMA
      CONTAINS
      
      PURE REAL FUNCTION ATG(S4,T4,P4)
!$omp declare target
      REAL, INTENT(IN) :: S4,T4,P4
      DS = S4 - 35.0
      ATG = (((-2.1687E-16*T4+1.8676E-14)*T4-4.6206E-13)*P4
     * +((2.7759E-12*T4-1.1351E-10)*DS+((-5.4481E-14*T4
     * +8.733E-12)*T4-6.7795E-10)*T4+1.8741E-8))*P4
     * +(-4.2393E-8*T4+1.8932E-6)*DS
     * +((6.6228E-10*T4-6.836E-8)*T4+8.5258E-6)*T4+3.5803E-5
      RETURN
      END FUNCTION ATG


      PURE REAL FUNCTION THETA(S4,T04,P04,PR)
!$omp declare target
      REAL, INTENT(IN) :: S4,T04,P04,PR
      P4=P04
      T4=T04
      H4 = PR - P4
      XK = H4*ATG(S4,T4,P4)
      T4 = T4 + 0.5*XK
      Q4 = XK
      P4 = P4 + 0.5*H4
      XK = H4*ATG(S4,T4,P4)
      T4 = T4 + 0.29289322*(XK-Q4)
      Q4 = 0.58578644*XK + 0.121320344*Q4
      XK = H4*ATG(S4,T4,P4)
      T4 = T4 + 1.707106781*(XK-Q4)
      Q4 = 3.414213562*XK - 4.121320344*Q4
      P4 = P4 + 0.5*H4
      XK = H4*ATG(S4,T4,P4)
      THETA = T4 + (XK-2.0*Q4)/6.0
      RETURN
      END FUNCTION THETA
      
      
c      REAL FUNCTION SVAN(S4,T4,P04,SIGMA)
      PURE REAL FUNCTION SIGMA(S4,T4,P04) !CBR: SVAN is not used here
!$omp declare target
      parameter R3500=1028.1063,RR4=4.8314E-4
      parameter DR350=28.106331
      REAL, INTENT(IN) :: S4,T4,P04
      REAL P4,SIG,SR,RR1,RR2,RR3
      REAL A4,B4,C4,D4,E4,AA1,BB1,AW,BW,KK,K0,KW,K35
      
      P4 = P04/10.
      SR = SQRT(ABS(S4))
      RR1 = ((((6.536332E-9*T4-1.120083E-6)*T4+1.001685E-4)*T4
     *   -9.095290E-3)*T4+6.793952E-2)*T4-28.263737
      RR2 = (((5.3875E-9*T4-8.2467E-7)*T4+7.6438E-5)*T4-4.0899E-3)*T4
     * +8.24493E-1
      RR3 = (-1.6546E-6*T4+1.0227E-4)*T4-5.72466E-3
      SIG = (RR4*S4 + RR3*SR + RR2)*S4 + RR1
      V350P = 1.0/R3500
      SVA = -SIG*V350P/(R3500+SIG)
      SIGMA = SIG+DR350
c      SVAN = SVA*1.0E+8 !CBR: SVAN is not used here
      if (P4/=0) then
      
      E4 = (9.1697E-10*T4+2.0816E-8)*T4-9.9348E-7
      BW = (5.2787E-8*T4-6.12293E-6)*T4+3.47718E-5
      B4 = BW + E4*S4
      D4 = 1.91075E-4
      C4 = (-1.6078E-6*T4-1.0981E-5)*T4+2.2838E-3
      AW = ((-5.77905E-7*T4+1.16092E-4)*T4+1.43713E-3)*T4
     *   -0.1194975
      A4 = (D4*SR + C4)*S4 + AW
      BB1 = (-5.3009E-4*T4+1.6483E-2)*T4+7.944E-2
      AA1 = ((-6.1670E-5*T4+1.09987E-2)*T4-0.603459)*T4+54.6746
      KW = (((-5.155288E-5*T4+1.360477E-2)*T4-2.327105)*T4
     *   +148.4206)*T4-1930.06
      K0 = (BB1*SR + AA1)*S4 + KW
      DK = (B4*P4 + A4)*P4 + K0
      K35  = (5.03217E-5*P4+3.359406)*P4+21582.27
      GAM = P4/K35
      PK = 1.0 - GAM
      SVA = SVA*PK + (V350P+SVA)*P4*DK/(K35*(K35+DK))

c      SVAN = SVA*1.0E+8 !CBR: SVAN is not used here
      V350P = V350P*PK
      DR35P = GAM/V350P
      DVAN = SVA/(V350P*(V350P+SVA))
      SIGMA = DR350+DR35P-DVAN
      
      endif
      RETURN
      END FUNCTION SIGMA
      
      END MODULE DENSFMA

! =================== DENS CPU version ===================      
      
      SUBROUTINE DENS
      USE MOD_GLOBAL
      use DENSFMA
      PR = 0.0
c      RHO=0.
      DO K = 1,KBM1
      DO I = 1,N_CTRD
      IF(FSM(I).GT.0.0) THEN   !WUHUI
          RZU = -GRAV*1.025*(ZZ(K)*H(I)+ELF(I)*(ZZ(K)+1.))*0.01
          PT = THETA(S(I,K),25.0,RZU,PR)
c          SVA = SVAN(S(I,K),PT,RZU,SIGMA)
c          RHO(I,K) = SIGMA*1.E-3
          RHO(I,K)=SIGMA(S(I,K),PT,RZU)*1.E-3
#ifdef MODULE_SED
          IF (THOUR.GE.SED_BEG) THEN
              SEDMT = AMIN1(5.0,SED(I,K))*1.E-3
              RHO(I,K) = RHO(I,K)+(1.-(RHO(I,K)+1.025)/2.65)*SEDMT !LXY 
          ENDIF
#endif
      ELSE
          RHO(I,K) = 0.0
      ENDIF
      ENDDO
      ENDDO
      RETURN
      END SUBROUTINE DENS
      
! =================== DENS GPU version ===================      
      
      SUBROUTINE DENS_GPU
      USE MOD_GLOBAL
      use DENSFMA
      PR = 0.0
c      RHO=0.
      
      
#ifdef MODULE_SED
      IF (THOUR.GE.SED_BEG) THEN
#endif

c!$OMP TARGET DEFAULTMAP(present: allocatable)
c!$OMP TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
c      DO K = 1,KBM1
c      DO I = 1,N_CTRD
      DO concurrent (I=1:N_CTRD,K=1:KBM1) local(RZU,PT) local_init(PR)
      IF(FSM(I).GT.0.0) THEN   !WUHUI
          RZU = -GRAV*1.025*(ZZ(K)*H(I)+ELF(I)*(ZZ(K)+1.))*0.01
          PT = THETA(S(I,K),25.0,RZU,PR)
c          SVA = SVAN(S(I,K),PT,RZU,SIGMA)
c          RHO(I,K) = SIGMA*1.E-3
          RHO(I,K)=SIGMA(S(I,K),PT,RZU)*1.E-3
#ifdef MODULE_SED
          RHO(I,K) = RHO(I,K)+
     *+(1.-(RHO(I,K)+1.025)/2.65)*AMIN1(5.0,SED(I,K))*1.E-3 !LXY
c          SEDMT = AMIN1(5.0,SED(I,K))*1.E-3
c          RHO(I,K) = RHO(I,K)+(1.-(RHO(I,K)+1.025)/2.65)*SEDMT !LXY 
#endif
      ELSE
          RHO(I,K) = 0.0
      ENDIF
      ENDDO
c      ENDDO
c!$OMP END TARGET 
      
#ifdef MODULE_SED
      ENDIF
#endif

      
      RETURN
      END SUBROUTINE DENS_GPU
      
      