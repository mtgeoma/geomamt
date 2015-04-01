      SUBROUTINE G13BAZ(Y,NY,MR,NMR,IB,W0,PAR,NPAR,CY,WA,NWA,B,NB,
     *                  IERROR)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11 REVISED. IER-443 (FEB 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G13BAZ SUPERVISES THE CALCULATIONS FOR G13BAF AND G13BBF
C
C     Y     - TIME SERIES TO BE FILTERED
C     NY    - LENGTH OF TIME SERIES Y
C     MR    - ORDERS ARRAY
C     IB    - DELAY TIME,-1 FOR ARIMA,0,1... FOR TF
C     W0    - OMEGA ZERO PARAMETER FOR TF FILTERING
C     PAR   - PARAMETER ARRAY
C     NPAR  - LENGTH OF PARAMETER ARRAY
C     WA    - WORKING ARRAY
C     NWA   - LENGTH OF WORKING ARRAY
C     B     - FILTERED SERIES AND WORKING ARRAY
C     NB    - LENGTH OF ARRAY B
C     IERROR- ERROR INDICATOR
C
C     SET UP SOME USEFUL INTEGER CONSTANTS
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CY, W0
      INTEGER           IB, IERROR, NB, NMR, NPAR, NWA, NY
C     .. Array Arguments ..
      DOUBLE PRECISION  B(NB), PAR(NPAR), WA(NWA), Y(NY)
      INTEGER           MR(NMR)
C     .. Local Scalars ..
      DOUBLE PRECISION  A0, B0, PHIDY1, RES1, RES2, RK, X, Z
      INTEGER           I, IDPLDS, IFAIL1, IPDD, IPDYD, IPYD, IQD,
     *                  IQPDYD, J, K, K1, K2, K3, K4, K5, L1, L2,
     *                  LENGTH, MAXQPD, ME, NBB, NFA
C     .. Local Arrays ..
      INTEGER           N(8)
C     .. External Functions ..
      DOUBLE PRECISION  G13BAX
      EXTERNAL          G13BAX
C     .. External Subroutines ..
      EXTERNAL          F04ARF, X03AAF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MOD
C     .. Executable Statements ..
      IQD = MR(3) + MR(6)*MR(7)
      IDPLDS = MR(2) + MR(5)
      IPDD = MR(1) + MR(2) + (MR(4)+MR(5))*MR(7)
      MAXQPD = MAX(IQD,IPDD)
      NBB = NY
      LENGTH = NY
      N(1) = 0
      N(2) = MR(1)
      N(3) = N(2) + MR(3)
      N(4) = N(3) + MR(4)
      IF (NMR.EQ.7) GO TO 20
      IPYD = MR(8) + MR(11)*MR(14)
      IPDYD = IPYD + MR(9) + MR(12)*MR(14)
      NBB = NBB + MAXQPD
      LENGTH = LENGTH + IPDD
      IQPDYD = IQD + IPDYD
      ME = IQPDYD*IQPDYD
      N(5) = N(4) + MR(6)
      N(6) = N(5) + MR(8)
      N(7) = N(6) + MR(10)
      N(8) = N(7) + MR(11)
   20 CONTINUE
C
C     SHIFT SERIES INTO B LEAVING SPACE FOR ADDITIONAL BACKFORECASTS
C
      K1 = 0
      IF (NMR.EQ.14) K1 = MAXQPD
      DO 40 I = 1, NY
         K1 = K1 + 1
         B(K1) = Y(I)
   40 CONTINUE
      IF (NMR.EQ.7) GO TO 240
C
C     FOR BACKFORECASTING GENERATE THE CONSTANT
C     PHI-DASHED-Y OF 1, AND THE B COEFFICIENTS
C
      K1 = ME + IQPDYD
      IF (K1.EQ.0) GO TO 80
      DO 60 I = 1, K1
         WA(I) = 0.0D0
   60 CONTINUE
   80 CONTINUE
      L1 = N(5)
      K1 = MR(8)
      L2 = N(7)
      K2 = MR(11)
      X = 1.0D0
      PHIDY1 = G13BAX(X,PAR,NPAR,L1,K1)*G13BAX(X,PAR,NPAR,L2,K2)
      K1 = K1 + 1
      K2 = K2 + 1
      B0 = 0.0D0
      DO 120 I = 1, K1
         X = -PAR(L1)
         IF (I.EQ.1) X = 1.0D0
         K4 = L2
         DO 100 J = 1, K2
            Z = -PAR(K4)
            IF (J.EQ.1) Z = 1.0D0
            K3 = ME + (I-1) + (J-1)*MR(14)
            Z = Z*X
            IF (K3.EQ.ME) B0 = B0 + Z
            IF (K3.GT.ME) WA(K3) = WA(K3) + Z
            K4 = K4 + 1
  100    CONTINUE
         L1 = L1 + 1
  120 CONTINUE
      K1 = MR(9) + MR(12)
      IF (K1.EQ.0) GO TO 180
      K2 = ME + IPYD
      DO 160 I = 1, K1
         K3 = 1
         IF (I.GT.MR(9)) K3 = MR(14)
         K2 = K2 + K3
         K4 = K2
         DO 140 J = 1, K2
            K5 = K4 - K3
            IF (K5.EQ.ME) WA(K4) = WA(K4) - B0
            IF (K5.GT.ME) WA(K4) = WA(K4) - WA(K5)
            K4 = K4 - 1
  140    CONTINUE
  160 CONTINUE
C
C     BACKFORECAST THE REQUIRED VALUES
C
  180 IF (MOD(K1,2).NE.0) CY = -CY
      IF (IPDD.EQ.0) GO TO 240
      K1 = MAXQPD
      DO 220 I = 1, IPDD
         B(K1) = PHIDY1*CY
         IF (IPDYD.EQ.0) GO TO 200
         K2 = K1 + 1
         K3 = ME + 1
         IF (K2+IPDYD-1.GT.NBB) GO TO 1080
         NFA = 0
         CALL X03AAF(B(K2),IPDYD,WA(K3),IPDYD,IPDYD,1,1,0.0D0,0.0D0,
     *               RES1,RES2,.FALSE.,NFA)
         B(K1) = B(K1) - RES1
  200    CONTINUE
         K1 = K1 - 1
  220 CONTINUE
  240 CONTINUE
C
C     DIFFERENCE THE SERIES
C
      IF (IDPLDS.EQ.0) GO TO 300
      DO 280 I = 1, IDPLDS
         K3 = 1
         IF (I.GT.MR(2)) K3 = MR(7)
         LENGTH = LENGTH - K3
         K4 = NBB
         IF (LENGTH.LT.1) GO TO 1080
         DO 260 J = 1, LENGTH
            K2 = K4 - K3
            B(K4) = B(K4) - B(K2)
            K4 = K4 - 1
  260    CONTINUE
  280 CONTINUE
  300 CONTINUE
C
C     APPLY SEASONAL THEN NON-SEASONAL AR FILTER
C
      IF (IB.GE.0) GO TO 380
      DO 360 I = 1, 2
         L1 = N(1)
         IF (I.EQ.1) L1 = N(3)
         K1 = 1
         IF (I.EQ.1) K1 = 4
         K1 = MR(K1)
         IF (K1.EQ.0) GO TO 360
         K2 = 1
         IF (I.EQ.1) K2 = MR(7)
         LENGTH = LENGTH - K1*K2
         IF (LENGTH.LT.1) GO TO 1080
         K4 = NBB
         DO 340 J = 1, LENGTH
            L2 = L1
            K3 = K4
            DO 320 K = 1, K1
               K3 = K3 - K2
               L2 = L2 + 1
               B(K4) = B(K4) - PAR(L2)*B(K3)
  320       CONTINUE
            K4 = K4 - 1
  340    CONTINUE
  360 CONTINUE
      GO TO 460
  380 CONTINUE
C
C     FUDGED CODE FOR AR-LIKE PART OF TF FILTERING
C
      L1 = N(1) + IB
      K1 = MR(1)
      LENGTH = LENGTH - K1
      K1 = K1 - IB
      IF (LENGTH.LT.1) GO TO 1080
      K4 = NBB
      DO 440 J = 1, LENGTH
         K3 = K4 - IB
         B(K4) = W0*B(K3)
         IF (K1.EQ.0) GO TO 420
         L2 = L1
         DO 400 K = 1, K1
            K3 = K3 - 1
            L2 = L2 + 1
            B(K4) = B(K4) - PAR(L2)*B(K3)
  400    CONTINUE
  420    K4 = K4 - 1
  440 CONTINUE
  460 CONTINUE
      IF (NMR.EQ.7) GO TO 780
      IF (IQD.EQ.0) GO TO 960
C
C     GENERATE THE A COEFFICIENTS
C
      IF (IPDYD.EQ.0) GO TO 600
      L1 = N(2)
      L2 = N(4)
      K1 = MR(3) + 1
      K2 = MR(6) + 1
      K4 = ME + IPDYD
      A0 = 0.0D0
      DO 500 I = 1, K1
         X = 1.0D0
         IF (I.GT.1) X = -PAR(L1)
         K5 = L2
         DO 480 J = 1, K2
            Z = 1.0D0
            IF (J.GT.1) Z = -PAR(K5)
            K3 = K4 + (I-1) + (J-1)*MR(7)
            Z = X*Z
            IF (K3.EQ.K4) A0 = A0 + Z
            IF (K3.GT.K4) WA(K3) = WA(K3) + Z
            K5 = K5 + 1
  480    CONTINUE
         L1 = L1 + 1
  500 CONTINUE
C
C     FORM MATRIX OF A AND B COEFFICIENTS
C     IN WA(1) ONWARDS
C
      K3 = ME
      K4 = 1
      K1 = IPDYD + 1
      DO 540 I = 1, K1
         X = B0
         IF (I.GT.1) X = WA(K3)
         K2 = K4
         DO 520 J = 1, IQD
            WA(K2) = X
            K2 = K2 + IQPDYD + 1
  520    CONTINUE
         K4 = K4 + IQPDYD
         K3 = K3 + 1
  540 CONTINUE
      K3 = ME + IQPDYD
      K1 = IQD + 1
      K4 = K1
      DO 580 I = 1, K1
         X = WA(K3)
         IF (I.EQ.K1) X = A0
         K2 = K4
         DO 560 J = 1, IPDYD
            WA(K2) = X
            K2 = K2 + IQPDYD + 1
  560    CONTINUE
         K4 = K4 + IQPDYD
         K3 = K3 - 1
  580 CONTINUE
  600 CONTINUE
C
C     CALCULATE PHI-DASHED OF 1 AND THETA-DASHED OF 1
C     AND HENCE CONSTANT K.
C
      IF (IDPLDS.GT.0) GO TO 620
      L1 = N(1)
      K1 = MR(1)
      L2 = N(3)
      K2 = MR(4)
      RK = 1.0D0
      X = G13BAX(W0,PAR,NPAR,L1,K1)*G13BAX(RK,PAR,NPAR,L2,K2)
      L1 = N(2)
      K1 = MR(3)
      L2 = N(4)
      K2 = MR(6)
      Z = G13BAX(RK,PAR,NPAR,L1,K1)*G13BAX(RK,PAR,NPAR,L2,K2)
      RK = (X/Z)*PHIDY1*CY
      GO TO 640
  620 RK = 0.0D0
  640 CONTINUE
C
C     PUT K AND ZS INTO RHS OF EQUATION
C
      K1 = ME
      DO 660 I = 1, IQD
         K1 = K1 + 1
         WA(K1) = RK
  660 CONTINUE
      IF (IPDYD.EQ.0) GO TO 700
      IF (IPDYD.GT.LENGTH) GO TO 1080
      K2 = MAXQPD
      DO 680 I = 1, IPDYD
         K1 = K1 + 1
         K2 = K2 + 1
         WA(K1) = B(K2)
  680 CONTINUE
C
C     SOLVE EQUATIONS FOR UNKNOWN BETAS
C
      K1 = ME + 1
      K2 = ME + IQPDYD + 1
      IFAIL1 = 1
      CALL F04ARF(WA(1),IQPDYD,WA(K1),IQPDYD,WA(K1),WA(K2),IFAIL1)
      IF (IFAIL1.NE.0) GO TO 1100
C
C     MOVE INITIAL BETAS INTO B ARRAY
C
  700 CONTINUE
      K1 = ME + 1
      K2 = MAXQPD - IQD
      DO 720 I = 1, IQD
         K2 = K2 + 1
         B(K2) = WA(K1)
         K1 = K1 + 1
  720 CONTINUE
C
C     CALCULATE INITIAL GAMMAS
C
      K1 = MR(6)*MR(7)
      IF (K1.EQ.0) GO TO 820
      K4 = MR(3)
      IF (K4.EQ.0) GO TO 820
      K2 = MAXQPD
      L1 = N(2)
      DO 760 I = 1, K1
         K3 = K2
         L2 = L1
         DO 740 J = 1, K4
            K3 = K3 - 1
            L2 = L2 + 1
            B(K2) = B(K2) - B(K3)*PAR(L2)
  740    CONTINUE
         K2 = K2 - 1
  760 CONTINUE
      GO TO 820
C
C     ZERO START OF B ARRAY IF NO ARIMA FOR Y
C
  780 CONTINUE
      K1 = NBB - LENGTH
      IF (K1.EQ.0) GO TO 820
      DO 800 I = 1, K1
         B(I) = 0.0D0
  800 CONTINUE
C
C     APPLY SEASONAL MA FILTER
C
  820 CONTINUE
      K5 = MR(6)
      IF (K5.EQ.0) GO TO 900
      K1 = NBB - LENGTH
      K2 = K1
      L1 = N(4)
      DO 860 I = 1, LENGTH
         K2 = K2 + 1
         K4 = K2
         L2 = L1
         DO 840 J = 1, K5
            K4 = K4 - MR(7)
            X = 0.0D0
            IF (K4.GE.1) X = B(K4)
            L2 = L2 + 1
            B(K2) = B(K2) + PAR(L2)*X
  840    CONTINUE
  860 CONTINUE
      IF (NMR.EQ.7) GO TO 900
C
C     AGAIN MOVE INITIAL BETAS INTO B ARRAY
C     OVERWRITING INITIAL GAMMAS
C
      K1 = MR(3)
      IF (K1.EQ.0) GO TO 900
      K2 = MAXQPD
      K3 = ME + IQD
      DO 880 I = 1, K1
         B(K2) = WA(K3)
         K2 = K2 - 1
         K3 = K3 - 1
  880 CONTINUE
C
C     APPLY NON SEASONAL MA FILTER
C
  900 CONTINUE
      K5 = MR(3)
      IF (K5.EQ.0) GO TO 960
      K1 = NBB - LENGTH
      L1 = N(2)
      K2 = K1
      DO 940 I = 1, LENGTH
         K2 = K2 + 1
         K4 = K2
         L2 = L1
         DO 920 J = 1, K5
            K4 = K4 - 1
            X = 0.0D0
            IF (K4.GE.1) X = B(K4)
            L2 = L2 + 1
            B(K2) = B(K2) + PAR(L2)*X
  920    CONTINUE
  940 CONTINUE
  960 CONTINUE
      IF (NMR.EQ.7) GO TO 1020
C
C     MOVE UP FILTERED SERIES IF ARIMA FOR Y
C
      K1 = MAXQPD
      DO 980 I = 1, LENGTH
         K1 = K1 + 1
         B(I) = B(K1)
  980 CONTINUE
      K2 = MR(9) + MR(12)
      IF (MOD(K2,2).NE.0) CY = -CY
C
C     PAD B WITH ZEROS
C
      K1 = LENGTH + 1
      IF (K1.GT.NBB) GO TO 1020
      DO 1000 I = K1, NBB
         B(I) = 0.0D0
 1000 CONTINUE
 1020 CONTINUE
      K1 = NBB + 1
      IF (K1.GT.NB) GO TO 1060
      DO 1040 I = K1, NB
         B(I) = 0.0D0
 1040 CONTINUE
 1060 CONTINUE
C
C     SUCCESSFUL RETURN
      RETURN
C
C     UNSUCCESSFUL RETURN
C     SERIES TOO SHORT
 1080 IERROR = 4
      GO TO 1120
C     MATRIX SOLN FOR INITIAL BETAS IS SINGULAR
 1100 IERROR = 9
 1120 RETURN
      END
