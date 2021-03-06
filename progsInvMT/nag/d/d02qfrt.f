      SUBROUTINE D02QFR(X,Y,XOUT,YOUT,YPOUT,NINT,NEQN,KOLD,PHI,IVC,IV,
     *                  KGI,GI,ALPHA,OG,OW,OX,OY)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C
C     MODIFIED TO INTERPOLATE FIRST NINT COMPONENTS
C                             ----------
C
C
C     PURPOSE
C     =======
C            APPROXIMATES THE SOLUTION AT XOUT BY EVALUATING THE
C            POLYNOMIAL COMPUTED IN D02QFQ AT XOUT.  MUST BE USED IN
C            CONJUNCTION WITH D02QFQ.
C     DESCRIPTION
C     ===========
C
C     WRITTEN BY L. F. SHAMPINE AND M. K. GORDON
C
C     ABSTRACT
C
C
C     THE METHODS IN SUBROUTINE D02QFQ  APPROXIMATE THE SOLUTION NEAR  X
C     BY A POLYNOMIAL.  SUBROUTINE  D02QFR  APPROXIMATES THE SOLUTION AT
C     XOUT BY EVALUATING THE POLYNOMIAL THERE. INFORMATION DEFINING THIS
C     POLYNOMIAL IS PASSED FROM D02QFQ  SO D02QFR  CANNOT BE USED ALONE.
C
C     THIS CODE IS COMPLETELY EXPLAINED AND DOCUMENTED IN THE TEXT,
C     COMPUTER SOLUTION OF ORDINARY DIFFERENTIAL EQUATIONS, THE INITIAL
C     VALUE PROBLEM  BY L. F. SHAMPINE AND M. K. GORDON.
C     FURTHER DETAILS ON USE OF THIS CODE ARE AVAILABLE IN *SOLVING
C     ORDINARY DIFFERENTIAL EQUATIONS WITH ODE, STEP, AND INTRP*,
C     BY L. F. SHAMPINE AND M. K. GORDON, SLA-73-1060.
C
C     INPUT TO D02QFR --
C
C     THE USER PROVIDES STORAGE IN THE CALLING PROGRAM FOR THE ARRAYS IN
C     THE CALL LIST
C      DIMENSION Y(NEQN),YOUT(NEQN),YPOUT(NEQN),PHI(NEQN,16),OY(NEQN)
C                AND ALPHA(12),OG(13),OW(12),GI(11),IV(10)
C     AND DEFINES
C      XOUT -- POINT AT WHICH SOLUTION IS DESIRED.
C     THE REMAINING PARAMETERS ARE DEFINED IN  D02QFQ  AND PASSED TO
C     D02QFR  FROM THAT SUBROUTINE
C
C     OUTPUT FROM  D02QFR  --
C
C      YOUT(*) -- SOLUTION AT  XOUT
C      YPOUT(*) -- DERIVATIVE OF SOLUTION AT  XOUT
C
C     THE REMAINING PARAMETERS ARE RETURNED UNALTERED FROM THEIR INPUT
C     VALUES.  INTEGRATION WITH  D02QFQ  MAY BE CONTINUED.
C
C     REFERENCES
C     ==========
C               SHAMPINE L.F., GORDON M.K., *SOLVING ORDINARY
C                 DIFFERENTIAL EQUATIONS WITH ODE, STEP, AND INTRP*,
C                 SLA-73-1060, SANDIA LABORATORIES, 1973.
C               WATTS H.A., SHAMPINE L.F., *A SMOOTHER INTERPOLANT FOR
C                 DE/STEP,INTRP : II*, SAND84-0293, SANDIA LABORATORIES,
C                 1984.
C
C     DATE WRITTEN   740101   (YYMMDD)
C     REVISION DATE  840201   (YYMMDD)
C     AUTHOR  SHAMPINE, L.F.,  SNLA
C           GORDON, M.K.
C             MODIFIED BY H.A. WATTS
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  OX, X, XOUT
      INTEGER           IVC, KGI, KOLD, NEQN, NINT
C     .. Array Arguments ..
      DOUBLE PRECISION  ALPHA(12), GI(11), OG(13), OW(12), OY(NEQN),
     *                  PHI(NEQN,16), Y(NEQN), YOUT(NINT), YPOUT(NINT)
      INTEGER           IV(10)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALP, GAMMA, GDI, GDIF, H, HI, HMU, RMU, SIGMA,
     *                  TEMP1, TEMP2, TEMP3, XI, XIM1, XIQ
      INTEGER           I, IQ, IW, J, JQ, KP1, KP2, L, M
C     .. Local Arrays ..
      DOUBLE PRECISION  C(13), G(13), W(13)
C     .. Executable Statements ..
C
      KP1 = KOLD + 1
      KP2 = KOLD + 2
C
      HI = XOUT - OX
      H = X - OX
      XI = HI/H
      XIM1 = XI - 1.D0
C
C     INITIALIZE W(*) FOR COMPUTING G(*)
C
      XIQ = XI
      DO 20 IQ = 1, KP1
         XIQ = XI*XIQ
         TEMP1 = IQ*(IQ+1)
         W(IQ) = XIQ/TEMP1
   20 CONTINUE
C
C     COMPUTE THE DOUBLE INTEGRAL TERM GDI
C
      IF (KOLD.LE.KGI) GO TO 100
      IF (IVC.GT.0) GO TO 40
      GDI = 1.0D0/TEMP1
      M = 2
      GO TO 60
   40 IW = IV(IVC)
      GDI = OW(IW)
      M = KOLD - IW + 3
   60 IF (M.GT.KOLD) GO TO 120
      DO 80 I = M, KOLD
         GDI = OW(KP2-I) - ALPHA(I)*GDI
   80 CONTINUE
      GO TO 120
  100 GDI = GI(KOLD)
C
C     COMPUTE G(*) AND C(*)
C
  120 G(1) = XI
      G(2) = 0.5D0*XI*XI
      C(1) = 1.0D0
      C(2) = XI
      IF (KOLD.LT.2) GO TO 180
      DO 160 I = 2, KOLD
         ALP = ALPHA(I)
         GAMMA = 1.0D0 + XIM1*ALP
         L = KP2 - I
         DO 140 JQ = 1, L
            W(JQ) = GAMMA*W(JQ) - ALP*W(JQ+1)
  140    CONTINUE
         G(I+1) = W(1)
         C(I+1) = GAMMA*C(I)
  160 CONTINUE
C
C     DEFINE INTERPOLATION PARAMETERS
C
  180 SIGMA = (W(2)-XIM1*W(1))/GDI
      RMU = XIM1*C(KP1)/GDI
      HMU = RMU/H
C
C     INTERPOLATE FOR THE SOLUTION -- YOUT
C     AND FOR THE DERIVATIVE OF THE SOLUTION -- YPOUT
C
C     DO 200 L = 1, NEQN
      DO 200 L = 1, NINT
         YOUT(L) = 0.0D0
         YPOUT(L) = 0.0D0
  200 CONTINUE
      DO 240 J = 1, KOLD
         I = KP2 - J
         GDIF = OG(I) - OG(I-1)
         TEMP2 = (G(I)-G(I-1)) - SIGMA*GDIF
         TEMP3 = (C(I)-C(I-1)) + RMU*GDIF
C        DO 220 L = 1, NEQN
         DO 220 L = 1, NINT
            YOUT(L) = YOUT(L) + TEMP2*PHI(L,I)
            YPOUT(L) = YPOUT(L) + TEMP3*PHI(L,I)
  220    CONTINUE
  240 CONTINUE
C     DO 260 L = 1, NEQN
      DO 260 L = 1, NINT
         YOUT(L) = ((1.0D0-SIGMA)*OY(L)+SIGMA*Y(L)) + H*(YOUT(L)+(G(1)
     *             -SIGMA*OG(1))*PHI(L,1))
         YPOUT(L) = HMU*(OY(L)-Y(L)) + (YPOUT(L)+(C(1)+RMU*OG(1))*PHI(L,
     *              1))
  260 CONTINUE
C
      RETURN
C
C
C     END OF D02QFR (SINTRP)
C
C
      END
