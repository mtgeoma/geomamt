      SUBROUTINE D02SAV(P,M,A,B,N,N1,PE,PF,E,DP,NPOINT,WP,IWP,X,X1,
     *                  ICOUNT,RANGE,RANGE1,BC,BC1,FCN,FCN1,FCN2,EQN,
     *                  CONSTR,YMAX,MONIT,PRSOL,C,M1,W,IW1,IW2,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 8F REVISED. IER-299 (APR 1981).
C     MARK 9B REVISED. IER-363 (JAN 1982)
C     MARK 11 REVISED. IER-419 (FEB 1984).
C     MARK 11D REVISED. IER-468 (NOV 1985).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     SOLVES A TWO POINT BVP
C     BC1, BC, EQN, FCN1, FCN2, FCN, MONIT, PRSOL, RANGE1, RANGE
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X, X1, YMAX
      INTEGER           ICOUNT, IFAIL, IW1, IW2, IWP, M, M1, N, N1,
     *                  NPOINT
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IW1,2), B(IW1,2), C(IW1,M1), DP(M), E(N),
     *                  P(M), PE(M), PF(M), W(IW1,IW2), WP(IWP,6)
C     .. Function Arguments ..
      LOGICAL           CONSTR
      EXTERNAL          CONSTR
C     .. Subroutine Arguments ..
      EXTERNAL          BC, BC1, EQN, FCN, FCN1, FCN2, MONIT, PRSOL,
     *                  RANGE, RANGE1
C     .. Scalars in Common ..
      DOUBLE PRECISION  DM, EPS, EPSFAC, MACHEP, PNORM, PNORM1, SQEPS
      INTEGER           COUNT, ICASE, IEPS, IFAIL1, IFLAG, ISTATE, IZZ
C     .. Arrays in Common ..
      DOUBLE PRECISION  COUT(2)
C     .. Local Scalars ..
      DOUBLE PRECISION  S
      INTEGER           I, IF1, IW3, J
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          D02SAT, D02SAU, D02SAW, D02SAX, D02SAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN, SIGN, SQRT
C     .. Common blocks ..
      COMMON            /AD02SA/PNORM, PNORM1, EPS, MACHEP, SQEPS, DM,
     *                  EPSFAC, ISTATE, IFLAG, COUNT, IFAIL1, IEPS,
     *                  ICASE
      COMMON            /BD02SA/COUT, IZZ
C     .. Executable Statements ..
      I = 3*M + 23
      J = 4*M + 12
      WP(NPOINT,2) = 0.D0
      IF (N1.GT.0 .AND. N.GE.N1 .AND. M.GE.N1 .AND. IWP.GE.NPOINT .AND.
     *    ICOUNT.GE.0 .AND. NPOINT.GE.2 .AND. IW1.GE.M .AND. IW1.GE.
     *    N .AND. IW2.GE.I .AND. IW2.GE.J .AND. YMAX.GE.0.D0 .AND.
     *    M1.GE.1) GO TO 40
C     INPUT ERROR
   20 IFAIL1 = 1
      GO TO 820
   40 DO 60 I = 1, N
         IF (E(I).LE.0.D0) GO TO 20
   60 CONTINUE
      DO 80 I = 1, M
         IF (PE(I).LE.0.D0) GO TO 20
   80 CONTINUE
      IW3 = MAX(M+15,2*M+4)
      MACHEP = X02AJF()
      SQEPS = MIN(1.0D-5,SQRT(MACHEP))
      DO 100 I = 1, M
         IF (PF(I).LE.0.0D0) PF(I) = SQEPS
  100 CONTINUE
      ISTATE = 0
      COUNT = 0
C     CHECK CONSTRAINTS
  120 IF (CONSTR(P,M)) GO TO 140
C     CONSTRAINTS VIOLATED
      IFAIL1 = 2
      GO TO 660
C     SET UP RANGE
  140 GO TO (160,180,200) ICASE
  160 CALL RANGE1(WP(1,4),NPOINT,P,M)
      GO TO 220
  180 CALL RANGE(WP(1,4),WP(2,4),P)
      NPOINT = 2
      GO TO 220
  200 WP(1,4) = X
      WP(2,4) = X1
      NPOINT = 2
      CALL D02SAT(P,A,B,IW1,0)
  220 IF (WP(1,4).EQ.WP(NPOINT,4)) GO TO 260
      IF (NPOINT.EQ.2) GO TO 280
      S = SIGN(1.D0,WP(NPOINT,4)-WP(1,4))
      DO 240 I = 2, NPOINT
         IF (S*(WP(I,4)-WP(I-1,4)).LE.0.D0) GO TO 260
  240 CONTINUE
      GO TO 280
C     RANGE INVALID
  260 IFAIL1 = 3
      GO TO 660
C     CALCULATE RESIDUAL
  280 IF1 = 1
      CALL D02SAZ(W(1,1),P,M,A,E,N,N1,WP(1,4),NPOINT,WP(1,1),WP(1,2)
     *            ,WP(1,3),0,FCN,FCN1,FCN2,EQN,BC,BC1,YMAX,W(1,2*M+9)
     *            ,IW1,IF1)
      IFAIL1 = IF1
      COUNT = COUNT + 1
      IF (IFAIL1.EQ.0) GO TO 420
C     ERROR IN D02SAZ
      GO TO (300,320,340,360,300) IFAIL1
C     IMPOSS ERROR IN D02SAZ
  300 IFAIL1 = 16
      GO TO 820
C     STEP TOO SMALL
  320 IFAIL1 = 4
      GO TO 380
C     STEP TOO SMALL INITIALLY
  340 IFAIL1 = 5
      GO TO 380
C     SOLUTION TOO LARGE
  360 IFAIL1 = 6
  380 WP(1,5) = W(1,2*M+10)
      DO 400 I = 1, N
         W(I,1) = W(I,2*M+9)
  400 CONTINUE
      IF (ISTATE.EQ.0) GO TO 820
      GO TO 660
  420 IF (ISTATE.NE.0) GO TO 680
C     CALCULATE JACOBIAN AND SVD
  440 IF1 = IFAIL1
      CALL D02SAW(W(1,8),W(1,9),W(1,M+9),W(1,1)
     *            ,P,PF,M,A,B,DP,E,N,N1,WP(1,5),WP(1,6)
     *            ,NPOINT,X,X1,WP(1,1),WP(1,2),WP(1,3)
     *            ,FCN,FCN1,FCN2,EQN,BC,BC1,RANGE,RANGE1,CONSTR,YMAX,
     *            W(1,2*M+9),IW1,IW3,IF1)
      IF (ICASE.EQ.3) CALL D02SAT(P,A,B,IW1,1)
      IFAIL1 = IF1
      IF (IFAIL1.EQ.0) GO TO 680
      GO TO (460,480,500,300,520,540,560,300,480,500,
     *       620,640) IFAIL1
C     IMPOSS EXIT FROM D02SAW
  460 IFAIL1 = 17
      GO TO 820
C     CONSTRAINTS VIOLATED IN D02SAW
  480 IFAIL1 = 7
      GO TO 660
C     RANGE VIOLATION IN D02SAW
  500 IFAIL1 = 8
      GO TO 660
C     STEP TOO SMALL IN D02SAW
  520 IFAIL1 = 9
      GO TO 580
C     INITIAL STEP TOO SMALL IN D02SAW
  540 IFAIL1 = 10
      GO TO 580
C     SOLUTION TOO LARGE
  560 IFAIL1 = 11
  580 WP(1,5) = W(1,2*M+14)
      DO 600 I = 1, N
         W(I,1) = W(1,2*M+13)
  600 CONTINUE
      WP(1,5) = W(1,2*M+14)
      GO TO 660
C     COLUMN OF JACOBIAN INSIGNIFICANT
  620 IFAIL1 = 12
      GO TO 660
C     FAIL IN F02SZF
  640 IFAIL1 = 13
  660 IF (ISTATE.EQ.0) ISTATE = 7
      ISTATE = -ISTATE
C     CALL NEWTON.
  680 CALL D02SAX(P,M,A,B,W(1,9),W(1,M+9),W(1,8),IW1,PE,PF,W(1,1),W(1,2)
     *            ,MONIT)
      IF (ISTATE.EQ.7) GO TO 820
      IF (IFLAG.GE.0) GO TO 740
      IF (IFLAG.EQ.-2) GO TO 720
C     IMPOSS ERROR IN D02SAX
  700 IFAIL1 = 18
      GO TO 820
C     ITERATION HAS NOT CONVERGED
  720 IFAIL1 = 14
      GO TO 820
  740 IF (IFLAG.EQ.0) GO TO 780
C     CHECK NUMBER OF RESIDUAL CALCULATIONS
      IF (ICOUNT.EQ.0) GO TO 760
      IF (COUNT.LT.ICOUNT) GO TO 760
      IFAIL1 = 15
      GO TO 820
C     ITERATE
  760 IF (ISTATE.LT.6) GO TO 120
      IF (IFLAG.EQ.2) GO TO 700
      GO TO 440
C     ITERATION CONVERGED
C     COMPUTE SOLUTION
  780 IF1 = IFAIL1
      ISTATE = 0
      CALL D02SAU(P,M,A,E,N,WP(1,4),NPOINT,WP(1,1),WP(1,2),WP(1,3)
     *            ,FCN,FCN1,FCN2,BC,BC1,YMAX,W(1,2*M+9)
     *            ,IW1,PRSOL,C,M1,IF1)
      IFAIL1 = IF1
      IF (IFAIL1.EQ.0) GO TO 820
      GO TO (800,800,320,340,360,300) IFAIL1
C     IMPOSS ERROR IN D02SAU
  800 IFAIL1 = 19
  820 IFAIL = IFAIL1
      RETURN
      END
