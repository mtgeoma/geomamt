      SUBROUTINE D02HBF(P,N1,PE,E,N,SOLN,M1,FCN,BC,RANGE,W,IW,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 9 REVISED. IER-306 (SEP 1981).
C     MARK 9C REVISED. IER-366 (MAY 1982).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. IER-619 (APR 1988).
C     BC, FCN, RANGE
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02HBF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IW, M1, N, N1
C     .. Array Arguments ..
      DOUBLE PRECISION  E(N), P(N1), PE(N1), SOLN(N,M1), W(N,IW)
C     .. Subroutine Arguments ..
      EXTERNAL          BC, FCN, RANGE
C     .. Scalars in Common ..
      INTEGER           ICASE, IFAIL2, IFAIL3, IW2
C     .. Arrays in Common ..
      DOUBLE PRECISION  W1(7)
      INTEGER           IW1(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  X, YMAX
      INTEGER           I, ICOUNT, IFAIL1, IWP, NPOINT
C     .. Local Arrays ..
      DOUBLE PRECISION  WP(2,6)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      LOGICAL           D02HBY
      EXTERNAL          P01ABF, D02HBY
C     .. External Subroutines ..
      EXTERNAL          D02HBR, D02HBS, D02HBT, D02HBU, D02HBV, D02HBW,
     *                  D02HBX, D02HBZ, D02SAV
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MOD
C     .. Common blocks ..
      COMMON            /AD02HB/IFAIL2, IFAIL3
      COMMON            /AD02SA/W1, IW1, IW2, ICASE
C     .. Executable Statements ..
      IFAIL2 = MOD(IFAIL/10,10)
      IFAIL3 = MOD(IFAIL/100,10)
C     TEST PARAMETERS
      IF (N1.LE.0 .OR. M1.LT.1 .OR. N.LT.N1 .OR. IW.LT.3*N+14+MAX(11,N))
     *    GO TO 420
      DO 40 I = 1, N1
         IF (PE(I).LE.0.0D0) GO TO 420
   40 CONTINUE
      DO 60 I = 1, N
         IF (E(I).LE.0.0D0) GO TO 420
   60 CONTINUE
C     SET UP ADDITIONAL PARAMETERS FOR D02SAV
C     AND COMMON VARIABLES
      DO 80 I = 1, 3
         WP(1,I) = 0.0D0
   80 CONTINUE
C     LET W(.,1) BE PF
      DO 100 I = 1, N
         W(I,1) = 1.0D0
  100 CONTINUE
C     LET W(.,2) BE DP
      DO 120 I = 1, N
         W(I,2) = 0.0D0
  120 CONTINUE
      X = 0.0D0
      YMAX = 0.0D0
      NPOINT = 2
      ICOUNT = 0
      IFAIL1 = 1
      ICASE = 2
      IWP = 2
      CALL D02SAV(P,N1,W,W,N,N1,PE,W(1,1),E,W(1,2)
     *            ,NPOINT,WP,IWP,X,X,ICOUNT,RANGE,D02HBU,BC,D02HBT,
     *            D02HBS,FCN,D02HBR,D02HBZ,D02HBY,YMAX,D02HBX,D02HBW,
     *            SOLN,M1,W(1,3),N,IW-2,IFAIL1)
      IF (IFAIL1.EQ.0) GO TO 460
      GO TO (360,360,360,140,160,360,360,360,180,
     *       200,360,220,240,260,360,280,300,320,
     *       340) IFAIL1
C     STEP TOO SHORT
  140 IFAIL1 = 2
      GO TO 380
C     NO INITIAL STEP POSSIBLE
  160 IFAIL1 = 3
      GO TO 380
C     STEP TOO SMALL IN JACOBIAN CALCULATION
  180 IFAIL1 = 4
      GO TO 380
C     NO INITIAL STEP IN JACOBIAN CALCULATION
  200 IFAIL1 = 5
      GO TO 380
C     INSIGNIFICANT COLUMN IN JACOBIAN
  220 IFAIL1 = 6
      GO TO 440
C     SVD FAIL
  240 IFAIL1 = 7
      GO TO 440
C     NEWTON FAILS TO CONVERGE
  260 IFAIL1 = 8
      GO TO 440
C     SERIOUS ERROR IN D02SAZ
  280 IFAIL1 = 9
      GO TO 440
C     SERIOUS ERROR IN D02SAW
  300 IFAIL1 = 10
      GO TO 440
C     SERIOUS ERROR IN D02SAX
  320 IFAIL1 = 11
      GO TO 440
C     SERIOUS ERROR IN D02SAU
  340 IFAIL1 = 12
      GO TO 440
C     SERIOUS ERROR IN D02SAV
  360 IFAIL1 = 13
      GO TO 440
  380 W(1,2) = WP(1,5)
      DO 400 I = 1, N
         W(I,1) = W(I,3)
  400 CONTINUE
      GO TO 440
  420 IFAIL1 = 1
  440 CALL D02HBV(IFAIL1)
  460 IFAIL = P01ABF(IFAIL,IFAIL1,SRNAME,0,P01REC)
      RETURN
      END
