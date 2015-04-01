      SUBROUTINE D02HAF(A,B,N,X,X1,TOL,FCN,SOLN,M1,W,IW,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     FCN
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02HAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL, X, X1
      INTEGER           IFAIL, IW, M1, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N,2), B(N,2), SOLN(N,M1), W(N,IW)
C     .. Subroutine Arguments ..
      EXTERNAL          FCN
C     .. Scalars in Common ..
      INTEGER           ICASE, IFAIL2, IFAIL3, IW2
C     .. Arrays in Common ..
      DOUBLE PRECISION  W1(7)
      INTEGER           IW1(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  YMAX
      INTEGER           I, ICOUNT, IFAIL1, IWP, J, NPOINT
C     .. Local Arrays ..
      DOUBLE PRECISION  WP(2,6)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      LOGICAL           D02HBY
      EXTERNAL          P01ABF, D02HBY
C     .. External Subroutines ..
      EXTERNAL          D02HAX, D02HAY, D02HAZ, D02HBR, D02HBT, D02HBU,
     *                  D02HBV, D02HBW, D02HBX, D02HBZ, D02SAV
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MOD
C     .. Common blocks ..
      COMMON            /AD02HB/IFAIL2, IFAIL3
      COMMON            /AD02SA/W1, IW1, IW2, ICASE
C     .. Executable Statements ..
      IFAIL2 = MOD(IFAIL/10,10)
      IFAIL3 = MOD(IFAIL/100,10)
C     TEST PARAMETERS
      IF (N.LE.0 .OR. TOL.LE.0.0D0 .OR. M1.LT.1 .OR. IW.LT.3*N+17+
     *    MAX(11,N)) GO TO 360
      J = 0
      ICOUNT = 0
      DO 20 I = 1, N
         IF (B(I,1).EQ.0.0D0) J = J + 1
         IF (B(I,2).EQ.0.0D0) ICOUNT = ICOUNT + 1
   20 CONTINUE
      IF ((J.EQ.N) .OR. (J.EQ.0) .OR. ((J+ICOUNT).NE.N)) GO TO 360
C     SET UP ADDITIONAL PARAMETERS FOR D02SAV
C     AND COMMON VARIABLES.LET W(.,1) BE P,
C     W(.,2) BE PE,W(.,3) BE E,W(.,4) BE PF AND
C     W(.,5) BE DP.
      DO 40 I = 1, 3
         WP(1,I) = 0.0D0
   40 CONTINUE
      DO 60 I = 1, N
         W(I,2) = TOL
         W(I,3) = TOL
         W(I,4) = 1.0D0
         W(I,5) = 0.0D0
   60 CONTINUE
      YMAX = 0.0D0
      NPOINT = 2
      ICOUNT = 0
      IFAIL1 = 1
      ICASE = 3
      IWP = 2
      CALL D02SAV(W(1,1),N,A,B,N,N,W(1,2),W(1,4),W(1,3),W(1,5)
     *            ,NPOINT,WP,IWP,X,X1,ICOUNT,D02HAY,D02HBU,D02HAZ,
     *            D02HBT,FCN,D02HAX,D02HBR,D02HBZ,D02HBY,YMAX,D02HBX,
     *            D02HBW,SOLN,M1,W(1,6),N,IW-5,IFAIL1)
      IF (IFAIL1.EQ.0) GO TO 400
      GO TO (300,300,300,80,100,300,300,300,120,
     *       140,300,160,180,200,300,220,240,260,
     *       280) IFAIL1
C     STEP TOO SHORT
   80 IFAIL1 = 2
      GO TO 320
C     NO INITIAL STEP POSSIBLE
  100 IFAIL1 = 3
      GO TO 320
C     STEP TOO SMALL IN JACOBIAN CALCULATION
  120 IFAIL1 = 4
      GO TO 320
C     NO INITIAL STEP IN JACOBIAN CALCULATION
  140 IFAIL1 = 5
      GO TO 320
C     INSIGNIFICANT COLUMN IN JACOBIAN
  160 IFAIL1 = 6
      GO TO 380
C     SVD FAIL
  180 IFAIL1 = 7
      GO TO 380
C     NEWTON FAILS TO CONVERGE
  200 IFAIL1 = 8
      GO TO 380
C     SERIOUS ERROR IN D02SAZ
  220 IFAIL1 = 9
      GO TO 380
C     SERIOUS ERROR IN D02SAW
  240 IFAIL1 = 10
      GO TO 380
C     SERIOUS ERROR IN D02SAX
  260 IFAIL1 = 11
      GO TO 380
C     SERIOUS ERROR IN D02SAU
  280 IFAIL1 = 12
      GO TO 380
C     SERIOUS ERROR IN D02SAV
  300 IFAIL1 = 13
      GO TO 380
  320 W(1,2) = WP(1,5)
      DO 340 I = 1, N
         W(I,1) = W(I,6)
  340 CONTINUE
      GO TO 380
  360 IFAIL1 = 1
  380 CALL D02HBV(IFAIL1)
  400 IFAIL = P01ABF(IFAIL,IFAIL1,SRNAME,0,P01REC)
      RETURN
      END
