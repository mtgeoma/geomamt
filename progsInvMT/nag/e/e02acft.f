      SUBROUTINE E02ACF(X,Y,N,AA,M1,REF)
C     MARK 1 RELEASE.  NAG COPYRIGHT 1971
C     MARK 4.5 REVISED
C     MARK 5C REVISED
C     MARK 9B REVISED. IER-361 (JAN 1982)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14C REVISED. IER-877 (NOV 1990).
C     CALCULATES A MINIMAX POLYNOMIAL FIT TO A SET OF DATA POINTS
C     AS A
C     SERIES OF CHEBYSHEV POLYNOMIALS.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E02ACF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  REF
      INTEGER           M1, N
C     .. Array Arguments ..
      DOUBLE PRECISION  AA(M1), X(N), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ABSHI, AI, AI1, D, DENOM, H, HI, HIMAX, HMAX,
     *                  ONE, P1, P5, PREVH, RHI, RHI1, XI, XJ, XNEXTH,
     *                  ZERO
      INTEGER           I, I1, I2, IFAIL, IJ1, IMAX, IRI, IRJ, J, J1, K,
     *                  M, M2
C     .. Local Arrays ..
      DOUBLE PRECISION  A(100), RH(100), RX(100)
      INTEGER           IR(100)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, INT
C     .. Data statements ..
      DATA              ZERO/0.0D0/, ONE/1.0D0/, P5/0.5D0/, P1/0.1D0/
C     .. Executable Statements ..
C     ENFORCED HARD FAIL FOR OUT-OF-BOUNDS PARAMETERS
      IF (M1.GE.N .OR. M1.GE.100) IFAIL = P01ABF(0,1,SRNAME,0,P01REC)
      DO 20 I = 2, N
         IF (X(I).LE.X(I-1)) IFAIL = P01ABF(0,2,SRNAME,0,P01REC)
   20 CONTINUE
      M2 = M1 + 1
      M = M1 - 1
      PREVH = -ONE
      IR(1) = 1
      IR(M2) = N
      D = DBLE(N-1)/DBLE(M1)
      H = D
      IF (M.EQ.0) GO TO 60
      DO 40 I = 2, M1
         IR(I) = INT(H+P5) + 1
         H = H + D
   40 CONTINUE
   60 H = -ONE
      DO 80 I = 1, M2
         IRI = IR(I)
         RX(I) = X(IRI)
         A(I) = Y(IRI)
         RH(I) = -H
         H = -H
   80 CONTINUE
      DO 120 J = 1, M1
         I1 = M2
         AI1 = A(I1)
         RHI1 = RH(I1)
         I = M2
  100    I = I - 1
         IJ1 = I - J + 1
         DENOM = RX(I1) - RX(IJ1)
         AI = A(I)
         RHI = RH(I)
         A(I1) = (AI1-AI)/DENOM
         RH(I1) = (RHI1-RHI)/DENOM
         I1 = I
         AI1 = AI
         RHI1 = RHI
         IF (I-J) 120, 120, 100
  120 CONTINUE
      H = -A(M2)/RH(M2)
      DO 140 I = 1, M2
         A(I) = A(I) + RH(I)*H
  140 CONTINUE
      IF (M.EQ.0) GO TO 200
      J = M1
  160 J = J - 1
      XJ = RX(J)
      I = J
      AI = A(I)
      J = J + 1
      DO 180 I1 = J, M1
         AI1 = A(I1)
         A(I) = AI - XJ*AI1
         AI = AI1
         I = I1
  180 CONTINUE
      J = J - 1
      IF (J-1) 200, 200, 160
  200 CONTINUE
      HMAX = ABS(H)
      IF (HMAX.GT.PREVH) GO TO 220
      A(M2) = -HMAX
      GO TO 480
  220 A(M2) = HMAX
      PREVH = HMAX
      IMAX = IR(1)
      HIMAX = H
      J = 1
      IRJ = IR(J)
      DO 300 I = 1, N
         IF (I.EQ.IRJ) GO TO 280
         XI = X(I)
         HI = ZERO
         K = M2
  240    K = K - 1
         HI = HI*XI + A(K)
         IF (K-1) 260, 260, 240
  260    HI = HI - Y(I)
         ABSHI = ABS(HI)
         IF (ABSHI.LE.HMAX) GO TO 300
         HMAX = ABSHI
         HIMAX = HI
         IMAX = I
         GO TO 300
  280    IF (J.GE.M2) GO TO 300
         J = J + 1
         IRJ = IR(J)
  300 CONTINUE
      IF (IMAX.EQ.IR(1)) GO TO 480
      DO 320 I = 1, M2
         IF (IMAX.LT.IR(I)) GO TO 340
  320 CONTINUE
      I = M2
  340 I2 = INT(DBLE(I)*P5)
      I2 = I - 2*I2
      XNEXTH = H
      IF (I2.EQ.0) XNEXTH = -H
      IF (HIMAX*XNEXTH.LT.0.0D0) GO TO 360
      IR(I) = IMAX
      GO TO 60
  360 IF (IMAX.GE.IR(1)) GO TO 420
      J1 = M2
      J = M2
  380 J = J - 1
      IR(J1) = IR(J)
      J1 = J
      IF (J-1) 400, 400, 380
  400 IR(1) = IMAX
      GO TO 60
  420 IF (IMAX.LE.IR(M2)) GO TO 460
      J = 1
      DO 440 J1 = 2, M2
         IR(J) = IR(J1)
         J = J1
  440 CONTINUE
      IR(M2) = IMAX
      GO TO 60
  460 IR(I-1) = IMAX
      GO TO 60
  480 CONTINUE
      DO 500 I = 1, M1
         AA(I) = A(I)
  500 CONTINUE
      REF = A(M2)
      RETURN
      END
