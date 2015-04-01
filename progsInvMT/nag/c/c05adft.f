      SUBROUTINE C05ADF(A,B,EPS,ETA,F,X,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 10C REVISED. IER-422 (JUL 1983).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     DRIVER FOR C05AZF
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C05ADF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B, EPS, ETA, X
      INTEGER           IFAIL
C     .. Function Arguments ..
      DOUBLE PRECISION  F
      EXTERNAL          F
C     .. Local Scalars ..
      DOUBLE PRECISION  FX, Y
      INTEGER           IFAIL1, IND, IR
C     .. Local Arrays ..
      DOUBLE PRECISION  C(17)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          C05AZF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SIGN
C     .. Executable Statements ..
      IFAIL1 = 1
C     INPUT ERROR
      IF (A.EQ.B .OR. EPS.LE.0.D0) GO TO 120
      X = A
      FX = F(X)
      IFAIL1 = 0
C     ZERO AT INITIAL POINT
      IF (ABS(FX).LT.ETA .OR. FX.EQ.0.D0) GO TO 120
      Y = X
      C(1) = FX
      X = B
      FX = F(X)
C     ZERO AT INITIAL POINT
      IF (ABS(FX).LT.ETA .OR. FX.EQ.0.D0) GO TO 120
      IFAIL1 = 1
C     NO ROOT IN RANGE
      IF (SIGN(1.D0,FX).EQ.SIGN(1.D0,C(1))) GO TO 120
      IR = 1
      IND = -1
   20 CALL C05AZF(X,Y,FX,EPS,IR,C,IND,IFAIL1)
      IF (IND.EQ.0) GO TO 40
      FX = F(X)
      IF (ABS(FX).GE.ETA .AND. FX.NE.0.D0) GO TO 20
C     ZERO HIT EXACTLY
      IFAIL1 = 0
      GO TO 120
   40 IF (IFAIL1.EQ.0) GO TO 120
      GO TO (60,60,60,100,80) IFAIL1
C     IMPOSSIBLE EXIT
   60 IFAIL1 = 4
      GO TO 120
C     TOO MUCH ACCURACY REQUESTED
   80 IFAIL1 = 2
      GO TO 120
C     PROBABLY A POLE
  100 IFAIL1 = 3
  120 IFAIL = P01ABF(IFAIL,IFAIL1,SRNAME,0,P01REC)
      RETURN
      END
