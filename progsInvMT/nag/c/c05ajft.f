      SUBROUTINE C05AJF(X,EPS,ETA,F,NFMAX,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     DRIVER FOR C05AXF
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C05AJF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS, ETA, X
      INTEGER           IFAIL, NFMAX
C     .. Function Arguments ..
      DOUBLE PRECISION  F
      EXTERNAL          F
C     .. Local Scalars ..
      DOUBLE PRECISION  FX, SCALE
      INTEGER           FCOUNT, I, IFAIL1, IND, IR
C     .. Local Arrays ..
      DOUBLE PRECISION  C(26)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          C05AXF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
C     .. Executable Statements ..
      IF (EPS.GT.0.0D0 .AND. NFMAX.GT.0) GO TO 20
C     INPUT ERROR
      I = 1
      GO TO 180
   20 SCALE = SQRT(SQRT(X02AJF()))*MAX(1.D0,ABS(X))
      FCOUNT = 0
      IR = 0
      IFAIL1 = 1
      IND = -1
   40 FX = F(X)
      FCOUNT = FCOUNT + 1
      IF (ABS(FX).GE.ETA .AND. FX.NE.0.0D0) GO TO 60
C     ROOT FOUND EXACTLY
      I = 0
      GO TO 180
   60 IF (FCOUNT.LT.NFMAX) GO TO 80
C     TOO MANY F CALLS
      I = 4
      GO TO 180
   80 CALL C05AXF(X,FX,EPS,IR,SCALE,C,IND,IFAIL1)
      IF (IND.NE.0) GO TO 40
      IF (IFAIL1.EQ.0) GO TO 160
      GO TO (100,100,120,140,140,140) IFAIL1
C     INPUT ERROR -IMPOSSIBLE
  100 I = 5
      GO TO 180
C     WRONG SCALE
  120 I = 2
      GO TO 180
C     TOO MUCH ACCURACY REQUESTED
  140 I = 3
      GO TO 180
  160 I = 0
  180 IFAIL = P01ABF(IFAIL,I,SRNAME,0,P01REC)
      RETURN
      END
