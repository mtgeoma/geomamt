      SUBROUTINE C05AGF(X,HH,EPS,ETA,F,A,B,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 9 REVISED. IER-300 (SEP 1981).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     DRIVER FOR C05AVF AND C05AZF
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C05AGF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B, EPS, ETA, HH, X
      INTEGER           IFAIL
C     .. Function Arguments ..
      DOUBLE PRECISION  F
      EXTERNAL          F
C     .. Local Scalars ..
      DOUBLE PRECISION  FX, H, LB, UB, Y
      INTEGER           I, IFAIL1, IND, IR
C     .. Local Arrays ..
      DOUBLE PRECISION  C(17)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          C05AVF, C05AZF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      H = HH
      IF (EPS.GT.0.0D0 .AND. X+H.NE.X) GO TO 20
C     INPUT ERROR
      I = 1
      GO TO 320
   20 FX = F(X)
      IF (ABS(FX).GE.ETA) GO TO 40
C     ROOT AT INITIAL POINT WITHIN ETA
      I = 0
      A = X
      B = X
      GO TO 320
C     LOOK FOR INTERVAL
   40 IND = -1
      IFAIL1 = 1
      LB = X - 256.0D0*ABS(H)
      UB = X + 256.0D0*ABS(H)
   60 CALL C05AVF(X,FX,H,LB,UB,A,C,IND,IFAIL1)
      IF (IND.EQ.0) GO TO 80
      FX = F(X)
      IF (ABS(FX).GE.ETA) GO TO 60
C     ROOT WITHIN ETA
      A = X
      B = X
      I = 0
      GO TO 320
   80 IF (IFAIL1.EQ.0) GO TO 140
      GO TO (100,100,100,120) IFAIL1
C     IMPOSSIBLE EXIT FROM C05AVF
  100 I = 5
      GO TO 320
C     C05AVF COULD NOT FIND AN INTERVAL
  120 I = 2
      A = C(8)
      B = C(10)
      GO TO 320
  140 IF (A.NE.X) GO TO 160
      IF (ETA.GT.0.0D0) GO TO 100
C     ROOT FOUND EXACTLY
      I = 0
      A = X
      B = X
      GO TO 320
C     INTERVAL FOUND - NOW FIND ROOT
  160 B = X
      Y = A
      IND = -1
      IR = 0
      IFAIL1 = 1
  180 CALL C05AZF(X,Y,FX,EPS,IR,C,IND,IFAIL1)
      IF (IND.EQ.0) GO TO 200
      IF (IND.NE.4) GO TO 220
      FX = F(X)
      IF (ABS(FX).GE.ETA) GO TO 180
C     ROOT WITHIN ETA
      A = X
      B = X
      I = 0
      GO TO 320
  200 IF (IFAIL1.EQ.0) GO TO 280
      GO TO (220,220,220,240,260) IFAIL1
C     IMPOSSIBLE EXIT FROM C05AZF
  220 I = 6
      GO TO 320
C     PROBABLY A POLE
  240 I = 3
      GO TO 320
C     TOO MUCH ACCURACY REQUESTED
  260 I = 4
      GO TO 320
  280 IF (X.NE.Y) GO TO 300
      IF (ETA.GT.0.0D0) GO TO 220
C     ROOT FOUND EXACTLY
      A = X
      B = X
      I = 0
      GO TO 320
C     SUCCESS
  300 I = 0
  320 IF (A.LT.B) GO TO 340
      Y = A
      A = B
      B = Y
  340 IFAIL = P01ABF(IFAIL,I,SRNAME,0,P01REC)
      RETURN
      END
