      SUBROUTINE C05AVF(X,FX,H,BOUNDL,BOUNDU,A,C,IND,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 10B REVISED. IER-398 (JAN 1983).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     FINDS INTERVAL [A,X] IN WHICH LIES ROOT OF F.
C     [BOUNDL,BOUNDU] CONTAINS [A,X],X IS INITIAL GUESS
C     AND H IS STARTING STEP.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C05AVF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, BOUNDL, BOUNDU, FX, H, X
      INTEGER           IFAIL, IND
C     .. Array Arguments ..
      DOUBLE PRECISION  C(11)
C     .. Local Scalars ..
      INTEGER           I
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SIGN, INT
C     .. Executable Statements ..
      IF (IND.NE.1 .AND. IND.NE.-1) GO TO 100
      IF (BOUNDL.LT.BOUNDU) GO TO 20
C     INPUT ERROR
      I = 1
      GO TO 480
   20 IF (X.GE.BOUNDL .AND. X.LE.BOUNDU) GO TO 40
C     INPUT ERROR
      I = 1
      GO TO 480
   40 IF (X+H.NE.X) GO TO 60
C     INPUT ERROR - H TOO SMALL
      I = 2
      GO TO 480
   60 IF ((X+H.GE.BOUNDL .AND. X+H.LE.BOUNDU)
     *    .OR. (X-H.GE.BOUNDL .AND. X-H.LE.BOUNDU)) GO TO 80
C     INPUT ERROR - H TOO LARGE
      I = 1
      GO TO 480
   80 C(3) = X
      C(5) = H
      C(6) = 0.D0
      C(7) = 0.D0
      A = X
      C(8) = X
      C(10) = X
      IF (IND.EQ.-1) GO TO 120
      IND = 2
      RETURN
  100 IF (IND.EQ.2 .OR. IND.EQ.3) GO TO 120
C     INPUT ERROR - IND INCORRECTLY SET
      I = 3
      GO TO 480
  120 IF (FX.NE.0.0D0) GO TO 140
C     HIT ROOT EXACTLY
      A = X
      GO TO 460
  140 IF (IND.EQ.3) GO TO 280
      C(4) = FX
      C(1) = FX
      C(9) = FX
      C(11) = FX
  160 IF (X.EQ.BOUNDL .AND. X+H.LT.BOUNDL) H = -H
      IF (X.EQ.BOUNDU .AND. X+H.GT.BOUNDU) H = -H
      IF (X.EQ.BOUNDL) C(7) = -1.0D0
      IF (X.EQ.BOUNDU) C(7) = 1.0D0
  180 X = X + H
C     TEST FOR ENDS OF INTERVAL
      IF (X.GT.BOUNDL) GO TO 220
      X = BOUNDL
      IF (C(7).EQ.-1.0D0 .OR. ABS(C(7)).GT.1.5D0) GO TO 360
      IF (C(7).EQ.0.0D0) GO TO 200
      C(7) = -2.0D0
      GO TO 260
  200 C(7) = -1.0D0
      GO TO 260
  220 IF (X.LT.BOUNDU) GO TO 260
      X = BOUNDU
      IF (C(7).EQ.1.0D0 .OR. ABS(C(7)).GT.1.5D0) GO TO 360
      IF (C(7).EQ.0.0D0) GO TO 240
      C(7) = 2.0D0
      GO TO 260
  240 C(7) = 1.0D0
  260 IND = 3
      RETURN
  280 C(2) = FX
      IF (SIGN(1.D0,C(1)).NE.SIGN(1.D0,C(2))) GO TO 440
      IF (ABS(X-BOUNDL).GT.ABS(C(8)-BOUNDL)) GO TO 300
      C(8) = X
      C(9) = FX
      GO TO 320
  300 IF (ABS(X-BOUNDU).GT.ABS(C(10)-BOUNDU)) GO TO 320
      C(10) = X
      C(11) = FX
  320 IF (ABS(C(1)).LE.ABS(C(2))) GO TO 340
C     CARRY ON IN SAME DIRECTION
      H = 2.0D0*H
      A = X
      C(1) = C(2)
      GO TO 180
C     BACK TRACK
  340 H = -2.0D0*H
      GO TO 180
  360 IF (C(6).LT.2.5D0) GO TO 420
C     GIVE UP
  380 I = 4
      IF (C(8).GT.C(10)) GO TO 400
      X = C(8)
      FX = C(9)
      A = C(10)
      C(1) = C(11)
      GO TO 480
  400 X = C(10)
      FX = C(11)
      A = C(8)
      C(1) = C(9)
      GO TO 480
C     TRY AGAIN WITH SHORTER STEP
  420 C(6) = C(6) + 1.0D0
      I = INT(C(6)+0.1D0)
      H = (-0.1D0)**I*C(5)
      IF (C(3).EQ.C(3)+H) GO TO 380
      X = C(3)
      A = X
      C(8) = A
      C(10) = A
      C(9) = C(4)
      C(11) = C(4)
      C(1) = C(4)
      C(7) = 0.D0
      GO TO 160
C     INTERVAL FOUND
  440 FX = C(2)
      A = C(8)
      C(1) = C(9)
      IF (ABS(X-C(10)).GE.ABS(X-C(8))) GO TO 460
      A = C(10)
      C(1) = C(11)
  460 I = 0
C     ERROR EXITS
  480 IND = 0
      IFAIL = P01ABF(IFAIL,I,SRNAME,0,P01REC)
      RETURN
      END
