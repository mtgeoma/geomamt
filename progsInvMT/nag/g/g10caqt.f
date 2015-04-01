      SUBROUTINE G10CAQ(Y,N,CHANGE)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     find 2-flats in Y() and apply splitting algorithm
C
C     .. Scalar Arguments ..
      INTEGER           N
      LOGICAL           CHANGE
C     .. Array Arguments ..
      DOUBLE PRECISION  Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  Y1
      INTEGER           I, I1
C     .. Local Arrays ..
      DOUBLE PRECISION  W(6)
C     .. External Subroutines ..
      EXTERNAL          G10CAP
C     .. Executable Statements ..
C
C     W() is a window 6 points wide which is slid along Y()
C
      DO 20 I = 1, 4
         W(I+2) = Y(I)
   20 CONTINUE
C
C     if Y(1)=Y(2) .NE. Y(3), treat first 2 like a 2-flat
C     with end pt rule
C
      W(2) = Y(3)
      I1 = 1
   40 IF (W(3).NE.W(4)) GO TO 80
      IF ((W(3)-W(2))*(W(5)-W(4)).GE.0.0D0) GO TO 80
C
C     W(3) and W(4) form a 2-flat
C
      IF (I1.LT.3) GO TO 60
C
C     apply right end pt rule at I1
C
      Y1 = 3.0D0*W(2) - 2.0D0*W(1)
      CALL G10CAP(Y1,W(3),W(2),Y(I1),CHANGE)
   60 IF (I1.GE.N-2) GO TO 80
C
C     apply left end pt rule at I1+1
C
      Y1 = 3.0D0*W(5) - 2.0D0*W(6)
      CALL G10CAP(Y1,W(4),W(5),Y(I1+1),CHANGE)
C
C     slide window
C
   80 DO 100 I = 1, 5
         W(I) = W(I+1)
  100 CONTINUE
      I1 = I1 + 1
      IF (I1.GE.N-2) GO TO 120
      W(6) = Y(I1+3)
      GO TO 40
C
C     apply rule to last 2 points if needed
C
  120 W(6) = W(3)
      IF (I1.LT.N) GO TO 40
      RETURN
      END
