      SUBROUTINE G01BKY(X,A,IAX)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Computes I(A,X) = Integral from t=0 to t=X of
C                        (t/X)**(A-1) * EXP(X-t)
C
C     Intended ranges of the input arguments X and A:
C
C     0 .LT. X .LE. A .LE. 1E6
C
C     .. Parameters ..
      DOUBLE PRECISION  EPS
      PARAMETER         (EPS=0.5D-6)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, IAX, X
C     .. Local Scalars ..
      DOUBLE PRECISION  B, FAC
      INTEGER           I, NSTEP
C     .. Executable Statements ..
C
C     Determine NSTEP, the number of steps in the backward recursion:
C
      NSTEP = 0
      FAC = 1.0D0
      B = A
   20 CONTINUE
      FAC = FAC*X/B
      B = B + 1.0D0
      NSTEP = NSTEP + 1
      IF (FAC.GT.EPS) GO TO 20
C
C     Backward recursion with NSTEP iteration steps:
C
      IAX = 0.0D0
      DO 40 I = 1, NSTEP
         B = B - 1.0D0
         IAX = (1.0D0+IAX)*X/B
   40 CONTINUE
      RETURN
      END
