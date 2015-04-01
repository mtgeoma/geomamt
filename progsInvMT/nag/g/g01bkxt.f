      SUBROUTINE G01BKX(X,A,JAX)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Computes J(A,X) = Integral from t=X to t=infinity of
C                        (t/X)**(A-1) * EXP(X-t)
C
C     Intended ranges of the input arguments X and A:
C
C     0 .LT. A .LE. X .LE. 1E6
C
C     In this version, which is required for the Poisson
C     distribution, A is assumed to be an integer.
C
C     .. Parameters ..
      DOUBLE PRECISION  EPS
      PARAMETER         (EPS=0.5D-6)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, JAX, X
C     .. Local Scalars ..
      DOUBLE PRECISION  B, FAC
      INTEGER           I, NSTEP
C     .. Executable Statements ..
C
C     Determine NSTEP, the number of steps in the forward recursion:
C
      NSTEP = 0
      FAC = 1.0D0
      B = A
   20 CONTINUE
      B = B - 1.0D0
      FAC = FAC*B/X
      NSTEP = NSTEP + 1
      IF (FAC.GT.EPS) GO TO 20
C
C     Forward recursion with NSTEP iteration steps:
C
      JAX = 1.0D0
      DO 40 I = 2, NSTEP
         B = B + 1.0D0
         JAX = 1.0D0 + JAX*B/X
   40 CONTINUE
      RETURN
      END
