      SUBROUTINE G01BJX(X,A,B,JABX)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Computes J(A,B,X) = Integral from t=X to t=1 of
C                          (t/X)**(A-1) * ((1-t)/(1-X))**(B-1)
C
C     Intended ranges of the input arguments:
C
C     1 .LT. A,B
C     0 .LT. X .LE. 0.5
C     A/(A+B) .LT. X
C
C     In this version, which is required for the binomial
C     distribution, A and B are assumed to be integers.
C
C     .. Parameters ..
      DOUBLE PRECISION  EPS
      PARAMETER         (EPS=0.5D-6)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B, JABX, X
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, B1, FAC, R
      INTEGER           I, NSTEP
C     .. Executable Statements ..
C
C     Determine NSTEP, the number of steps in the forward recursion
C
      NSTEP = 0
      FAC = 1.0D0
      R = (1.0D0-X)/X
      A1 = A
      B1 = B
   20 CONTINUE
      FAC = FAC*R*A1/B1
      A1 = A1 - 1.0D0
      B1 = B1 + 1.0D0
      NSTEP = NSTEP + 1
      IF (FAC.GT.EPS) GO TO 20
C
C     Forward recursion with NSTEP iteration steps
C
      JABX = 0.0D0
      DO 40 I = 1, NSTEP
         A1 = A1 + 1.0D0
         B1 = B1 - 1.0D0
         JABX = (1.0D0+JABX)*R*A1/B1
   40 CONTINUE
C
C     Final normalization
C
      JABX = JABX*X/A
      RETURN
      END
