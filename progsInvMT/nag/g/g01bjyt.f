      SUBROUTINE G01BJY(X,A,B,IABX)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Computes I(A,B,X) = Integral from t=0 to t=X of
C                          (t/X)**(A-1) * ((1-t)/(1-X))**(B-1)
C
C     Intended ranges of the input arguments:
C
C     1 .LT. A,B
C     0 .LT. X .LE. 0.5
C     X .LE. A/(A+B)
C
C     In this version which is required for the binomial
C     distribution, A and B are assumed to be integers.
C
C     .. Parameters ..
      DOUBLE PRECISION  EPS
      PARAMETER         (EPS=0.5D-6)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B, IABX, X
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, B1, FAC, R
      INTEGER           I, NSTEP
C     .. Executable Statements ..
C
C     Determine NSTEP, the number of steps in the backward recursion
C
      NSTEP = 0
      FAC = 1.0D0
      R = X/(1.0D0-X)
      A1 = A
      B1 = B
   20 CONTINUE
      FAC = FAC*R*B1/A1
      A1 = A1 + 1.0D0
      B1 = B1 - 1.0D0
      NSTEP = NSTEP + 1
      IF (FAC.GT.EPS) GO TO 20
C
C     Backward recursion with NSTEP iteration steps
C
      IABX = 0.0D0
      DO 40 I = 1, NSTEP
         A1 = A1 - 1.0D0
         B1 = B1 + 1.0D0
         IABX = (1.0D0+IABX)*R*B1/A1
   40 CONTINUE
C
C     Final normalization
C
      IABX = IABX*(1.0D0-X)/B
      RETURN
      END
