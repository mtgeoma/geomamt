      SUBROUTINE G01BKW(X,A,LPAX)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Computes log(P(A,X)) where
C
C     P(A,X) = EXP(-X) * X**(A-1) / GAMMA(A)
C
C     Intended ranges of the input arguments X and A:
C
C     1.0 .LT. A .LE. 1E6
C     0.0 .LT. X .LE. 1E6
C
C     G01BJU returns log(1+X) for X.ge.0
C     G01BJV returns log(gamma(X+1)*exp(X)/X**X)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, LPAX, X
C     .. Local Scalars ..
      DOUBLE PRECISION  B, TERM
C     .. External Functions ..
      DOUBLE PRECISION  G01BJU, G01BJV
      EXTERNAL          G01BJU, G01BJV
C     .. Executable Statements ..
C
      B = A - 1.0D0
C
C     Compute TERM = log(X/(A-1))
C
      IF (B.GE.X) THEN
         TERM = G01BJU((B-X)/X)
      ELSE
         TERM = -G01BJU((X-B)/B)
      END IF
C
      LPAX = (B-X) - B*TERM - G01BJV(B)
      RETURN
      END
