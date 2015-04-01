      SUBROUTINE G01BJW(X,A,B,LPABX)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Computes log(P(A,B,X)), where
C
C     P(A,B,X) = Gamma(A+B)/(Gamma(A)*Gamma(B)) * X**(A-1)*(1-X)**(B-1)
C
C     Intended range of the input arguments X, A, and B:
C
C     1 .LT. A,B
C     0.0 .LT. X .LT. 1.0
C     A*B/(A+B) .LE. 1E6
C
C     G01BJU returns log(1+X) for X.ge.0
C     G01BJV returns log(gamma(X+1)*exp(X)/X**X) for X.ge.0
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B, LPABX, X
C     .. Local Scalars ..
      DOUBLE PRECISION  D, DD, TERM1, TERM2
      INTEGER           IFAIL
C     .. Local Arrays ..
      DOUBLE PRECISION  AV(2), BV(2)
C     .. External Functions ..
      DOUBLE PRECISION  G01BJU, G01BJV
      EXTERNAL          G01BJU, G01BJV
C     .. External Subroutines ..
      EXTERNAL          X03AAF
C     .. Intrinsic Functions ..
      INTRINSIC         LOG
C     .. Executable Statements ..
C
C     Compute (A+B)*X - A, using additional precision in implementations
C     of low precision. In the call of X03AAF, the argument SW may be
C     changed to .FALSE. in double precision implementations.
C
      AV(1) = A
      AV(2) = B
      BV(1) = X
      BV(2) = X
      IFAIL = 0
      CALL X03AAF(AV,2,BV,2,2,1,1,-A,0.0D0,D,DD,.FALSE.,IFAIL)
C
C     Compute TERM1 = log((A+B)*X/A)
C             TERM2 = log((A+B)*(1-X)/B)
C
      IF (D.GE.0.0D0) THEN
         TERM1 = G01BJU(D/A)
         TERM2 = -G01BJU(D/((A+B)*(1.0D0-X)))
      ELSE
         TERM1 = -G01BJU(-D/((A+B)*X))
         TERM2 = G01BJU(-D/B)
      END IF
C
      LPABX = LOG(A+B) + G01BJV(A+B) - G01BJV(A) - G01BJV(B)
      LPABX = LPABX + (A-1.0D0)*TERM1 + (B-1.0D0)*TERM2
      RETURN
      END
