      DOUBLE PRECISION FUNCTION G01BJV(X)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Computes log(gamma(X+1)*exp(X)/X**X)
C
C     Intended range of the argument X:
C
C     0.LE.X
C
C     .. Parameters ..
      DOUBLE PRECISION                 A1, A2, A3, A4, A5, A6, A7
      PARAMETER                        (A1=1.0D0/12.0D0,
     *                                 A2=-1.0D0/360.0D0,
     *                                 A3=1.0D0/1260.0D0,
     *                                 A4=-1.0D0/1680.0D0,
     *                                 A5=1.0D0/1188.0D0,
     *                                 A6=-691.0D0/360360.0D0,
     *                                 A7=1.0D0/156.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
C     .. Local Scalars ..
      DOUBLE PRECISION                 RX, RX2
      INTEGER                          IFAIL
C     .. External Functions ..
      DOUBLE PRECISION                 S14ABF, X01AAF
      EXTERNAL                         S14ABF, X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC                        LOG
C     .. Executable Statements ..
      IF (X.EQ.0.0D0) THEN
         G01BJV = 0.0D0
      ELSE IF (X.LE.2.0D0) THEN
C
C        Use log gamma function (S14ABF)
C
         IFAIL = 0
         G01BJV = S14ABF(X+1.0D0,IFAIL) + X - X*LOG(X)
      ELSE
C
C        Use asymptotic expansion
C
         RX = 1.0D0/X
         RX2 = RX*RX
         G01BJV = 0.5D0*LOG(2.0D0*X01AAF(0.0D0)*X) + RX*((((((A7*RX2+A6)
     *            *RX2+A5)*RX2+A4)*RX2+A3)*RX2+A2)*RX2+A1)
      END IF
      RETURN
      END
