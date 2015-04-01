      DOUBLE PRECISION FUNCTION S14BAZ(A,X,GETP,EPS,UNDFL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      DOUBLE PRECISION                 TWO, ONE, HALF, ZERO
      PARAMETER                        (TWO=2.0D0,ONE=1.0D0,HALF=0.5D0,
     *                                 ZERO=0.0D0)
      DOUBLE PRECISION                 RT2PI
      PARAMETER                        (RT2PI=2.5066282746310005024D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A, EPS, UNDFL, X
      LOGICAL                          GETP
C     .. Local Scalars ..
      DOUBLE PRECISION                 DIF, ETA, U, V, Y
      INTEGER                          IFAIL, S
C     .. External Functions ..
      DOUBLE PRECISION                 S14BAX, S14BAY, S15ADF
      EXTERNAL                         S14BAX, S14BAY, S15ADF
C     .. Intrinsic Functions ..
      INTRINSIC                        EXP, SQRT
C     .. Executable Statements ..
C
      IF (GETP) THEN
         S = -ONE
      ELSE
         S = ONE
      END IF
      DIF = (X-A)/A
      Y = A*S14BAX(DIF)
      IF (Y.LT.ZERO) Y = ZERO
      ETA = SQRT(TWO*Y/A)
      V = SQRT(Y)
      IF (X.LT.A) THEN
         ETA = -ETA
         V = -V
      END IF
      IFAIL = 0
      U = HALF*S15ADF(S*V,IFAIL)
      IF (-Y.GE.UNDFL) THEN
         V = S*EXP(-Y)*S14BAY(ETA,A,EPS)/(RT2PI*SQRT(A))
      ELSE
C        exp(-Y) underflows.
         V = ZERO
      END IF
      S14BAZ = U + V
      RETURN
      END
