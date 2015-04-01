      DOUBLE PRECISION FUNCTION G01EMZ(X)
C
C     ASYMPTOTIC FUNCTION USED FOR PROBABILITY INTEGRAL.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
C     .. Scalars in Common ..
      DOUBLE PRECISION                 C, EPS, Q, R, UFLOW, V
      INTEGER                          IR
C     .. Local Scalars ..
      INTEGER                          IFAULT
C     .. External Functions ..
      DOUBLE PRECISION                 G01MAZ, S15ABF
      EXTERNAL                         G01MAZ, S15ABF
C     .. Common blocks ..
      COMMON                           /AG01EM/Q, R, V, C, EPS, UFLOW,
     *                                 IR
C     .. Executable Statements ..
C
      IFAULT = 0
      G01EMZ = R*G01MAZ(X)*((S15ABF(X,IFAULT)-S15ABF(X-Q,IFAULT))**IR)
C
      RETURN
      END
