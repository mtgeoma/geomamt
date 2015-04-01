      DOUBLE PRECISION FUNCTION G01EMW(Y)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     COMPUTES UPPER LIMIT OF INNER INTEGRAL.
C
C     .. Parameters ..
      DOUBLE PRECISION                 L2PI
      PARAMETER                        (L2PI=
     *                              1.83787706640934548356065947281139D0
     *                                 )
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 Y
C     .. Scalars in Common ..
      DOUBLE PRECISION                 C, EPSPR1, Q, R, UFLOW, V
      INTEGER                          IR
C     .. Local Scalars ..
      DOUBLE PRECISION                 QY, X2, XU
C     .. Intrinsic Functions ..
      INTRINSIC                        LOG, MIN, SQRT
C     .. Common blocks ..
      COMMON                           /AG01EM/Q, R, V, C, EPSPR1,
     *                                 UFLOW, IR
C     .. Executable Statements ..
C
      QY = Q*Y
      XU = QY
      X2 = -2.0D0*LOG(EPSPR1/QY)
      IF (X2.GT.L2PI) XU = SQRT(X2-L2PI) + QY
      G01EMW = MIN(XU,7.0D0)
      RETURN
      END
