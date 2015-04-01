      DOUBLE PRECISION FUNCTION D02PAY(X,H,XEND)
C     MARK 14 RELEASE.  NAG COPYRIGHT 1989.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 H, X, XEND
C     .. Local Scalars ..
      DOUBLE PRECISION                 TEMP
C     .. External Functions ..
      DOUBLE PRECISION                 D02PAZ, X02AJF
      EXTERNAL                         D02PAZ, X02AJF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, MAX, SIGN
C     .. Executable Statements ..
      TEMP = D02PAZ(X+H) - XEND
      IF (ABS(TEMP).LE.2.0D0*X02AJF()*MAX(ABS(X),ABS(XEND)))
     *    TEMP = X02AJF()*SIGN(1.0D0,XEND-X)
      D02PAY = TEMP
      RETURN
      END
