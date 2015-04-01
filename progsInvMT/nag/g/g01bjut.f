      DOUBLE PRECISION FUNCTION G01BJU(X)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Computes log(1+X) with high relative accuracy
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
C     .. Local Scalars ..
      DOUBLE PRECISION                 Y
      INTEGER                          IFAIL
C     .. External Functions ..
      DOUBLE PRECISION                 S21BAF
      EXTERNAL                         S21BAF
C     .. Intrinsic Functions ..
      INTRINSIC                        LOG
C     .. Executable Statements ..
      IF (X.GE.1.0D0) THEN
         G01BJU = LOG(1.0D0+X)
      ELSE
         Y = 1.0D0 + 0.5D0*X
         IFAIL = 0
         G01BJU = X*S21BAF(Y*Y,1.0D0+X,IFAIL)
      END IF
      RETURN
      END
