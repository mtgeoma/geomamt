      DOUBLE PRECISION FUNCTION G02GCV(FV,Y,T)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15 REVISED. IER-937 (APR 1991).
C
C     CALCULATES DEVIANCE FOR POISSON GLM
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 FV, T, Y
C     .. Local Scalars ..
      DOUBLE PRECISION                 DEV
C     .. Intrinsic Functions ..
      INTRINSIC                        LOG
C     .. Executable Statements ..
      DEV = 0.0D0
      IF (FV.LE.0.0D0) THEN
         T = -1.0D0
      ELSE IF (Y.GT.0.0D0) THEN
         DEV = Y*LOG(Y/FV) - (Y-FV)
      ELSE
         DEV = FV
      END IF
      IF (DEV.LT.0.0D0) DEV = 0.0D0
      G02GCV = 2.0D0*DEV
      RETURN
      END
