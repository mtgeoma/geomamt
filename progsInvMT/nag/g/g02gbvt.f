      DOUBLE PRECISION FUNCTION G02GBV(FV,Y,T)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15 REVISED. IER-937 (APR 1991).
C     MARK 15B REVISED. IER-955 (NOV 1991).
C
C     CALCULATES DEVIANCE FOR BINOMIAL GLM
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 FV, T, Y
C     .. Local Scalars ..
      DOUBLE PRECISION                 DEV
C     .. Intrinsic Functions ..
      INTRINSIC                        LOG
C     .. Executable Statements ..
      DEV = 0.0D0
      IF (Y.GT.0.0D0 .AND. Y.LT.T) THEN
         DEV = (Y*LOG(Y/FV)+(T-Y)*LOG((T-Y)/(T-FV)))
      ELSE IF (Y.GT.0.0D0) THEN
         DEV = Y*LOG(Y/FV)
      ELSE IF (T.GT.0.0D0) THEN
         DEV = (T-Y)*LOG((T-Y)/(T-FV))
      END IF
      IF (DEV.LT.0.0D0) DEV = 0.0D0
      G02GBV = 2.0D0*DEV
      RETURN
      END
