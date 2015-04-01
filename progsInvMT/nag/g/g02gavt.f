      DOUBLE PRECISION FUNCTION G02GAV(FV,Y,T)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     CALCULATES DEVIANCE FOR NORMAL GLM
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 FV, T, Y
C     .. Executable Statements ..
      G02GAV = (Y-FV)**2
      END
