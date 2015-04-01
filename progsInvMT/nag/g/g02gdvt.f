      DOUBLE PRECISION FUNCTION G02GDV(FV,Y,T)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 FV, T, Y
C     .. Intrinsic Functions ..
      INTRINSIC                        LOG
C     .. Executable Statements ..
      IF (FV.LE.0.0D0) THEN
         T = -1.0D0
         G02GDV = 0.0D0
      ELSE
         G02GDV = 2.0D0*(LOG(FV)+Y/FV)
      END IF
      RETURN
      END
