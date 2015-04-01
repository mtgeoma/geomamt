      DOUBLE PRECISION FUNCTION G01MAZ(X)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Parameters ..
      DOUBLE PRECISION                 FPI
      PARAMETER                        (FPI=
     *                             0.398942280401432677939946059934350D0
     *                                 )
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
C     .. Local Scalars ..
      DOUBLE PRECISION                 ARG, UFLOW
C     .. External Functions ..
      DOUBLE PRECISION                 X02AMF
      EXTERNAL                         X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC                        EXP, LOG
C     .. Executable Statements ..
C
      UFLOW = LOG(X02AMF())
      ARG = -0.5D0*(X**2)
      IF (ARG.LE.UFLOW) THEN
         G01MAZ = 0.0D0
      ELSE
         G01MAZ = FPI*EXP(ARG)
      END IF
      RETURN
      END
