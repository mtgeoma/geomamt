      DOUBLE PRECISION FUNCTION G01JDU(U)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     Function for integrator for Imhof's procedure for N.le.5 -
C     called by D01AJF.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 U
C     .. Scalars in Common ..
      DOUBLE PRECISION                 C
      INTEGER                          M
C     .. Arrays in Common ..
      DOUBLE PRECISION                 A(5)
C     .. Local Scalars ..
      DOUBLE PRECISION                 PU, QU, SGN, SUM1, SUM2, TEMP
      INTEGER                          I
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, ATAN, EXP, LOG, SIGN, SIN
C     .. Common blocks ..
      COMMON                           /AG01JD/A, C, M
C     .. Executable Statements ..
      SUM1 = 0.0D0
      SUM2 = 0.0D0
      IF (U.EQ.0.0D0) THEN
         DO 20 I = 1, M
            SUM1 = SUM1 + A(I)
   20    CONTINUE
         G01JDU = 0.5D0*SUM1 - 0.5D0*C
      ELSE
         DO 40 I = 1, M
            SUM1 = SUM1 + ATAN(A(I)*U)
            SUM2 = SUM2 + LOG(1.0D0+U*U*(A(I)**2))
   40    CONTINUE
         QU = 0.5D0*SUM1 - 0.5D0*C*U
         PU = 0.25D0*SUM2
         TEMP = SIN(QU)
         SGN = SIGN(1.0D0,TEMP)
         G01JDU = (SGN/U)*EXP(LOG(ABS(TEMP))-PU)
      END IF
      RETURN
      END
