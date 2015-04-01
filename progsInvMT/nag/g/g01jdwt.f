      DOUBLE PRECISION FUNCTION G01JDW(U,M,A)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     Function for integrator (G01JDX) for Imhof's procedure for
C     N.ge.6 - called by G01JDY.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 U
      INTEGER                          M
C     .. Array Arguments ..
      DOUBLE PRECISION                 A(M+1)
C     .. Local Scalars ..
      DOUBLE PRECISION                 PU, QU, SGN, SUM1, SUM2, TEMP
      INTEGER                          I
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, ATAN, EXP, LOG, SIGN, SIN
C     .. Executable Statements ..
C
      SUM1 = 0.0D0
      SUM2 = 0.0D0
C
      IF (U.EQ.0.0D0) THEN
         DO 20 I = 1, M
            SUM1 = SUM1 + A(I)
   20    CONTINUE
         G01JDW = 0.5D0*SUM1 - 0.5D0*A(M+1)
      ELSE
         DO 40 I = 1, M
            SUM1 = SUM1 + ATAN(A(I)*U)
            SUM2 = SUM2 + LOG(1.0D0+(A(I)**2)*(U**2))
   40    CONTINUE
         QU = 0.5D0*SUM1 - 0.5D0*A(M+1)*U
         PU = 0.25D0*SUM2
         TEMP = SIN(QU)
         SGN = SIGN(1.0D0,TEMP)
         G01JDW = (SGN/U)*EXP(LOG(ABS(TEMP))-PU)
      END IF
      RETURN
      END
