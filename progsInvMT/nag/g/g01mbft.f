      DOUBLE PRECISION FUNCTION G01MBF(X)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     G01MBF returns the reciprocal of Mills ratio
C     i.e. it returns the value of Z(X)/Q(X).
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
C     .. Local Scalars ..
      DOUBLE PRECISION                 A, A0, A1, A2, AA, B0, B1, B2, C,
     *                                 R, S, T, TOL, Z
C     .. External Functions ..
      DOUBLE PRECISION                 S15ADZ, X01AAF, X02AJF, X02AMF
      EXTERNAL                         S15ADZ, X01AAF, X02AJF, X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, EXP, LOG, SQRT
C     .. Executable Statements ..
C
      IF (X.EQ.0.0D0) THEN
         G01MBF = SQRT(2.0D0/X01AAF(0.0D0))
      ELSE IF (X.LT.0.0D0) THEN
         C = 0.5D0*X*X
         IF (-C.LE.LOG(X02AMF()*2.0D0)) THEN
            G01MBF = 0.0D0
         ELSE
            AA = X/(SQRT(2.0D0))
            Z = (2.0D0*EXP(C)) - S15ADZ(AA)
            G01MBF = SQRT(2.0D0/X01AAF(0.0D0))/Z
         END IF
      ELSE IF (X.LE.9.0D0) THEN
         AA = X/(SQRT(2.0D0))
         G01MBF = SQRT(2.0D0/X01AAF(0.0D0))/S15ADZ(AA)
      ELSE
C
C        SWAN - Applied Statistics, 1969
C
         TOL = X02AJF()
         A = 2.0D0
         R = X
         S = X
         B1 = X
         A1 = X*X + 1.0D0
         A2 = X*(A1+2.0D0)
         B2 = A1 + 1.0D0
         T = A2/B2
   20    CONTINUE
         A = A + 1.0D0
         A0 = A1
         A1 = A2
         A2 = X*A1 + A*A0
         B0 = B1
         B1 = B2
         B2 = X*B1 + A*B0
         R = S
         S = T
         T = A2/B2
         IF (ABS(T-R)/T.GT.TOL .AND. ABS(T-S)/T.GT.TOL) GO TO 20
         G01MBF = T
      END IF
      RETURN
      END
