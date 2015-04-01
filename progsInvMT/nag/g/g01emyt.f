      DOUBLE PRECISION FUNCTION G01EMY(X,Y)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     FUNCTION USED FOR PROBABILITY INTEGRAL. NOTE :- X AND Y SWOPPED.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X, Y
C     .. Scalars in Common ..
      DOUBLE PRECISION                 C, DUMMY3, Q, R, UFLOW, V
      INTEGER                          IR
C     .. Local Scalars ..
      DOUBLE PRECISION                 TEMP
      INTEGER                          IFAULT
C     .. External Functions ..
      DOUBLE PRECISION                 G01MAZ, S15ABF
      EXTERNAL                         G01MAZ, S15ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        EXP, LOG
C     .. Common blocks ..
      COMMON                           /AG01EM/Q, R, V, C, DUMMY3,
     *                                 UFLOW, IR
C     .. Executable Statements ..
C
      IF (Y.EQ.0.0D0) THEN
         G01EMY = 0.0D0
      ELSE
         TEMP = (V-1.0D0)*LOG(Y) - 0.5D0*V*(Y**2) + C
         IF (TEMP.LE.UFLOW) THEN
            G01EMY = 0.0D0
         ELSE
            IFAULT = 0
            G01EMY = EXP(TEMP)*R*(G01MAZ(X)*((S15ABF(X,IFAULT)
     *               -S15ABF(X-Q*Y,IFAULT))**IR))
         END IF
      END IF
      RETURN
      END
