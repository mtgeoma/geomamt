      DOUBLE PRECISION FUNCTION D02KDU(X,COEFFN)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     COEFFN
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
C     .. Subroutine Arguments ..
      EXTERNAL                         COEFFN
C     .. Scalars in Common ..
      DOUBLE PRECISION                 BP, LAMDA, MINSC, ONE, PI, PSIGN,
     *                                 TWO, ZER
      INTEGER                          JINT
C     .. Arrays in Common ..
      DOUBLE PRECISION                 YL(3), YR(3)
C     .. Local Scalars ..
      DOUBLE PRECISION                 DQDL, P, Q
C     .. Intrinsic Functions ..
      INTRINSIC                        SQRT
C     .. Common blocks ..
      COMMON                           /AD02KD/ZER, ONE, TWO, PI, LAMDA,
     *                                 PSIGN, MINSC, BP, YL, YR, JINT
C     .. Executable Statements ..
      CALL COEFFN(P,Q,DQDL,X,LAMDA,JINT)
      Q = P*Q
      P = MINSC*P
      P = P*P
      D02KDU = SQRT(SQRT(P*P+Q*Q))
      RETURN
      END
