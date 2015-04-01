      DOUBLE PRECISION FUNCTION D02KDS(X)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     EXP AVOIDING UNDERFLOW ERROR
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
C     .. External Functions ..
      DOUBLE PRECISION                 X02AMF
      EXTERNAL                         X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC                        EXP, LOG
C     .. Executable Statements ..
      D02KDS = 0.D0
      IF (X.GE.LOG(X02AMF())) D02KDS = EXP(X)
      RETURN
      END
