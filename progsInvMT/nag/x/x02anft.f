      DOUBLE PRECISION FUNCTION X02ANF()
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     Returns the 'safe range' parameter for complex numbers,
C     i.e. the smallest positive model number Z such that
C     for any X which satisfies X.ge.Z and X.le.1/Z
C     the following can be computed without overflow, underflow or other
C     error
C
C        -W
C        1.0/W
C        SQRT(W)
C        LOG(W)
C        EXP(LOG(W))
C        Y**(LOG(W)/LOG(Y)) for any Y
C        ABS(W)
C
C     where W is any of cmplx(X,0), cmplx(0,X), cmplx(X,X),
C                   cmplx(1/X,0), cmplx(0,1/X), cmplx(1/X,1/X).
C
      DOUBLE PRECISION X02CON
      DATA X02CON /2.22507385850721D-308 /
C     .. Executable Statements ..
      X02ANF = X02CON
      RETURN
      END
