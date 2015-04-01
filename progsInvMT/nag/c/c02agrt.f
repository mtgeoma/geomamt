      DOUBLE PRECISION FUNCTION C02AGR(X,EXP)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     BASED ON THE ROUTINE  SETEXP, WRITTEN BY BRIAN T. SMITH
C
C     THIS FUNCTION COMPUTES THE VALUE  FRACTION(X) * B ** EXP
C     WHERE  B  IS THE BASE OF ENTITIES OF TYPE DOUBLE PRECISION.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
      INTEGER                          EXP
C     .. External Functions ..
      DOUBLE PRECISION                 C02AGY
      INTEGER                          C02AGX
      EXTERNAL                         C02AGY, C02AGX
C     .. Executable Statements ..
      C02AGR = C02AGY(X,EXP-C02AGX(X))
      RETURN
      END
