      DOUBLE PRECISION FUNCTION X02AKF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS  B**(EMIN-1)  (THE SMALLEST POSITIVE MODEL NUMBER)
C
      DOUBLE PRECISION X02CON
      DATA X02CON /2.22507385850721D-308 /
C     .. Executable Statements ..
      X02AKF = X02CON
      RETURN
      END
