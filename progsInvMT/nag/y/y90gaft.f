      DOUBLE PRECISION FUNCTION Y90GAF(K)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Parameters ..
      DOUBLE PRECISION                 ONE, TWO, THREE
      PARAMETER                        (ONE=1.0D0,TWO=2.0D0,THREE=3.0D0)
C     .. Scalar Arguments ..
      INTEGER                          K
C     .. Intrinsic Functions ..
      INTRINSIC                        MOD
C     .. Executable Statements ..
      IF (K.EQ.1) THEN
         Y90GAF = TWO
      ELSE IF (MOD(K,2).EQ.1) THEN
         Y90GAF = THREE
      ELSE
         Y90GAF = ONE
      END IF
      RETURN
      END
