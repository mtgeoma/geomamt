      INTEGER FUNCTION F11GBT(N,EPS,D,E2,X)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-----------------------------------------------------------------------
C
C     F11GBT - Utility for F11GBWP (Symmetric iterative solver suite)
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION        ZERO
      PARAMETER               (ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION        EPS, X
      INTEGER                 N
C     .. Array Arguments ..
      DOUBLE PRECISION        D(N), E2(N)
C     .. Local Scalars ..
      DOUBLE PRECISION        Q
      INTEGER                 I, K
C     .. Executable Statements ..
      Q = D(1) - X
      IF (Q.LT.ZERO) THEN
         K = 1
      ELSE
         K = 0
      END IF
      DO 20 I = 2, N
         IF (Q.EQ.ZERO) Q = EPS
         Q = D(I) - (X+E2(I-1)/Q)
         IF (Q.LT.ZERO) K = K + 1
   20 CONTINUE
      F11GBT = K
C
C     End of function F11GBT
C
      RETURN
      END
