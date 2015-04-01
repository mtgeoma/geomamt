      SUBROUTINE F02BKZ(N,J,EPS,U)
C     MARK 7 RELEASE. NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     W.PHILLIPS. OXFORD UNIVERSITY COMPUTING SERVICE. 1 JUN 1977
C
C     GUESSVEC
C
C     AUXILIARY ROUTINE CALLED BY F02BKF AND F02BLF.
C     SETS UP VECTOR FOR INVERSE ITERATION
C
C     THIS ROUTINE REPLACES F02ATZ
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS
      INTEGER           J, N
C     .. Array Arguments ..
      DOUBLE PRECISION  U(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  X, Y
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Executable Statements ..
      X = SQRT(DBLE(N))
      Y = EPS/(X+1.0D0)
      U(1) = EPS
      IF (N.LT.2) GO TO 40
      DO 20 I = 2, N
         U(I) = Y
   20 CONTINUE
   40 U(J) = U(J) - EPS*X
      RETURN
      END
