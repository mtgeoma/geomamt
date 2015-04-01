      SUBROUTINE G05EXZ(M,N,R,NR)
C     MARK 6 RELEASE  NAG COPYRIGHT 1976
C     MARK 7 REVISED IER-134 (DEC 1978)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C     THIS SETS UP THE INDEX PART OF THE REFERENCE VECTOR.
C     .. Scalar Arguments ..
      INTEGER           M, N, NR
C     .. Array Arguments ..
      DOUBLE PRECISION  R(NR)
C     .. Local Scalars ..
      DOUBLE PRECISION  HALF, ONE, X, Y, ZERO
      INTEGER           I, J, K
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, INT
C     .. Data statements ..
      DATA              HALF/0.5D0/, ZERO/0.0D0/, ONE/1.0D0/
C     .. Executable Statements ..
      R(1) = DBLE(N-3)
      R(2) = DBLE(M) + HALF
C     THIS IS IN CASE TRUNCATION IS TOWARDS -INFINITY.
      IF (INT(R(2)).GT.M) R(2) = R(2) - ONE
      X = R(NR)
      DO 20 I = N, NR
         R(I) = R(I)/X
   20 CONTINUE
      X = ONE/R(1)
      Y = ZERO
      J = N - 1
      K = N
      DO 80 I = 3, J
   40    IF (R(K).GT.Y) GO TO 60
         K = K + 1
         GO TO 40
   60    R(I) = DBLE(K) + HALF
         Y = Y + X
   80 CONTINUE
      RETURN
      END
