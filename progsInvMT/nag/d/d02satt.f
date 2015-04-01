      SUBROUTINE D02SAT(P,A,B,N,IFLAG)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     .. Scalar Arguments ..
      INTEGER           IFLAG, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N,2), B(N,2), P(N)
C     .. Local Scalars ..
      INTEGER           I, J
C     .. Executable Statements ..
      IF (IFLAG.EQ.1) GO TO 60
      J = 1
      DO 40 I = 1, N
         IF (B(I,1).EQ.0.0D0) GO TO 20
         P(J) = A(I,1)
         J = J + 1
   20    IF (B(I,2).EQ.0.0D0) GO TO 40
         P(J) = A(I,2)
         J = J + 1
   40 CONTINUE
      RETURN
   60 J = 1
      DO 100 I = 1, N
         IF (B(I,1).EQ.0.0D0) GO TO 80
         A(I,1) = P(J)
         J = J + 1
   80    IF (B(I,2).EQ.0.0D0) GO TO 100
         A(I,2) = P(J)
         J = J + 1
  100 CONTINUE
      RETURN
      END
