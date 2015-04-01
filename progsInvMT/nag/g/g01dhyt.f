      SUBROUTINE G01DHY(N,X,R,IWRK)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     .. Scalar Arguments ..
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  R(N), X(N)
      INTEGER           IWRK(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  AV, SUM
      INTEGER           I, J, K, NTIE
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
      I = 1
   20 CONTINUE
      IF (X(IWRK(I)).EQ.X(IWRK(I+1))) THEN
         NTIE = 2
         DO 40 J = I + 1, N - 1
            IF (X(IWRK(J)).LT.X(IWRK(J+1))) GO TO 60
            NTIE = NTIE + 1
   40    CONTINUE
         J = N
   60    SUM = 0.0D0
         DO 80 K = I, I - 1 + NTIE
            SUM = SUM + R(K)
   80    CONTINUE
         AV = SUM/DBLE(NTIE)
         DO 100 K = I, I - 1 + NTIE
            R(K) = AV
  100    CONTINUE
         I = J + 1
      ELSE
         I = I + 1
      END IF
      IF (I.LT.N) GO TO 20
      RETURN
      END
