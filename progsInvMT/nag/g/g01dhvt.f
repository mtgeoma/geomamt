      SUBROUTINE G01DHV(N,X,R,IWRK)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     .. Scalar Arguments ..
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  R(N), X(N)
      INTEGER           IWRK(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  AV
      INTEGER           I, J, K, NTIE
C     .. Intrinsic Functions ..
      INTRINSIC         MOD, DBLE
C     .. Executable Statements ..
C
      I = 1
   20 CONTINUE
      IF (X(IWRK(I)).EQ.X(IWRK(I+1))) THEN
         NTIE = 2
         DO 40 J = I + 1, N - 1
            IF (X(IWRK(J)).LT.X(IWRK(J+1))) GO TO 60
            NTIE = NTIE + 1
   40    CONTINUE
         J = N
   60    IF (MOD(NTIE,2).EQ.0) THEN
            AV = DBLE(I+NTIE/2) - 0.5D0
         ELSE
            AV = DBLE(I+(NTIE-1)/2)
         END IF
         DO 80 K = I, J
            R(K) = AV
   80    CONTINUE
         I = J + 1
      ELSE
         R(I) = DBLE(I)
         I = I + 1
      END IF
      IF (I.LT.N) THEN
         GO TO 20
      ELSE IF (I.EQ.N) THEN
         R(N) = DBLE(N)
      END IF
      RETURN
      END
