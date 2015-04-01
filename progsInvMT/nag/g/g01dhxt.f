      SUBROUTINE G01DHX(N,X,IWRK)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     .. Scalar Arguments ..
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N)
      INTEGER           IWRK(N)
C     .. Local Scalars ..
      INTEGER           I, IF2, J
C     .. External Subroutines ..
      EXTERNAL          G05EHF
C     .. Executable Statements ..
      I = 1
   20 CONTINUE
      IF (X(IWRK(I)).EQ.X(IWRK(I+1))) THEN
         DO 40 J = I + 1, N - 1
            IF (X(IWRK(J)).LT.X(IWRK(J+1))) GO TO 60
   40    CONTINUE
         J = N
   60    IF2 = 0
         CALL G05EHF(IWRK(I),J-I+1,IF2)
         I = J + 1
      ELSE
         I = I + 1
      END IF
      IF (I.LT.N) GO TO 20
      RETURN
      END
