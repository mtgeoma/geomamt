      SUBROUTINE G01DHW(TIES,N,X,R,IWRK)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     .. Scalar Arguments ..
      INTEGER           N
      CHARACTER         TIES
C     .. Array Arguments ..
      DOUBLE PRECISION  R(N), X(N)
      INTEGER           IWRK(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP
      INTEGER           I, J, K
C     .. Executable Statements ..
C
      IF (TIES.EQ.'L') THEN
         I = 1
         TEMP = R(1)
         DO 20 I = 2, N
            IF (X(IWRK(I-1)).EQ.X(IWRK(I))) THEN
               R(I) = TEMP
            ELSE
               TEMP = R(I)
            END IF
   20    CONTINUE
      ELSE IF (TIES.EQ.'H') THEN
         I = 1
   40    CONTINUE
         IF (X(IWRK(I)).EQ.X(IWRK(I+1))) THEN
            DO 60 J = I + 1, N - 1
               IF (X(IWRK(J)).LT.X(IWRK(J+1))) GO TO 80
   60       CONTINUE
            J = N
   80       DO 100 K = I, J - 1
               R(K) = R(J)
  100       CONTINUE
            I = J + 1
         ELSE
            I = I + 1
         END IF
         IF (I.LT.N) GO TO 40
      END IF
      RETURN
      END
