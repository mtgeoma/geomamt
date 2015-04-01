      SUBROUTINE G01DHU(TIES,SCORES,N,X,R,IWRK)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     .. Scalar Arguments ..
      INTEGER           N
      CHARACTER         SCORES, TIES
C     .. Array Arguments ..
      DOUBLE PRECISION  R(N), X(N)
      INTEGER           IWRK(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP
      INTEGER           I, J, K
C     .. External Functions ..
      DOUBLE PRECISION  G01DHT
      EXTERNAL          G01DHT
C     .. Executable Statements ..
C
      IF (TIES.EQ.'L') THEN
         I = 1
         TEMP = G01DHT(SCORES,I,N)
         R(1) = TEMP
         DO 20 I = 2, N
            IF (X(IWRK(I-1)).EQ.X(IWRK(I))) THEN
               R(I) = TEMP
            ELSE
               TEMP = G01DHT(SCORES,I,N)
               R(I) = TEMP
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
   80       TEMP = G01DHT(SCORES,J,N)
            DO 100 K = I, J
               R(K) = TEMP
  100       CONTINUE
            I = J + 1
         ELSE
            TEMP = G01DHT(SCORES,I,N)
            R(I) = TEMP
            I = I + 1
         END IF
         IF (I.LT.N) THEN
            GO TO 40
         ELSE IF (I.EQ.N) THEN
            R(N) = G01DHT(SCORES,I,N)
         END IF
      END IF
      RETURN
      END
