      SUBROUTINE G13DCN(K,Q,THETA,QQ,KR,TEMPK,TEMP)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Scalar Arguments ..
      INTEGER           K, KR, Q
C     .. Array Arguments ..
      DOUBLE PRECISION  QQ(K,K), TEMP(K,KR), TEMPK(K,KR), THETA(K,K*Q+1)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           I, I2, J, L, M
C     .. Executable Statements ..
C
C     This subroutine computes P(1/0)h for a pure moving average model
C
C     First calculate QQ * THETA(l)' for l = 1,2,...,q and store in TEMP
C
      DO 80 L = 1, Q
         DO 60 I = 1, K
            DO 40 J = 1, K
               SUM = 0.0D0
               DO 20 I2 = 1, K
                  SUM = SUM + QQ(I,I2)*THETA(J,(L-1)*K+I2)
   20          CONTINUE
               TEMP(I,(L-1)*K+J) = SUM
   40       CONTINUE
   60    CONTINUE
   80 CONTINUE
C
C     Now set up P(1/0)h
C
      DO 180 L = 1, Q
         DO 160 I = 1, K
            DO 140 J = 1, K
               SUM = 0.0D0
               DO 120 M = 0, Q - L
                  DO 100 I2 = 1, K
                     SUM = SUM + THETA(I,(M+L-1)*K+I2)*TEMP(I2,M*K+J)
  100             CONTINUE
  120          CONTINUE
               TEMPK(I,(L-1)*K+J) = SUM
  140       CONTINUE
  160    CONTINUE
  180 CONTINUE
      RETURN
C
      END
