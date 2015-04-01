      SUBROUTINE G13DSY(K,M,DEL,Y,IM)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER           IM, K, M
C     .. Array Arguments ..
      DOUBLE PRECISION  DEL(K,K), Y(IM,M*K*K)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           I, I2, J, J2, K2, L, MK2
C     .. Executable Statements ..
C
C     This subroutine sets Y = I ** (DEL ** DEL')
C     where ** denotes kronecker product
C
      K2 = K*K
      MK2 = M*K2
C
C     Initialise the Y array to zero
C
      DO 40 J = 1, MK2
         DO 20 I = 1, MK2
            Y(I,J) = 0.0D0
   20    CONTINUE
   40 CONTINUE
C
C     First set up DEL ** DEL' and put in the top left-hand
C     part of the array Y
C
      DO 120 J = 1, K
         DO 100 I = 1, K
            SUM = DEL(I,J)
            DO 80 J2 = 1, K
               DO 60 I2 = 1, K
                  Y((I-1)*K+I2,(J-1)*K+J2) = SUM*DEL(J2,I2)
   60          CONTINUE
   80       CONTINUE
  100    CONTINUE
  120 CONTINUE
C
C     Initialise the rest of Y
C
      DO 180 L = 2, M
         DO 160 J = 1, K2
            DO 140 I = 1, K2
               Y((L-1)*K2+I,(L-1)*K2+J) = Y(I,J)
  140       CONTINUE
  160    CONTINUE
  180 CONTINUE
      RETURN
C
      END
