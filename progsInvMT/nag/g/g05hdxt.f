      SUBROUTINE G05HDX(K,IP,IQ,PAR,NPAR,M,QQ,IK,MA,KW,Z,B,GAMMA,MAT,
     *                  WORK,IFAULT)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Scalar Arguments ..
      INTEGER           IFAULT, IK, IP, IQ, K, KW, M, NPAR
C     .. Array Arguments ..
      DOUBLE PRECISION  B(KW), GAMMA(K,(IQ+1)*K), MA(M,M), MAT(KW,KW),
     *                  PAR(NPAR), QQ(IK,K), WORK(KW), Z(KW)
C     .. Local Scalars ..
      INTEGER           I, I2, IJ, J, J2, K2, K4, KP
C     .. External Subroutines ..
      EXTERNAL          F06EFF, F06FBF, F06QFF, G05HDY
C     .. Executable Statements ..
C
      CALL G05HDY(K,IP,IQ,IK,PAR,NPAR,QQ,GAMMA,Z,KW,MAT,B,WORK,IFAULT)
      IF (IFAULT.GT.0) RETURN
C
C     Set up lower triangle of GAMMA(0) and store as MA
C
      DO 20 I = 1, M
         CALL F06FBF(M,0.0D0,MA(1,I),1)
   20 CONTINUE
C
C     Set up pure AR part of GAMMA(0)
C
      K4 = K*(K+1)/2
      KP = IP*K
      K2 = K*K
      DO 100 J = 1, IP
         DO 80 I = J, IP
            IJ = (I-J-1)*K2
            DO 60 I2 = 1, K
               DO 40 J2 = 1, K
                  IF (I.EQ.J) THEN
                     MA((I-1)*K+I2,(J-1)*K+J2) = Z(I2*(I2-1)/2+J2)
                  ELSE
                     MA((I-1)*K+I2,(J-1)*K+J2) = Z(K4+IJ+(I2-1)*K+J2)
                  END IF
   40          CONTINUE
   60       CONTINUE
   80    CONTINUE
  100 CONTINUE
C
C     set up pure MA part of GAMMA(0)
C
      DO 120 I = 1, IQ
         CALL F06QFF('G',K,K,QQ,IK,MA(KP+(I-1)*K+1,KP+(I-1)*K+1),M)
  120 CONTINUE
C
C     set up B matrix in GAMMA(0)
C
      DO 160 J = 1, IP
         DO 140 I = J, IQ
            CALL F06QFF('G',K,K,GAMMA(1,(I-J)*K+1),K,MA(KP+(I-1)
     *                  *K+1,(J-1)*K+1),M)
  140    CONTINUE
  160 CONTINUE
C
C     Set upper triangle of MA equal to lower triangle of MA
C
      DO 180 I = 1, M
         CALL F06EFF(M-I+1,MA(I,I),1,MA(I,I),M)
  180 CONTINUE
C
      RETURN
      END
