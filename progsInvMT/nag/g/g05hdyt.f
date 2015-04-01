      SUBROUTINE G05HDY(K,IP,IQ,IK,PAR,NPAR,QQ,GAMMA,Z,KW,MAT,B,WORK,
     *                  IERROR)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 17 REVISED. IER-1685 (JUN 1992).
C     .. Scalar Arguments ..
      INTEGER           IERROR, IK, IP, IQ, K, KW, NPAR
C     .. Array Arguments ..
      DOUBLE PRECISION  B(KW), GAMMA(K,(IQ+1)*K), MAT(KW,KW), PAR(NPAR),
     *                  QQ(IK,K), WORK(KW), Z(KW)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           I, I2, J, J2, K2, K3, K4, L, L2, L4, M, PK2
C     .. External Subroutines ..
      EXTERNAL          F04ARF, DCOPY, F06FBF, F06QFF, DGEMM
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C
C     This auxiliary routine calculates the theoretical cross
C     covariances between the W(t)'s and between the W(t)'s and the
C     E(t)'s
C
      IERROR = 0
      K3 = K*K
      PK2 = IP*K3
C
C     Generate the first q - 1 gamma's
C
C     First set GAMMA(0) = QQ
C
      CALL F06QFF('G',K,K,QQ,IK,GAMMA,K)
C
      DO 40 M = 1, IQ - 1
         CALL DGEMM('T','N',K,K,K,-1.0D0,PAR(PK2+(M-1)*K3+1),K,QQ,IK,
     *               0.0D0,GAMMA(1,M*K+1),K)
         DO 20 K2 = 1, MIN(IP,M)
            CALL DGEMM('T','N',K,K,K,1.0D0,PAR((K2-1)*K3+1),K,
     *                  GAMMA(1,(M-K2)*K+1),K,1.0D0,GAMMA(1,M*K+1),K)
   20    CONTINUE
   40 CONTINUE
C
C     Calculate C(0) (lower triangle), C(1), ... , C(p-1)
C
C     First initialise MAT to zero
C
      DO 60 I = 1, KW
         CALL F06FBF(KW,0.0D0,MAT(1,I),1)
   60 CONTINUE
      K4 = K*(K+1)/2
C
      DO 200 M = 1, IP - 1
         DO 180 I = 1, K
            DO 160 J = 1, K
               L = K4 + (M-1)*K3 + (I-1)*K + J
C
C              first calculate right-hand side vector Z
C
               SUM = 0.0D0
               DO 100 L4 = M, IQ
                  DO 80 I2 = 1, K
                     SUM = SUM - GAMMA(I,(L4-M)*K+I2)*PAR(PK2+(L4-1)
     *                     *K3+(J-1)*K+I2)
   80             CONTINUE
  100          CONTINUE
               Z(L) = SUM
C
C              now set up MAT array
C
               MAT(L,L) = 1.0D0
               DO 140 I2 = 1, IP
                  DO 120 K2 = 1, K
                     IF (M.EQ.I2) THEN
                        IF (I.GE.K2) THEN
                           L4 = (I-1)*I/2 + K2
                        ELSE
                           L4 = K2*(K2-1)/2 + I
                        END IF
                     END IF
                     IF (M.GT.I2) L4 = K4 + (M-I2-1)*K3 + (I-1)*K + K2
                     IF (M.LT.I2) L4 = K4 + (I2-M-1)*K3 + (K2-1)*K + I
                     MAT(L,L4) = MAT(L,L4) - PAR((I2-1)*K3+(J-1)*K+K2)
  120             CONTINUE
  140          CONTINUE
C
  160       CONTINUE
  180    CONTINUE
  200 CONTINUE
C
C     now unravel the equation for C(0) (lower triangle)
C
      DO 540 I = 1, K
         DO 520 J = 1, I
            SUM = QQ(I,J)
            L = (I-1)*I/2 + J
C
C           first set up right-hand side vector Z
C
            DO 260 M = 1, IQ
               DO 240 I2 = 1, K
                  DO 220 J2 = 1, K
                     SUM = SUM + PAR(PK2+(M-1)*K3+(I-1)*K+I2)*QQ(I2,J2)
     *                     *PAR(PK2+(M-1)*K3+(J-1)*K+J2)
  220             CONTINUE
  240          CONTINUE
  260       CONTINUE
C
            DO 340 M = 1, IP
               DO 320 L2 = M, IQ
                  DO 300 I2 = 1, K
                     DO 280 J2 = 1, K
                        SUM = SUM - PAR((M-1)*K3+(I-1)*K+I2)*GAMMA(I2,
     *                        (L2-M)*K+J2)*PAR(PK2+(L2-1)*K3+(J-1)*K+J2)
  280                CONTINUE
  300             CONTINUE
  320          CONTINUE
  340       CONTINUE
C
            DO 420 M = 1, IP
               DO 400 L2 = M, IQ
                  DO 380 I2 = 1, K
                     DO 360 J2 = 1, K
                        SUM = SUM - PAR((M-1)*K3+(J-1)*K+I2)*GAMMA(I2,
     *                        (L2-M)*K+J2)*PAR(PK2+(L2-1)*K3+(I-1)*K+J2)
  360                CONTINUE
  380             CONTINUE
  400          CONTINUE
  420       CONTINUE
C
            Z(L) = SUM
C
            MAT(L,L) = 1.0D0
            DO 500 M = 1, IP
               DO 480 L2 = 1, IP
                  DO 460 I2 = 1, K
                     DO 440 J2 = 1, K
                        IF (M.GT.L2) L4 = K4 + (M-L2-1)*K3 + (I2-1)*K +
     *                                    J2
                        IF (M.LT.L2) L4 = K4 + (L2-M-1)*K3 + (J2-1)*K +
     *                                    I2
                        IF (M.EQ.L2) THEN
                           IF (I2.GE.J2) THEN
                              L4 = (I2-1)*I2/2 + J2
                           ELSE
                              L4 = J2*(J2-1)/2 + I2
                           END IF
                        END IF
                        MAT(L,L4) = MAT(L,L4) - PAR((M-1)*K3+(I-1)*K+I2)
     *                              *PAR((L2-1)*K3+(J-1)*K+J2)
  440                CONTINUE
  460             CONTINUE
  480          CONTINUE
  500       CONTINUE
C
  520    CONTINUE
  540 CONTINUE
C
      IF (IP.GE.1) THEN
         IERROR = 1
         CALL F04ARF(MAT,KW,Z,KW,B,WORK,IERROR)
         IF (IERROR.EQ.0) CALL DCOPY(KW,B,1,Z,1)
      END IF
C
      RETURN
      END
