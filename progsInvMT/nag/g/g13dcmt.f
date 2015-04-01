      SUBROUTINE G13DCM(K,P,Q,PHI,THETA,QQ,GAMMA,KR,Z,KW,MAT,B,WORK,
     *                  IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Scalar Arguments ..
      INTEGER           IFAIL, K, KR, KW, P, Q
C     .. Array Arguments ..
      DOUBLE PRECISION  B(KW), GAMMA(K,KR), MAT(KW,KW), PHI(K,P*K+1),
     *                  QQ(K,K), THETA(K,Q*K+1), WORK(KW), Z(KW)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           I, I2, J, J2, K2, K3, K4, L, L2, L4, M
C     .. External Subroutines ..
      EXTERNAL          F04ARF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C
C     This auxiliary routine calculates the theoretical cross
C     covariances between the W(t)'s and between the W(t)'s and the
C     E(t)'s
C
      K3 = K*K
C
C     generate the first q - 1 gamma's
C
C     first set GAMMA(0) = QQ
C
      DO 40 J = 1, K
         DO 20 I = 1, K
            GAMMA(I,J) = QQ(I,J)
   20    CONTINUE
   40 CONTINUE
C
      DO 160 M = 1, Q - 1
         DO 140 I = 1, K
            DO 120 J = 1, K
               SUM = 0.0D0
               DO 60 I2 = 1, K
                  SUM = SUM - THETA(I,(M-1)*K+I2)*QQ(I2,J)
   60          CONTINUE
               DO 100 K2 = 1, MIN(P,M)
                  DO 80 J2 = 1, K
                     SUM = SUM + PHI(I,(K2-1)*K+J2)*GAMMA(J2,(M-K2)*K+J)
   80             CONTINUE
  100          CONTINUE
               GAMMA(I,M*K+J) = SUM
  120       CONTINUE
  140    CONTINUE
  160 CONTINUE
C
C     calculate C(0) (lower triangle), C(1), ... , C(p-1)
C
C     first initialise MAT to zero
C
      DO 200 J = 1, KW
         DO 180 I = 1, KW
            MAT(I,J) = 0.0D0
  180    CONTINUE
  200 CONTINUE
      K4 = K*(K+1)/2
C
      DO 340 M = 1, P - 1
         DO 320 I = 1, K
            DO 300 J = 1, K
               L = K4 + (M-1)*K3 + (I-1)*K + J
C
C              first calculate right-hand side vector Z
C
               SUM = 0.0D0
               DO 240 L4 = M, Q
                  DO 220 I2 = 1, K
                     SUM = SUM - GAMMA(I,(L4-M)*K+I2)*THETA(J,(L4-1)
     *                     *K+I2)
  220             CONTINUE
  240          CONTINUE
               Z(L) = SUM
C
C              now set up MAT array
C
               MAT(L,L) = 1.0D0
               DO 280 I2 = 1, P
                  DO 260 K2 = 1, K
                     IF (M.EQ.I2) THEN
                        IF (I.GE.K2) THEN
                           L4 = (I-1)*I/2 + K2
                        ELSE
                           L4 = K2*(K2-1)/2 + I
                        END IF
                     END IF
                     IF (M.GT.I2) L4 = K4 + (M-I2-1)*K3 + (I-1)*K + K2
                     IF (M.LT.I2) L4 = K4 + (I2-M-1)*K3 + (K2-1)*K + I
                     MAT(L,L4) = MAT(L,L4) - PHI(J,(I2-1)*K+K2)
  260             CONTINUE
  280          CONTINUE
C
  300       CONTINUE
  320    CONTINUE
  340 CONTINUE
C
C     now unravel the equation for C(0) (lower triangle)
C
      DO 680 I = 1, K
         DO 660 J = 1, I
            SUM = QQ(I,J)
            L = (I-1)*I/2 + J
C
C           first set up right-hand side vector Z
C
            DO 400 M = 1, Q
               DO 380 I2 = 1, K
                  DO 360 J2 = 1, K
                     SUM = SUM + THETA(I,(M-1)*K+I2)*QQ(I2,J2)*THETA(J,
     *                     (M-1)*K+J2)
  360             CONTINUE
  380          CONTINUE
  400       CONTINUE
C
            DO 480 M = 1, P
               DO 460 L2 = M, Q
                  DO 440 I2 = 1, K
                     DO 420 J2 = 1, K
                        SUM = SUM - PHI(I,(M-1)*K+I2)*GAMMA(I2,(L2-M)
     *                        *K+J2)*THETA(J,(L2-1)*K+J2)
  420                CONTINUE
  440             CONTINUE
  460          CONTINUE
  480       CONTINUE
C
            DO 560 M = 1, P
               DO 540 L2 = M, Q
                  DO 520 I2 = 1, K
                     DO 500 J2 = 1, K
                        SUM = SUM - PHI(J,(M-1)*K+I2)*GAMMA(I2,(L2-M)
     *                        *K+J2)*THETA(I,(L2-1)*K+J2)
  500                CONTINUE
  520             CONTINUE
  540          CONTINUE
  560       CONTINUE
C
            Z(L) = SUM
C
            MAT(L,L) = 1.0D0
            DO 640 M = 1, P
               DO 620 L2 = 1, P
                  DO 600 I2 = 1, K
                     DO 580 J2 = 1, K
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
                        MAT(L,L4) = MAT(L,L4) - PHI(I,(M-1)*K+I2)*PHI(J,
     *                              (L2-1)*K+J2)
  580                CONTINUE
  600             CONTINUE
  620          CONTINUE
  640       CONTINUE
C
  660    CONTINUE
  680 CONTINUE
C
      IF (P.GE.1) THEN
         IFAIL = 1
         CALL F04ARF(MAT,KW,Z,KW,B,WORK,IFAIL)
         IF (IFAIL.EQ.0) THEN
            DO 700 I = 1, KW
               Z(I) = B(I)
  700       CONTINUE
         END IF
      ELSE
         IFAIL = 0
      END IF
      RETURN
C
      END
