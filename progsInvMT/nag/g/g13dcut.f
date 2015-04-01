      SUBROUTINE G13DCU(I6,P,X,N2,K,KR,KMAT,PP,QQ,IK,GAMMA,MAT,Z,KZ,B,
     *                  PHI,PHISTA,SIGMA,SIGSTA,L2,LSTAR,AA,BB,XX,
     *                  IFAILX)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     THIS ROUTINE COMPUTES THE INVERSE TRANSFORMATION FROM
C     ANSLEY AND KOHN'S PAPER
C
C     C.F. ANSLEY AND P. NEWBOLD (1980). 'MULTIVARIATE
C                                         PARTIAL AUTOCORRELATIONS'
C     .. Scalar Arguments ..
      INTEGER           I6, IFAILX, IK, K, KMAT, KR, KZ, N2, P
C     .. Array Arguments ..
      DOUBLE PRECISION  AA(K,K), B(K*K*(P+1)), BB(K,K), GAMMA(K,KR+K),
     *                  L2(K,K), LSTAR(K,K), MAT(KMAT,KMAT), PHI(KR,KR),
     *                  PHISTA(KR,KR), PP(K,P*K+1), QQ(IK,K),
     *                  SIGMA(K,K), SIGSTA(K,K), X(N2), XX(N2), Z(KZ)
C     .. Local Scalars ..
      DOUBLE PRECISION  D1, EPS, SM, SUM, SUM2
      INTEGER           I, I2, ID, IFAIL, IJAA, IK5, J, K2, KW, L, L3,
     *                  L4, M, M1, S
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          F03AEF, F04AFF, F04AGZ, F04ARF
C     .. Executable Statements ..
      IFAILX = 0
C     FIRST PUT X INTO PP
C
      DO 60 L = 1, P
         DO 40 I = 1, K
            DO 20 J = 1, K
               PP(I,(L-1)*K+J) = X(I6+(L-1)*K*K+(I-1)*K+J)
   20       CONTINUE
   40    CONTINUE
   60 CONTINUE
C
C     CALCULATE THE FIRST P C(J)'S.
C
C     INITIALISE MAT AND Z TO ZERO
C
      KW = K*K*(P+1)
      DO 100 I = 1, KW
         DO 80 J = 1, KW
            MAT(I,J) = 0.0D0
   80    CONTINUE
         Z(I) = 0.0D0
  100 CONTINUE
C
      M = -1
      DO 200 M1 = 1, P + 1
         M = M + 1
         DO 180 I = 1, K
            DO 160 J = 1, K
               L = M*K*K + (I-1)*K + J
C
C              FIRST FIND RIGHT HAND SIDE VECTOR Z
C
               IF (M.EQ.0) Z(L) = QQ(I,J)
C
C              SET UP MAT ARRAY
C
               IF (P.EQ.0) GO TO 160
               MAT(L,L) = 1.0D0
               DO 140 I2 = 1, P
                  DO 120 K2 = 1, K
                     IF (M.GE.I2) L4 = (M-I2)*K*K + (I-1)*K + K2
                     IF (M.LT.I2) L4 = (I2-M)*K*K + (K2-1)*K + I
                     MAT(L,L4) = MAT(L,L4) - PP(J,(I2-1)*K+K2)
  120             CONTINUE
  140          CONTINUE
C
  160       CONTINUE
  180    CONTINUE
  200 CONTINUE
C
      K2 = K*K*(P+1)
      IF (P.GT.0) THEN
         IFAIL = 1
         IJAA = 50
         CALL F04ARF(MAT,KMAT,Z,K2,Z,B,IFAIL)
         IF (IFAIL.NE.0) THEN
            IFAILX = 2
            RETURN
         END IF
      END IF
C
C     COPY C(J)'S ONTO GAMMA
C
      M = -1
      DO 260 M1 = 1, P + 1
         M = M + 1
         DO 240 I = 1, K
            DO 220 J = 1, K
               L = M*K*K + (I-1)*K + J
               GAMMA(I,M*K+J) = Z(L)
  220       CONTINUE
  240    CONTINUE
  260 CONTINUE
C
C     BEGINNING OF MAIN LOOP
C
      S = -1
      DO 900 L3 = 1, P
         S = S + 1
C
C        CALCULATE SIGMA AND SIGSTA
C
         DO 340 I = 1, K
            DO 320 J = 1, I
               SUM = GAMMA(I,J)
               SM = SUM
               IF (S.GT.0) THEN
                  DO 300 L = 1, S
                     DO 280 I2 = 1, K
                        SUM = SUM - PHI((S-1)*K+I,(L-1)*K+I2)*GAMMA(I2,
     *                        L*K+J)
                        SM = SM - PHISTA((S-1)*K+I,(L-1)*K+I2)*GAMMA(J,
     *                       L*K+I2)
  280                CONTINUE
  300             CONTINUE
               END IF
               SIGMA(I,J) = SUM
               SIGMA(J,I) = SUM
               SIGSTA(I,J) = SM
               SIGSTA(J,I) = SM
  320       CONTINUE
  340    CONTINUE
C
C        CALCULATE THE CHOLESKI DECOMPOSITIONS OF SIGMA AND SIGSTA
C
         IFAIL = 1
         CALL F03AEF(K,SIGMA,K,B,D1,ID,IFAIL)
         IF (IFAIL.NE.0) THEN
            IFAILX = 2
            RETURN
         END IF
         IFAIL = 1
         CALL F03AEF(K,SIGSTA,K,Z,D1,ID,IFAIL)
         IF (IFAIL.NE.0) THEN
            IFAILX = 2
            RETURN
         END IF
C
C        CALCULATE PHI(S+1,S+1) AND PHISTA(S+1,S+1)
C
C        SET UP MATRIX AA FIRST
C
         DO 420 I = 1, K
            DO 400 J = 1, K
               SUM = GAMMA(J,(S+1)*K+I)
               IF (S.GT.0) THEN
                  DO 380 L = 1, S
                     DO 360 I2 = 1, K
                        SUM = SUM - PHI((S-1)*K+I,(L-1)*K+I2)*GAMMA(J,
     *                        (S-L+1)*K+I2)
  360                CONTINUE
  380             CONTINUE
               END IF
               AA(J,I) = SUM
  400       CONTINUE
  420    CONTINUE
C
         EPS = X02AJF()
         IFAIL = 1
         CALL F04AFF(K,K,SIGSTA,K,Z,AA,K,EPS,BB,K,MAT,KMAT,IK5,IFAIL)
         IF (IFAIL.NE.0) THEN
            IFAILX = 2
            RETURN
         END IF
C
C        TRANSPOSE SOLUTION AND STORE AS PHI(S+1,S+1)
C
         DO 460 I = 1, K
            DO 440 J = 1, K
               PHI(S*K+I,S*K+J) = BB(J,I)
  440       CONTINUE
  460    CONTINUE
C
C        RESET MATRIX AA FOR CALCULATING PHISTA(S+1,S+1)
C
         DO 540 I = 1, K
            DO 520 J = 1, K
               SUM = GAMMA(I,(S+1)*K+J)
               IF (S.GT.0) THEN
                  DO 500 L = 1, S
                     DO 480 I2 = 1, K
                        SUM = SUM - PHISTA((S-1)*K+I,(L-1)*K+I2)
     *                        *GAMMA(I2,(S-L+1)*K+J)
  480                CONTINUE
  500             CONTINUE
               END IF
               AA(J,I) = SUM
  520       CONTINUE
  540    CONTINUE
C
         IFAIL = 1
         CALL F04AFF(K,K,SIGMA,K,B,AA,K,EPS,BB,K,MAT,KMAT,IK5,IFAIL)
         IF (IFAIL.NE.0) THEN
            IFAILX = 2
            RETURN
         END IF
C
C        TRANSPOSE SOLUTION AND STORE AS PHISTA(S+1,S+1)
C
         DO 580 I = 1, K
            DO 560 J = 1, K
               PHISTA(S*K+I,S*K+J) = BB(J,I)
  560       CONTINUE
  580    CONTINUE
C
C        CALCULATE PP
C
C        FIRST FORM PRODUCT PHI(S+1,S+1) * LSTAR AND STORE AS AA
C
C        FIRST SET UP LSTAR
C
         DO 620 I = 1, K
            DO 600 J = 1, K
               IF (I.LT.J) LSTAR(I,J) = 0.0D0
               IF (I.GT.J) LSTAR(I,J) = SIGSTA(I,J)
               IF (I.EQ.J) LSTAR(I,I) = 1.0D0/Z(I)
  600       CONTINUE
  620    CONTINUE
C
         DO 680 I = 1, K
            DO 660 J = 1, K
               SUM = 0.0D0
               DO 640 I2 = 1, K
                  SUM = SUM + PHI(S*K+I,S*K+I2)*LSTAR(I2,J)
  640          CONTINUE
               AA(I,J) = SUM
  660       CONTINUE
  680    CONTINUE
C
         DO 700 I = 1, K
            Z(I) = B(I)
  700    CONTINUE
         DO 760 J = 1, K
            DO 720 I = 1, K
               B(I) = AA(I,J)
  720       CONTINUE
            CALL F04AGZ(SIGMA,K,K,Z,B)
            DO 740 I = 1, K
               AA(I,J) = B(I)
  740       CONTINUE
  760    CONTINUE
C
C        COPY SOLUTION ONTO PP
C
         DO 800 I = 1, K
            DO 780 J = 1, K
               PP(I,S*K+J) = AA(I,J)
  780       CONTINUE
  800    CONTINUE
C
C        UPDATE THE PHI'S AND PHISTA'S
C
         IF (S.EQ.0) GO TO 900
C
         DO 880 L = 1, S
            DO 860 I = 1, K
               DO 840 J = 1, K
                  SUM = PHI((S-1)*K+I,(L-1)*K+J)
                  SUM2 = PHISTA((S-1)*K+I,(L-1)*K+J)
                  DO 820 I2 = 1, K
                     SUM = SUM - PHI(S*K+I,S*K+I2)*PHISTA((S-1)*K+I2,
     *                     (S-L)*K+J)
                     SUM2 = SUM2 - PHISTA(S*K+I,S*K+I2)*PHI((S-1)*K+I2,
     *                      (S-L)*K+J)
  820             CONTINUE
                  PHI(S*K+I,(L-1)*K+J) = SUM
                  PHISTA(S*K+I,(L-1)*K+J) = SUM2
  840          CONTINUE
  860       CONTINUE
  880    CONTINUE
C
  900 CONTINUE
C
C     CONVERT THE PP'S
C
      DO 1120 L = 1, P
C
         DO 960 I = 1, K
            DO 940 J = 1, I
               SUM = 0.0D0
               IF (I.EQ.J) SUM = 1.0D0
               DO 920 I2 = 1, K
                  SUM = SUM - PP(I,(L-1)*K+I2)*PP(J,(L-1)*K+I2)
  920          CONTINUE
               AA(I,J) = SUM
               AA(J,I) = SUM
  940       CONTINUE
  960    CONTINUE
C
         IFAIL = 1
         CALL F03AEF(K,AA,K,Z,D1,ID,IFAIL)
         IF (IFAIL.NE.0) THEN
            IFAILX = 2
            RETURN
         END IF
C
C        COPY PP ONTO L2
C
         DO 1000 I = 1, K
            DO 980 J = 1, K
               L2(I,J) = PP(I,(L-1)*K+J)
  980       CONTINUE
 1000    CONTINUE
C
         DO 1060 J = 1, K
            DO 1020 I = 1, K
               B(I) = L2(I,J)
 1020       CONTINUE
            CALL F04AGZ(AA,K,K,Z,B)
            DO 1040 I = 1, K
               L2(I,J) = B(I)
 1040       CONTINUE
 1060    CONTINUE
C
C        COPY L2 ONTO XX ARRAY
C
         DO 1100 I = 1, K
            DO 1080 J = 1, K
               XX(I6+(L-1)*K*K+(I-1)*K+J) = L2(I,J)
 1080       CONTINUE
 1100    CONTINUE
C
 1120 CONTINUE
C
      RETURN
      END
