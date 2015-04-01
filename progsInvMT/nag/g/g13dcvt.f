      SUBROUTINE G13DCV(X,P,PP,K,N2,SIG,PHI,TEMP,A,B,C,L2,LSTAR,W,SIGMA,
     *                  SIGSTA,PHISTA,KR,I9,IFAILX)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     THIS ROUTINE IMPLEMENTS THE REPARAMETERISATION TO ENFORCE
C     STATIONARITY AND INVERTIBILITY ON THE ARMA MODEL WHEN ALL
C     THE PARAMETERS ARE FREE FROM CONSTRAINTS
C
C     C.F. ANSLEY AND R. KOHN (1986). JOURNAL OF STATISTICAL
C     COMPUTATION AND SIMULATION. ' A NOTE ON REPARAMETERISING A
C     VECTOR AUTOREGRESSIVE MOVING AVERAGE MODEL TO ENFORCE
C     STATIONARITY '.
C
C     .. Scalar Arguments ..
      INTEGER           I9, IFAILX, K, KR, N2, P
C     .. Array Arguments ..
      DOUBLE PRECISION  A(K,K), B(K,K), C(K,K), L2(K,K), LSTAR(K,K),
     *                  PHI(KR,P*K), PHISTA(KR,P*K), PP(K,P*K+1),
     *                  SIG(K,K), SIGMA(K,K), SIGSTA(K,K), TEMP(K,K),
     *                  W(K), X(N2)
C     .. Local Scalars ..
      DOUBLE PRECISION  D1, SUM, SUM2
      INTEGER           I, I2, ID, IFAIL, J, J2, L, S
C     .. External Subroutines ..
      EXTERNAL          F03AEF, F04AGZ
C     .. Executable Statements ..
C
C     SET UP PP MATRICES
C
      IFAILX = 0
      DO 160 L = 1, P
C
C        CALCULATE A(L)
C
         DO 40 I = 1, K
            DO 20 J = 1, K
               A(I,J) = X(I9+(L-1)*K*K+(I-1)*K+J)
   20       CONTINUE
   40    CONTINUE
C
C        CONSTRUCT B(L)
C
         DO 100 I = 1, K
            DO 80 J = 1, K
               IF (I.EQ.J) THEN
                  SUM = 1.0D0
               ELSE
                  SUM = 0.0D0
               END IF
               DO 60 I2 = 1, K
                  SUM = SUM + A(I,I2)*A(J,I2)
   60          CONTINUE
               TEMP(I,J) = SUM
   80       CONTINUE
  100    CONTINUE
C
         IFAIL = 1
         CALL F03AEF(K,TEMP,K,W,D1,ID,IFAIL)
         IF (IFAIL.NE.0) THEN
            IFAILX = 6
            RETURN
         END IF
C
C        NOW SET UP PP(L)
C
         DO 140 J = 1, K
            CALL F04AGZ(TEMP,K,K,W,A(1,J))
            DO 120 I = 1, K
               PP(I,(L-1)*K+J) = A(I,J)
  120       CONTINUE
  140    CONTINUE
C
  160 CONTINUE
C
C     INITIALISE L2,LSTAR,SIGMA,SIGSTA
C
      DO 200 I = 1, K
         DO 180 J = 1, K
            IF (I.EQ.J) THEN
               L2(I,J) = 1.0D0
               SIGMA(I,J) = 1.0D0
               SIGSTA(I,J) = 1.0D0
               LSTAR(I,J) = 1.0D0
            ELSE
               L2(I,J) = 0.0D0
               SIGMA(I,J) = 0.0D0
               SIGSTA(I,J) = 0.0D0
               LSTAR(I,J) = 0.0D0
            END IF
  180    CONTINUE
  200 CONTINUE
C
      DO 220 I = 1, K
         SIG(I,1) = 1.0D0
         W(I) = 1.0D0
  220 CONTINUE
C
      DO 820 S = 0, P - 1
C
C        CALCULATE PHI(S+1,S+1)
C
         DO 260 I = 1, K
            DO 240 J = 1, K
               C(I,J) = LSTAR(I,J)
  240       CONTINUE
  260    CONTINUE
         DO 320 J = 1, K
            DO 280 I = 1, K
               A(I,1) = 0.0D0
  280       CONTINUE
            A(J,1) = 1.0D0
            CALL F04AGZ(SIGSTA,K,K,SIG(1,1),A(1,1))
            DO 300 I = 1, K
               TEMP(I,J) = A(I,1)
  300       CONTINUE
  320    CONTINUE
C
         DO 400 I = 1, K
            DO 380 J = 1, K
               SUM = 0.0D0
               DO 360 I2 = 1, K
                  DO 340 J2 = 1, K
                     SUM = SUM + L2(I,I2)*PP(I2,S*K+J2)*TEMP(J2,J)
  340             CONTINUE
  360          CONTINUE
               PHI(S*K+I,S*K+J) = SUM
  380       CONTINUE
  400    CONTINUE
C
C        CALCULATE PHISTA(S+1,S+1)
C
         DO 460 J = 1, K
            DO 420 I = 1, K
               A(I,1) = 0.0D0
  420       CONTINUE
            A(J,1) = 1.0D0
            CALL F04AGZ(SIGMA,K,K,W,A(1,1))
            DO 440 I = 1, K
               TEMP(I,J) = A(I,1)
  440       CONTINUE
  460    CONTINUE
C
         DO 540 I = 1, K
            DO 520 J = 1, K
               SUM = 0.0D0
               DO 500 I2 = 1, K
                  DO 480 J2 = 1, K
                     SUM = SUM + C(I,I2)*PP(J2,S*K+I2)*TEMP(J2,J)
  480             CONTINUE
  500          CONTINUE
               PHISTA(S*K+I,S*K+J) = SUM
  520       CONTINUE
  540    CONTINUE
C
C        UPDATE THE PHI'S AND PHISTA'S
C
         IF (S.EQ.0) GO TO 640
C
         DO 620 L = 1, S
            DO 600 I = 1, K
               DO 580 J = 1, K
                  SUM = PHI((S-1)*K+I,(L-1)*K+J)
                  SUM2 = PHISTA((S-1)*K+I,(L-1)*K+J)
                  DO 560 I2 = 1, K
                     SUM = SUM - PHI(S*K+I,S*K+I2)*PHISTA((S-1)*K+I2,
     *                     (S-L)*K+J)
                     SUM2 = SUM2 - PHISTA(S*K+I,S*K+I2)*PHI((S-1)*K+I2,
     *                      (S-L)*K+J)
  560             CONTINUE
                  PHI(S*K+I,(L-1)*K+J) = SUM
                  PHISTA(S*K+I,(L-1)*K+J) = SUM2
  580          CONTINUE
  600       CONTINUE
  620    CONTINUE
C
C        UPDATE SIGMA AND SIGSTA
C
  640    DO 680 I = 1, K
            DO 660 J = I, K
               TEMP(I,J) = SIGMA(I,J)
               C(I,J) = SIGSTA(I,J)
               TEMP(J,I) = SIGMA(I,J)
               C(J,I) = SIGSTA(I,J)
  660       CONTINUE
  680    CONTINUE
C
         DO 760 I = 1, K
            DO 740 J = 1, K
               SUM = TEMP(I,J)
               SUM2 = C(I,J)
               DO 720 I2 = 1, K
                  DO 700 J2 = 1, K
                     SUM = SUM - PHI(S*K+I,S*K+I2)*C(I2,J2)*PHI(S*K+J,
     *                     S*K+J2)
                     SUM2 = SUM2 - PHISTA(S*K+I,S*K+I2)*TEMP(I2,J2)
     *                      *PHISTA(S*K+J,S*K+J2)
  700             CONTINUE
  720          CONTINUE
C
               SIGMA(I,J) = SUM
               SIGSTA(I,J) = SUM2
  740       CONTINUE
  760    CONTINUE
C
C        CALCULATE CHOLESKI FACTORIZATIONS OF SIGMA AND SIGSTA
C
         IFAIL = 1
         CALL F03AEF(K,SIGMA,K,W,D1,ID,IFAIL)
         IF (IFAIL.NE.0) THEN
            IFAILX = 6
            RETURN
         END IF
         IFAIL = 1
         CALL F03AEF(K,SIGSTA,K,SIG(1,1),D1,ID,IFAIL)
         IF (IFAIL.NE.0) THEN
            IFAILX = 6
            RETURN
         END IF
C
C        SET L2 AND LSTAR
C
         DO 800 I = 1, K
            DO 780 J = 1, I
               IF (I.NE.J) THEN
                  L2(I,J) = SIGMA(I,J)
                  LSTAR(I,J) = SIGSTA(I,J)
                  L2(J,I) = 0.0D0
                  LSTAR(J,I) = 0.0D0
               ELSE
                  L2(I,I) = 1.0D0/W(I)
                  LSTAR(I,I) = 1.0D0/SIG(I,1)
               END IF
  780       CONTINUE
  800    CONTINUE
C
  820 CONTINUE
C
C     STANDARDISE THE PHI'S
C
      DO 940 L = 1, P
C
         DO 880 I = 1, K
            DO 860 J = 1, K
               SUM = 0.0D0
               DO 840 I2 = 1, K
                  SUM = SUM + PHI((P-1)*K+I,(L-1)*K+I2)*L2(I2,J)
  840          CONTINUE
               A(I,J) = SUM
  860       CONTINUE
  880    CONTINUE
C
         DO 920 J = 1, K
            CALL F04AGZ(SIGMA,K,K,W,A(1,J))
            DO 900 I = 1, K
               PP(I,(L-1)*K+J) = A(I,J)
  900       CONTINUE
  920    CONTINUE
C
  940 CONTINUE
C
C     RESET SIG TO THE IDENTITY MATRIX
C
      DO 980 I = 1, K
         DO 960 J = 1, K
            SIG(I,J) = 0.0D0
  960    CONTINUE
         SIG(I,I) = 1.0D0
  980 CONTINUE
C
      RETURN
      END
