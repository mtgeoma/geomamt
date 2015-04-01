      SUBROUTINE G13DCY(X,N2,PHI,THETA,GAMMA,TEMP,TEMPL,W,V,F,Z,MAT,
     *                  INVF,TEMPM,TEMPK,QQ,A,B,MT,MEANS,REZ,IMZ,W3,PA,
     *                  PB,INTGR,SUM,IFLAG)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C     SUBROUTINE NFUNCT CALCULATES THE LOG LIKELIHOOD FUNCTION
C     OF A VECTOR AUTOREGRESSIVE - MOVING AVERAGE MODEL USING
C     A KALMAN FILTER ALGORITHM DISCUSSED BY B.L. SHEA (1986) IN
C     JOURNAL OF TIME SERIES ANALYSIS,' ESTIMATION OF
C     MULTIVARIATE TIME SERIES'.
C
C     USER - SUPPLIED OPTION :
C
C     CONDS = .FALSE. (EXACT LIKELIHOOD)
C             .TRUE.  (CONDITIONAL LIKELIHOOD)
C
C
C     CORRECT W FOR THE MEAN AND PUT IN W3
C
C     .. Scalars in Common ..
      DOUBLE PRECISION  ADDLOG, CONDD, DETQQ, EXPP, FBIG, LMAX, NORM,
     *                  SMALL, SSQ, SUMSSQ, XTOL
      INTEGER           IFAILX, INTEG, K, K3, K4, K5, K6, KMAT, KR, KZ,
     *                  LEW6, LEW7, LP, LW1, LW11, LW12, LW13, LW14,
     *                  LW16, LW17, LW18, LW2, LW3, LW4, LW5, LW8, N,
     *                  NITER, P, Q, R
      LOGICAL           ANOTT, CONDS, FULLP, FULLQ, MEAN, NOPRIN
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SUM
      INTEGER           IFLAG, N2
C     .. Array Arguments ..
      DOUBLE PRECISION  A(K+1,K), B(K*K*(P+1)), F(K,K), GAMMA(K,KR),
     *                  IMZ(KR), INVF(K+1,K), MAT(KMAT,KMAT), MEANS(K),
     *                  MT(K,K), PA(KR,KR), PB(KR,KR), PHI(K,P*K+1),
     *                  QQ(K,K), REZ(KR), TEMP(K,KR), TEMPK(K,KR),
     *                  TEMPL(K,K), TEMPM(K,K), THETA(K,Q*K+1), V(K,N),
     *                  W(K,N), W3(K,N), X(N2), Z(KZ)
      INTEGER           INTGR(KR)
C     .. Local Scalars ..
      DOUBLE PRECISION  DETP, SIG, SM, TSIG
      INTEGER           I, I2, I3, I9, IFAIL, IFAILW, J, J2, K1, K2,
     *                  K44, KW, L, L8, M, T
      LOGICAL           DELTA, INVERT, PEARL, STAT
C     .. External Subroutines ..
      EXTERNAL          F01ADF, F03ABF, F04ARF, G13DCM, G13DCN, G13DCV,
     *                  G13DCW, G13DCX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, LOG, MIN, DBLE
C     .. Common blocks ..
      COMMON            /BG13DC/ADDLOG, LMAX, CONDD, NORM, K, P, Q, K4,
     *                  K5, LW14, LW18, LP, K6, NITER, LEW6, LEW7, LW2,
     *                  LW1, LW3, LW4, LW5, LW8, LW11, LW12, LW13, LW16,
     *                  LW17, KR, MEAN, NOPRIN, FULLP, FULLQ
      COMMON            /CG13DC/XTOL, SSQ, SUMSSQ, EXPP, DETQQ, FBIG,
     *                  SMALL, N, R, INTEG, K3, KZ, KMAT, IFAILX, ANOTT,
     *                  CONDS
C     .. Executable Statements ..
      IF (MEAN) THEN
         DO 40 J = 1, N
            DO 20 I = 1, K
               W3(I,J) = W(I,J) - MEANS(I) - X(K5+I)
   20       CONTINUE
   40    CONTINUE
      ELSE
         DO 80 J = 1, N
            DO 60 I = 1, K
               W3(I,J) = W(I,J)
   60       CONTINUE
   80    CONTINUE
      END IF
C
C     SET QQ TO THE IDENTITY MATRIX
C
      DO 120 I = 1, K
         DO 100 J = 1, K
            QQ(I,J) = 0.0D0
            IF (I.EQ.J) QQ(I,J) = 1.0D0
  100    CONTINUE
  120 CONTINUE
C
C      CALCULATE ORIGINAL AR PARAMETERS AND PUT IN PHI ARRAY
C
      IF (P.EQ.0) GO TO 200
      IF (FULLP) THEN
         I9 = 0
         CALL G13DCV(X,P,PHI,K,N2,QQ,PA,GAMMA,TEMP,TEMPL,V,F,MT,B,TEMPM,
     *               TEMPK,PB,KR,I9,IFAILW)
         IF (IFAILW.NE.0) THEN
            IF (INTEG.EQ.1) SUM = FBIG
            IF (INTEG.NE.1) IFLAG = -1
            RETURN
         END IF
      ELSE
C
         DO 180 L = 1, P
            DO 160 I = 1, K
               DO 140 J = 1, K
                  PHI(I,(L-1)*K+J) = X((L-1)*K3+(I-1)*K+J)
  140          CONTINUE
  160       CONTINUE
  180    CONTINUE
C
      END IF
C
C      CALCULATE ORIGINAL MA PARAMETERS AND PUT IN THETA ARRAY
C
  200 IF (Q.EQ.0) GO TO 280
      IF (FULLQ) THEN
         CALL G13DCV(X,Q,THETA,K,N2,QQ,PA,GAMMA,TEMP,TEMPL,V,F,MT,B,
     *               TEMPM,TEMPK,PB,KR,K4,IFAILW)
         IF (IFAILW.NE.0) THEN
            IF (INTEG.EQ.1) SUM = FBIG
            IF (INTEG.NE.1) IFLAG = -1
            RETURN
         END IF
      ELSE
C
         DO 260 L = 1, Q
            DO 240 I = 1, K
               DO 220 J = 1, K
                  THETA(I,(L-1)*K+J) = X(K4+(L-1)*K3+(I-1)*K+J)
  220          CONTINUE
  240       CONTINUE
  260    CONTINUE
C
      END IF
C
C     CONSTRUCT SIGMA MATRIX
C
  280 DO 340 I = 1, K
         DO 320 J = 1, I
            QQ(I,J) = 0.0D0
            DO 300 I3 = 1, K
               IF (MIN(I,J).GE.I3) QQ(I,J) = QQ(I,J) + X(K6+(I-1)
     *             *I/2+I3)*X(K6+(J-1)*J/2+I3)
  300       CONTINUE
            QQ(I,J) = LMAX*LMAX*QQ(I,J)
            QQ(J,I) = QQ(I,J)
C
  320    CONTINUE
  340 CONTINUE
C
C      SET TSIG TO THE DETERMINANT OF SIGMA
C
      TSIG = X(K6+1)
      IF (K.EQ.1) GO TO 380
      DO 360 I = 2, K
         TSIG = TSIG*X(K6+(I+1)*I/2)
  360 CONTINUE
  380 TSIG = (TSIG**2)*(LMAX**(2*K))
C
C     TEST FOR STATIONARITY
C
      STAT = .TRUE.
      IF (P.GT.0) CALL G13DCX(P,K,PHI,MAT,KMAT,R,REZ,IMZ,INTGR,KR,STAT)
C
C     TEST FOR INVERTIBILITY
C
      INVERT = .TRUE.
      IF (Q.GT.0) CALL G13DCX(Q,K,THETA,MAT,KMAT,R,REZ,IMZ,INTGR,KR,
     *                        INVERT)
C
C     IF THE CURRENT VALUE OF X IS NOT WITHIN THE 'ADMISSABILITY' REGION
C     THEN WE SET THE LOG LIKELIHOOD EQUAL TO A LARGE VALUE
C
      IF (STAT .AND. INVERT) GO TO 400
      IF (INTEG.GE.2) THEN
         IFLAG = -1
         RETURN
      END IF
      IF (INTEG.EQ.0) THEN
         IFAILX = 2
         RETURN
      END IF
      SUM = FBIG
      RETURN
C
  400 IF (CONDS) THEN
         DO 440 I = 1, K
            DO 420 J = 1, K
               TEMPK(I,J) = 0.0D0
  420       CONTINUE
  440    CONTINUE
         GO TO 700
      END IF
C
      IF (P.EQ.0) THEN
         CALL G13DCN(K,Q,THETA,QQ,KR,TEMPK,TEMP)
         GO TO 700
      END IF
C
      KW = K*(K+1)/2
      IF (P.GE.2) KW = KW + (P-1)*K*K
      CALL G13DCM(K,P,Q,PHI,THETA,QQ,GAMMA,KR,Z,KW,MAT,B,TEMP(1,1),
     *            IFAIL)
      IF (IFAIL.NE.0) THEN
         IF (INTEG.EQ.1) SUM = FBIG
         IF (INTEG.NE.1) IFLAG = -1
         RETURN
      END IF
C
C      CALCULATE FIRST K COLUMNS OF P(1/0) AND STORE AS TEMPK.
C
      K44 = K*(K+1)/2
      DO 580 K1 = 1, P
         DO 560 I = 1, K
            DO 540 J = 1, K
               IF (K1.EQ.1) THEN
                  IF (I.GE.J) THEN
                     SUM = Z(I*(I-1)/2+J)
                  ELSE
                     SUM = Z(J*(J-1)/2+I)
                  END IF
               ELSE
                  SUM = Z(K44+(K1-2)*K*K+(J-1)*K+I)
                  DO 460 K2 = 1, K
                     IF (K2.GE.J) THEN
                        SUM = SUM - PHI(I,(K1-2)*K+K2)*Z(K2*(K2-1)/2+J)
                     ELSE
                        SUM = SUM - PHI(I,(K1-2)*K+K2)*Z(J*(J-1)/2+K2)
                     END IF
  460             CONTINUE
               END IF
               DO 500 M = 1, K1 - 2
                  DO 480 K2 = 1, K
                     SUM = SUM - PHI(I,(M-1)*K+K2)*Z(K44+(K1-2-M)
     *                     *K*K+(J-1)*K+K2)
  480             CONTINUE
C
  500          CONTINUE
C
               IF (K1.EQ.1) THEN
                  SUM = SUM - QQ(I,J)
               ELSE IF ((K1.GT.1) .AND. (K1.LE.Q+1)) THEN
                  DO 520 K2 = 1, K
                     SUM = SUM + THETA(I,(K1-2)*K+K2)*QQ(K2,J)
  520             CONTINUE
               END IF
C
               TEMPK(I,(K1-1)*K+J) = SUM
  540       CONTINUE
  560    CONTINUE
  580 CONTINUE
C
      DO 680 K1 = P + 1, R
         DO 660 I = 1, K
            DO 640 J = 1, K
C
               SUM = 0.0D0
               DO 620 M = K1, Q
                  DO 600 K2 = 1, K
                     SUM = SUM - THETA(I,(M-1)*K+K2)*GAMMA(J,(M-K1+1)
     *                     *K+K2)
  600             CONTINUE
  620          CONTINUE
C
               TEMPK(I,(K1-1)*K+J) = SUM
  640       CONTINUE
  660    CONTINUE
  680 CONTINUE
C
C     INITIALISE A(1/0),V(1),SSQ,F(1) AND DETP
C
  700 DO 720 I = 1, KR
         Z(I) = 0.0D0
  720 CONTINUE
C
      DO 760 I = 1, K
         DO 740 J = 1, K
            F(I,J) = TEMPK(I,J) + QQ(I,J)
  740    CONTINUE
         V(I,1) = W3(I,1)
  760 CONTINUE
      DO 800 I = 1, K
         DO 780 J = 1, K
            A(I,J) = F(I,J)
  780    CONTINUE
  800 CONTINUE
C
C     CALCULATE CHOLESKI DECOMPOSITION OF QQ IF XTOL .LT. 0.0
C
      IF ((INTEG.EQ.2) .AND. ( .NOT. CONDS)) THEN
         CALL G13DCW(QQ,K,K,STAT,F)
         DO 840 I = 1, K
            DO 820 J = 1, K
               PA(I,J) = F(I,J)
  820       CONTINUE
  840    CONTINUE
         DO 880 I = 1, K
            DO 860 J = 1, K
               F(I,J) = A(I,J)
  860       CONTINUE
  880    CONTINUE
      END IF
C
      IFAIL = 1
      CALL F03ABF(A,K+1,K,DETP,B,IFAIL)
      IF (IFAIL.NE.0) THEN
         IF (INTEG.EQ.1) SUM = FBIG
         IF (INTEG.NE.1) IFLAG = -1
         RETURN
      END IF
C
      DO 920 I = 1, K
         DO 900 J = 1, K
            A(I,J) = F(I,J)
  900    CONTINUE
  920 CONTINUE
C
      IFAIL = 1
      CALL F01ADF(K,A,K+1,IFAIL)
      IF (IFAIL.NE.0) THEN
         IF (INTEG.EQ.1) SUM = FBIG
         IF (INTEG.NE.1) IFLAG = -1
         RETURN
      END IF
C
      DO 960 I = 1, K
         DO 940 J = 1, I
            INVF(I,J) = A(I+1,J)
            INVF(J,I) = INVF(I,J)
  940    CONTINUE
  960 CONTINUE
C
      SSQ = 0.0D0
      DO 1000 I = 1, K
         DO 980 J = 1, K
            SSQ = SSQ + V(I,1)*INVF(I,J)*V(J,1)
  980    CONTINUE
 1000 CONTINUE
C
      IF ((INTEG.EQ.2) .AND. ( .NOT. CONDS)) THEN
C
C        CALCULATE CHOLESKI DECOMPOSITION OF F(1)
C
C        FIRST COPY F ONTO MT
C
         DO 1020 I = 1, K
            IMZ(I) = V(I,1)
 1020    CONTINUE
         CALL G13DCW(F,K,K,STAT,MT)
         IFAIL = 1
         CALL F04ARF(MT,K,IMZ,K,IMZ,REZ,IFAIL)
C
C        PUT STANDARDIZED V(T) INTO REZ
C
         DO 1060 I = 1, K
            SUM = 0.0D0
            DO 1040 I2 = 1, K
               SUM = SUM + PA(I,I2)*IMZ(I2)
 1040       CONTINUE
            REZ(I) = SUM
 1060    CONTINUE
C
      END IF
C
      DETP = LOG(DETP/TSIG)
      IF (CONDS) GO TO 1460
C
C     INITIALISE M(1) TO - F(1) INVERSE
C
      DO 1100 I = 1, K
         DO 1080 J = 1, K
            MT(I,J) = -INVF(I,J)
 1080    CONTINUE
 1100 CONTINUE
C
C
C     CALCULATE L(1) AND M(1)
C
      DO 1140 I = 1, K
         DO 1120 J = 1, KR
            MAT(I,J) = 0.0D0
            GAMMA(I,J) = 0.0D0
 1120    CONTINUE
 1140 CONTINUE
C
      DO 1260 L = 1, R
         DO 1240 I = 1, K
            DO 1220 J = 1, K
               SUM = 0.0D0
               IF (L.GT.P) GO TO 1180
               DO 1160 K2 = 1, K
                  SUM = SUM + PHI(I,(L-1)*K+K2)*TEMPK(K2,J)
 1160          CONTINUE
 1180          IF (L.LT.R) SUM = SUM + TEMPK(I,L*K+J)
               DO 1200 K2 = 1, K
                  SM = 0.0D0
                  IF (L.LE.P) SM = PHI(I,(L-1)*K+K2)
                  IF (L.LE.Q) SM = SM - THETA(I,(L-1)*K+K2)
                  SUM = SUM + SM*QQ(K2,J)
 1200          CONTINUE
               MAT(I,(L-1)*K+J) = SUM
C
 1220       CONTINUE
 1240    CONTINUE
 1260 CONTINUE
C
      IF (ANOTT) GO TO 1340
C
      DO 1320 L = 1, R
         DO 1300 I = 1, K
            DO 1280 J = 1, K
               GAMMA(I,(L-1)*K+J) = MAT(I,(L-1)*K+J)
 1280       CONTINUE
 1300    CONTINUE
 1320 CONTINUE
C
      GO TO 1460
 1340 DO 1440 L = 1, R
         DO 1420 I = 1, K
            DO 1400 J = 1, K
               SUM = 0.0D0
               IF (L.GT.Q) GO TO 1380
               DO 1360 K2 = 1, K
                  SUM = SUM + THETA(I,(L-1)*K+K2)*TEMPK(K2,J)
 1360          CONTINUE
 1380          IF (L.LT.R) SUM = SUM + TEMPK(I,L*K+J)
               GAMMA(I,(L-1)*K+J) = SUM
 1400       CONTINUE
 1420    CONTINUE
 1440 CONTINUE
C
C     START THE RECURSIONS
C
 1460 DELTA = .FALSE.
      IF (CONDS) THEN
C
         DELTA = .TRUE.
C
C        SET TEMPK = R
C
         DO 1520 I = 1, K
            DO 1500 J = 1, K
               DO 1480 L = 1, R
                  SUM = 0.0D0
                  IF (L.LE.P) SUM = PHI(I,(L-1)*K+J)
                  IF (L.LE.Q) SUM = SUM - THETA(I,(L-1)*K+J)
                  TEMPK(I,(L-1)*K+J) = SUM
 1480          CONTINUE
 1500       CONTINUE
 1520    CONTINUE
C
      END IF
      PEARL = .FALSE.
C
      DO 3240 T = 2, N
C
         IF (ANOTT .AND. (T.GT.P-Q)) PEARL = .TRUE.
C
C        CALCULATE TEMPK
C
         IF (DELTA) GO TO 1620
C
         DO 1600 L = 1, R
            DO 1580 I = 1, K
               DO 1560 J = 1, K
                  SUM = 0.0D0
                  DO 1540 K2 = 1, K
                     SUM = SUM + GAMMA(I,(L-1)*K+K2)*INVF(K2,J)
 1540             CONTINUE
                  TEMPK(I,(L-1)*K+J) = SUM
 1560          CONTINUE
 1580       CONTINUE
 1600    CONTINUE
C
C        CALCULATE TEMP
C
 1620    DO 1640 L = 1, KR
            TEMP(1,L) = 0.0D0
 1640    CONTINUE
C
         IF (ANOTT) GO TO 1740
C
C        TEMP = T*A(T-1/T-2)
C
         IF (T.EQ.2) GO TO 1920
         DO 1720 L = 1, R
            DO 1700 I = 1, K
               SUM = 0.0D0
               IF (L.GT.P) GO TO 1680
               DO 1660 K2 = 1, K
                  SUM = SUM + PHI(I,(L-1)*K+K2)*Z(K2)
 1660          CONTINUE
 1680          IF (L.LT.R) SUM = SUM + Z(L*K+I)
               TEMP(1,(L-1)*K+I) = SUM
 1700       CONTINUE
 1720    CONTINUE
         GO TO 1920
 1740    IF (T.EQ.2) GO TO 1840
C
C        TEMP = A*A(T-1/T-2) + R*W(T-1)
C
         DO 1820 L = 1, R
            DO 1800 I = 1, K
               SUM = 0.0D0
               IF (L.GT.Q) GO TO 1780
               DO 1760 K2 = 1, K
                  SUM = SUM + THETA(I,(L-1)*K+K2)*Z(K2)
 1760          CONTINUE
 1780          IF (L.LT.R) SUM = SUM + Z(L*K+I)
               TEMP(1,(L-1)*K+I) = SUM
 1800       CONTINUE
 1820    CONTINUE
C
 1840    DO 1900 L = 1, R
            DO 1880 I = 1, K
               SUM = 0.0D0
               DO 1860 K2 = 1, K
                  SM = 0.0D0
                  IF (L.LE.P) SM = PHI(I,(L-1)*K+K2)
                  IF (L.LE.Q) SM = SM - THETA(I,(L-1)*K+K2)
                  SUM = SUM + SM*W3(K2,T-1)
 1860          CONTINUE
               TEMP(1,(L-1)*K+I) = TEMP(1,(L-1)*K+I) + SUM
 1880       CONTINUE
 1900    CONTINUE
C
C        A(T/T-1) = TEMP+TEMPK*V(T-1)
C
 1920    DO 2000 L = 1, R
            DO 1980 L8 = 1, K
               SUM = TEMP(1,(L-1)*K+L8)
               IF (DELTA .AND. ANOTT) GO TO 1960
               DO 1940 K2 = 1, K
                  SUM = SUM + TEMPK(L8,(L-1)*K+K2)*V(K2,T-1)
 1940          CONTINUE
 1960          Z((L-1)*K+L8) = SUM
 1980       CONTINUE
 2000    CONTINUE
C
         IF ((INTEG.EQ.2) .AND. ( .NOT. CONDS)) THEN
            DO 2020 I = 1, K
               V(I,T-1) = REZ(I)
 2020       CONTINUE
         END IF
C
         IF (DELTA) GO TO 3080
C
C        CALCULATE TEMPL
C
         DO 2060 I = 1, K
            DO 2040 J = 1, K
               TEMPL(I,J) = MAT(I,J)
 2040       CONTINUE
 2060    CONTINUE
C
C        RECALCULATE TEMP
C
         IF (ANOTT) GO TO 2180
C
C        TEMP = T*L(T-1)
C
         DO 2160 L = 1, R
            DO 2140 I = 1, K
               DO 2120 J = 1, K
                  SUM = 0.0D0
                  IF (L.GT.P) GO TO 2100
                  DO 2080 K2 = 1, K
                     SUM = SUM + PHI(I,(L-1)*K+K2)*MAT(K2,J)
 2080             CONTINUE
 2100             IF (L.LT.R) SUM = SUM + MAT(I,L*K+J)
                  TEMP(I,(L-1)*K+J) = SUM
C
 2120          CONTINUE
 2140       CONTINUE
 2160    CONTINUE
C
         GO TO 2300
C
C        TEMP = A*L(T-1)
C
 2180    DO 2280 L = 1, R
            DO 2260 I = 1, K
               DO 2240 J = 1, K
                  SUM = 0.0D0
                  IF (L.GT.Q) GO TO 2220
                  DO 2200 K2 = 1, K
                     SUM = SUM + THETA(I,(L-1)*K+K2)*MAT(K2,J)
 2200             CONTINUE
 2220             IF (L.LT.R) SUM = SUM + MAT(I,L*K+J)
                  TEMP(I,(L-1)*K+J) = SUM
C
 2240          CONTINUE
 2260       CONTINUE
 2280    CONTINUE
C
C        CALCULATE TEMPM
C
 2300    DO 2360 I = 1, K
            DO 2340 J = 1, K
               TEMPM(I,J) = 0.0D0
               DO 2320 K2 = 1, K
                  TEMPM(I,J) = TEMPM(I,J) + MT(I,K2)*TEMPL(J,K2)
 2320          CONTINUE
 2340       CONTINUE
 2360    CONTINUE
C
C        UPDATE K(T)
C
         DO 2440 L = 1, R
            DO 2420 I = 1, K
               DO 2400 J = 1, K
                  SUM = 0.0D0
                  IF (PEARL .AND. (L.GE.Q+1)) GAMMA(I,(L-1)*K+J) = SUM
                  IF (PEARL .AND. (L.GE.Q+1)) GO TO 2400
                  DO 2380 K2 = 1, K
                     SUM = SUM + TEMP(I,(L-1)*K+K2)*TEMPM(K2,J)
 2380             CONTINUE
                  GAMMA(I,(L-1)*K+J) = GAMMA(I,(L-1)*K+J) + SUM
 2400          CONTINUE
 2420       CONTINUE
 2440    CONTINUE
C
C        UPDATE F(T)
C
         DO 2500 I = 1, K
            DO 2480 J = I, K
               DO 2460 K2 = 1, K
                  F(I,J) = F(I,J) + TEMPL(I,K2)*TEMPM(K2,J)
 2460          CONTINUE
               F(J,I) = F(I,J)
 2480       CONTINUE
 2500    CONTINUE
C
C        TEST FOR CONVERGENCE OF F(T)'S
C
         SUM = 0.0D0
         DO 2540 I = 1, K
            DO 2520 J = 1, I
               IF (ABS(QQ(I,J)).GT.SMALL) THEN
                  SM = ABS((F(I,J)-QQ(I,J))/QQ(I,J))
               ELSE
                  SM = ABS(F(I,J)-QQ(I,J))
               END IF
               IF (SM.GT.SUM) SUM = SM
 2520       CONTINUE
 2540    CONTINUE
         IF (SUM.LT.XTOL) DELTA = .TRUE.
C
         IF ( .NOT. DELTA) GO TO 2700
C
C        SET TEMPK = R AND INVF = SIGMA INVERSE
C
         DO 2600 I = 1, K
            DO 2580 J = 1, K
               DO 2560 L = 1, R
                  SUM = 0.0D0
                  IF (L.LE.P) SUM = PHI(I,(L-1)*K+J)
                  IF (L.LE.Q) SUM = SUM - THETA(I,(L-1)*K+J)
                  TEMPK(I,(L-1)*K+J) = SUM
 2560          CONTINUE
 2580       CONTINUE
 2600    CONTINUE
C
         IFAIL = 1
         DO 2640 I = 1, K
            DO 2620 J = 1, K
               A(I,J) = QQ(I,J)
 2620       CONTINUE
 2640    CONTINUE
         CALL F01ADF(K,A,K+1,IFAIL)
         IF (IFAIL.NE.0) THEN
            IF (INTEG.EQ.1) SUM = FBIG
            IF (INTEG.NE.1) IFLAG = -1
            RETURN
         END IF
C
         DO 2680 I = 1, K
            DO 2660 J = 1, I
               INVF(I,J) = A(I+1,J)
               INVF(J,I) = INVF(I,J)
 2660       CONTINUE
 2680    CONTINUE
C
         IF (DELTA) GO TO 3080
C
C        CALCULATE INVERSE OF F(T) AND DET(F(T))
C
 2700    DO 2740 I = 1, K
            DO 2720 J = 1, K
               A(I,J) = F(I,J)
 2720       CONTINUE
 2740    CONTINUE
C
         IFAIL = 1
         CALL F03ABF(A,K+1,K,SIG,B,IFAIL)
         IF (IFAIL.NE.0) THEN
            IF (INTEG.EQ.1) SUM = FBIG
            IF (INTEG.NE.1) IFLAG = -1
            RETURN
         END IF
C
         DO 2780 I = 1, K
            DO 2760 J = 1, K
               A(I,J) = F(I,J)
 2760       CONTINUE
 2780    CONTINUE
C
         IFAIL = 1
         CALL F01ADF(K,A,K+1,IFAIL)
         IF (IFAIL.NE.0) THEN
            IF (INTEG.EQ.1) SUM = FBIG
            IF (INTEG.NE.1) IFLAG = -1
            RETURN
         END IF
         DO 2820 I = 1, K
            DO 2800 J = 1, I
               INVF(I,J) = A(I+1,J)
               INVF(J,I) = INVF(I,J)
 2800       CONTINUE
 2820    CONTINUE
C
         IF ((INTEG.EQ.2) .AND. ( .NOT. CONDS)) THEN
C
C           CALCULATE CHOLESKI DECOMPOSITION OF F(T)
C
            DO 2860 I = 1, K
               DO 2840 J = 1, K
                  PB(I,J) = INVF(I,J)
                  INVF(I,J) = F(I,J)
 2840          CONTINUE
 2860       CONTINUE
C
            CALL G13DCW(INVF,K+1,K,STAT,A)
C
            DO 2900 I = 1, K
               DO 2880 J = 1, K
                  INVF(I,J) = PB(I,J)
 2880          CONTINUE
 2900       CONTINUE
C
         END IF
C
C        UPDATE M(T)
C
         DO 2980 I = 1, K
            DO 2960 J = I, K
               SUM = 0.0D0
               DO 2940 K2 = 1, K
                  DO 2920 J2 = 1, K
                     SUM = SUM + TEMPM(I,K2)*INVF(K2,J2)*TEMPM(J,J2)
 2920             CONTINUE
 2940          CONTINUE
               MT(I,J) = MT(I,J) - SUM
               MT(J,I) = MT(I,J)
 2960       CONTINUE
 2980    CONTINUE
C
C        UPDATE L(T)
C
         DO 3060 L = 1, R
            DO 3040 I = 1, K
               DO 3020 J = 1, K
                  SUM = 0.0D0
                  IF (PEARL .AND. (L.GE.Q+1)) MAT(I,(L-1)*K+J) = SUM
                  IF (PEARL .AND. (L.GE.Q+1)) GO TO 3020
                  DO 3000 K2 = 1, K
                     SUM = SUM + TEMPK(I,(L-1)*K+K2)*TEMPL(K2,J)
 3000             CONTINUE
                  MAT(I,(L-1)*K+J) = TEMP(I,(L-1)*K+J) - SUM
 3020          CONTINUE
 3040       CONTINUE
 3060    CONTINUE
C
C        UPDATE V(T),SSQ, AND DETP
C
 3080    DO 3100 I = 1, K
            V(I,T) = W3(I,T) - Z(I)
 3100    CONTINUE
C
         IF ((INTEG.EQ.2) .AND. ( .NOT. CONDS)) THEN
C
            IF (DELTA) THEN
               DO 3120 I = 1, K
                  REZ(I) = V(I,T)
 3120          CONTINUE
            ELSE
C
               DO 3140 I = 1, K
                  IMZ(I) = V(I,T)
 3140          CONTINUE
               IFAIL = 1
               CALL F04ARF(A,K+1,IMZ,K,IMZ,REZ,IFAIL)
C
C              PUT STANDARDIZED V(T) INTO REZ
C
               DO 3180 I = 1, K
                  SUM = 0.0D0
                  DO 3160 I2 = 1, K
                     SUM = SUM + PA(I,I2)*IMZ(I2)
 3160             CONTINUE
                  REZ(I) = SUM
 3180          CONTINUE
C
            END IF
         END IF
C
         DO 3220 I = 1, K
            DO 3200 J = 1, K
               SSQ = SSQ + V(I,T)*INVF(I,J)*V(J,T)
 3200       CONTINUE
 3220    CONTINUE
C
         IF (DELTA) GO TO 3240
C
C        UPDATE PRODUCT OF DETERMINANTS
C
         DETP = DETP + LOG(SIG/TSIG)
C
 3240 CONTINUE
C
      SUM = (SSQ+DBLE(N)*LOG(TSIG)+DETP)/2.0D0
      IF (INTEG.EQ.0) FBIG = SUM
C
      RETURN
      END
