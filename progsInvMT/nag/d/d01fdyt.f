      SUBROUTINE D01FDY(MAP,IR,NDIM,F,REGION)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     REGION
C
C     REPROCESS MAP(*)
C     NUM(I) WILL CONTAIN THE NUMBER OF CO-ORDINATE PARAMETERS
C     EQUAL TO IY(I),I=1,..,M
C
C     .. Scalar Arguments ..
      INTEGER           IR, NDIM
C     .. Array Arguments ..
      INTEGER           MAP(30)
C     .. Function Arguments ..
      DOUBLE PRECISION  F
      EXTERNAL          F
C     .. Subroutine Arguments ..
      EXTERNAL          REGION
C     .. Scalars in Common ..
      DOUBLE PRECISION  FACT, QF, RMOD, TMAX, TOT
      INTEGER           ICOUNT, ISG
C     .. Arrays in Common ..
      DOUBLE PRECISION  YY(30)
C     .. Local Scalars ..
      DOUBLE PRECISION  Q, TS, VAL, YL
      INTEGER           I, IF, ILIM, IQ, IYY, J, J1, J2, JP, K, KK, L1,
     *                  LC, M, MI
C     .. Local Arrays ..
      DOUBLE PRECISION  B(30), CHI(30)
      INTEGER           IY(30), NUM(30)
C     .. External Functions ..
      DOUBLE PRECISION  S10AAF
      EXTERNAL          S10AAF
C     .. External Subroutines ..
      EXTERNAL          D01FDX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AD01FD/YY, RMOD, TMAX, TOT, QF, FACT, ISG,
     *                  ICOUNT
C     .. Executable Statements ..
      M = 0
      IF (NDIM.EQ.IR) GO TO 20
      IY(1) = 1
      NUM(1) = NDIM - IR
      M = 1
   20 L1 = MAP(1)
      LC = 0
      DO 60 J = 1, IR
         IF (L1.NE.MAP(J)) GO TO 40
         LC = LC + 1
         GO TO 60
   40    M = M + 1
         NUM(M) = LC
         IY(M) = L1
         L1 = MAP(J)
         LC = 1
   60 CONTINUE
      M = M + 1
      NUM(M) = LC
      IY(M) = L1
      IF (M.EQ.1) GO TO 120
C
C     SORT DESCENDING
C
      ILIM = M - 1
      DO 100 I = 1, ILIM
         IQ = 0
         MI = M - I
         DO 80 J = 1, MI
            JP = J + 1
            J1 = NUM(J)
            J2 = NUM(JP)
            IF (J1.GE.J2) GO TO 80
            IQ = 1
            NUM(J) = J2
            NUM(JP) = J1
            J1 = IY(J)
            IY(J) = IY(JP)
            IY(JP) = J1
   80    CONTINUE
         IF (IQ.EQ.0) GO TO 120
  100 CONTINUE
  120 IF (ISG.EQ.0) GO TO 160
C
C     CALCULATE EVALUATION POINTS FOR SPHERE TO SPHERE
C
      DO 140 I = 1, M
         IYY = (IY(I)+1)/2
         CHI(I) = QF*YY(IYY)
  140 CONTINUE
      GO TO 200
C
C     COMPUTE LAYER JACOBIAN AND EVALUATION POINTS FOR PRODUCT
C     REGION
C
  160 FACT = RMOD**(NDIM+1)
      DO 180 I = 1, M
         IYY = (IY(I)+1)/2
         Q = YY(IYY)
         TS = RMOD*Q
         IF (ABS(TS).GE.TMAX) RETURN
         IF = 0
         Q = S10AAF(TS,IF)
         CHI(I) = Q
         FACT = FACT*(1.0D0-Q*Q)**NUM(I)
  180 CONTINUE
C
C     PERMUTATIONS
C
  200 IF (M.NE.1) GO TO 220
      YL = CHI(1)
      K = 0
      GO TO 240
  220 YL = CHI(M-1)
      K = NUM(M-1)
  240 KK = NUM(M)
      DO 260 I = 1, KK
         B(I) = CHI(M)
  260 CONTINUE
C
C     CALCULATE CONTRIBUTION FROM ALL PERMUTATION OF POINTS
C
      CALL D01FDX(YL,K,M,B,NUM,CHI,VAL,F,REGION)
      TOT = TOT + VAL*FACT
      RETURN
      END
