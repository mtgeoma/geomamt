      SUBROUTINE S14ACZ(X,N,KODE,M,ANS,IER)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Based on the routine PSIFN by D.E. Amos, June, 1982.
C
C     References
C         Handbook of Mathematical Functions, AMS 55, National Bureau
C         of Standards by M. Abramowitz and I.A. Stegun, 1964, pp.
C         258-260, eqtns. 6.3.5, 6.3.18, 6.4.6, 6.4.9, 6.4.10
C
C         ACM Trans. Math Software, 1983.
C
C     Abstract
C         S14ACZ computes M member sequences of scaled derivatives of
C         the PSI function
C
C              W(K,X)=(-1)**(K+1)*PSI(K,X)/GAMMA(K+1)
C
C         K=N,...,N+M-1 where PSI(K,X) is the K-th derivative of the PSI
C         function.  On KODE=1, S14ACZ returns the scaled derivatives
C         as described.  KODE=2 is operative only when K=0 and S14ACZ
C         returns -PSI(X) + LN(X).  That is, the logarithmic behavior
C         for large X is removed when KODE=2 and K=0. When sums or
C         differences of PSI functions are computed, the logarithmic
C         terms can be combined analytically and computed separately
C         to help retain significant digits.
C
C         The basic method of evaluation is the asymptotic expansion
C         for large X.ge.XMIN followed by backward recursion on a two
C         term recursion relation
C
C                  W(X+1) + X**(-N-1) = W(X).
C
C         This is supplemented by a series
C
C                  SUM( (X+K)**(-N-1) , K=0,1,2,... )
C
C         which converges rapidly for large N. Both XMIN and the
C         number of terms of the series are calculated from the unit
C         round off of the machine environment.
C
C         The nominal computational accuracy is the maximum of unit
C         roundoff and 1.0E-18 since critical constants are given to
C         only 18 digits.
C
C     Description of arguments
C
C         Input
C           X      - argument, X.gt.0.0E0
C           N      - first member of the sequence, 0.le.N.le.100
C                    N=0 gives ANS(1) = -PSI(X)       on KODE=1
C                                       -PSI(X)+LN(X) on KODE=2
C           KODE   - selection parameter
C                    KODE=1 returns scaled derivatives of the PSI
C                    function
C                    KODE=2 returns scaled derivatives of the PSI
C                    function except when N=0. In this case,
C                    ANS(1) = -PSI(X) + LN(X) is returned.
C           M      - number of members of the sequence, M.ge.1
C
C         Output
C           ANS    - a vector of length at least M whose first M
C                    components contain the sequence of derivatives
C                    scaled according to KODE
C
C     .. Parameters ..
      INTEGER           NMAX
      PARAMETER         (NMAX=100)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X
      INTEGER           IER, KODE, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  ANS(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  ARG, DEN, ELIM, EPS, FLN, FN, FNP, FNS, FX,
     *                  R1M4, R1M5, RLN, RXSQ, S, SLOPE, T, T1, T2, TA,
     *                  TK, TOL, TOLS, TSS, TST, TT, WDTOL, XDMLN, XDMY,
     *                  XINC, XLN, XM, XMIN, XQ, YINT
      INTEGER           I, J, K, MX, NN, NP, NX
C     .. Local Arrays ..
      DOUBLE PRECISION  B(22), TRM(22), TRMR(NMAX)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      INTEGER           X02BHF, X02BJF, X02BKF, X02BLF
      EXTERNAL          X02AJF, X02AMF, X02BHF, X02BJF, X02BKF, X02BLF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, EXP, INT, LOG, LOG10, MAX, MIN, DBLE
C     .. Data statements ..
C     Bernoulli numbers:
      DATA              B(1), B(2), B(3), B(4), B(5), B(6), B(7), B(8),
     *                  B(9), B(10), B(11)/1.00000000000000000D+00,
     *                  -5.00000000000000000D-01,
     *                  1.66666666666666667D-01,
     *                  -3.33333333333333333D-02,
     *                  2.38095238095238095D-02,
     *                  -3.33333333333333333D-02,
     *                  7.57575757575757576D-02,
     *                  -2.53113553113553114D-01,
     *                  1.16666666666666667D+00,
     *                  -7.09215686274509804D+00,
     *                  5.49711779448621554D+01/
      DATA              B(12), B(13), B(14), B(15), B(16), B(17), B(18),
     *                  B(19), B(20), B(21), B(22)
     *                  /-5.29124242424242424D+02,
     *                  6.19212318840579710D+03,
     *                  -8.65802531135531136D+04,
     *                  1.42551716666666667D+06,
     *                  -2.72982310678160920D+07,
     *                  6.01580873900642368D+08,
     *                  -1.51163157670921569D+10,
     *                  4.29614643061166667D+11,
     *                  -1.37116552050883328D+13,
     *                  4.88332318973593167D+14,
     *                  -1.92965793419400681D+16/
C     .. Executable Statements ..
C
      IER = 0
      IF (X.LE.0.0D0) THEN
         IER = 1
      ELSE IF (N.LT.0) THEN
         IER = 2
      ELSE IF (KODE.NE.1 .AND. KODE.NE.2) THEN
         IER = 3
      ELSE IF (M.LT.1) THEN
         IER = 4
      ELSE
         NN = N + M - 1
         FN = NN
         FNP = FN + 1.0D0
         NX = MIN(-X02BKF(),X02BLF())
         R1M5 = LOG10(DBLE(X02BHF()))
         R1M4 = X02AJF()
         WDTOL = MAX(R1M4,0.5D-18)
C        ELIM = approximate exponential over and underflow limit.
C        ELIM = 2.302*(NX*R1M5-3.0)
         ELIM = -LOG(X02AMF())
         XLN = LOG(X)
         T = FNP*XLN
C        Overflow and underflow test for small and large X.
         IF (ABS(T).GT.ELIM) THEN
            IF (T.GT.0.0D+0) THEN
               IER = 5
            ELSE
               IER = 6
            END IF
         ELSE IF (X.LT.WDTOL) THEN
C           Small X.lt.unit round off.
            ANS(1) = X**(-N-1)
            DO 20 K = 2, M
               ANS(K) = ANS(K-1)/X
   20       CONTINUE
            IF (N.EQ.0 .AND. KODE.EQ.2) ANS(1) = ANS(1) + XLN
         ELSE
C           Compute XMIN and the number of terms of the series, FLN+1.
            RLN = R1M5*X02BJF()
            RLN = MIN(RLN,18.06D0)
            FLN = MAX(RLN,3.0D0) - 3.0D0
            YINT = 3.50D0 + 0.40D0*FLN
            SLOPE = 0.21D0 + FLN*(0.0006038D0*FLN+0.008677D0)
            XM = YINT + SLOPE*FN
            MX = INT(XM) + 1
            XMIN = DBLE(MX)
            IF (N.NE.0) THEN
               XM = -2.302D0*RLN - MIN(0.0D0,XLN)
               FNS = N
               ARG = XM/FNS
               ARG = MIN(0.0D0,ARG)
               EPS = EXP(ARG)
               XM = 1.0D0 - EPS
               IF (ABS(ARG).LT.1.0D-3) XM = -ARG
               FLN = X*XM/EPS
               XM = XMIN - X
               IF (XM.GT.7.0D0 .AND. FLN.LT.15.0D0) THEN
C                 Compute by series (X+K)**(-(N+1)) , K=0,1,2,...
                  NN = INT(FLN) + 1
                  NP = N + 1
                  T1 = (FNS+1.0D0)*XLN
                  T = EXP(-T1)
                  S = T
                  DEN = X
                  DO 40 I = 1, NN
                     DEN = DEN + 1.0D0
                     TRM(I) = DEN**(-NP)
                     S = S + TRM(I)
   40             CONTINUE
                  ANS(1) = S
                  IF (N.EQ.0) THEN
                     IF (KODE.EQ.2) ANS(1) = S + XLN
                  END IF
                  IF (M.GT.1) THEN
C                    Generate higher derivatives, J.gt.N .
                     TOL = WDTOL/5.0D0
                     DO 100 J = 2, M
                        T = T/X
                        S = T
                        TOLS = T*TOL
                        DEN = X
                        DO 60 I = 1, NN
                           DEN = DEN + 1.0D0
                           TRM(I) = TRM(I)/DEN
                           S = S + TRM(I)
                           IF (TRM(I).LT.TOLS) GO TO 80
   60                   CONTINUE
   80                   ANS(J) = S
  100                CONTINUE
                  END IF
                  GO TO 320
               END IF
            END IF
            XDMY = X
            XDMLN = XLN
            XINC = 0.0D0
            IF (X.LT.XMIN) THEN
               NX = INT(X)
               XINC = XMIN - NX
               XDMY = X + XINC
               XDMLN = LOG(XDMY)
            END IF
C           Generate W(N+M-1,X) by the asymptotic expansion.
            T = FN*XDMLN
            T1 = XDMLN + XDMLN
            T2 = T + XDMLN
            TK = MAX(ABS(T),ABS(T1),ABS(T2))
            IF (TK.GT.ELIM) THEN
               IER = 5
            ELSE
               TSS = EXP(-T)
               TT = 0.5D0/XDMY
               T1 = TT
               TST = WDTOL*TT
               IF (NN.NE.0) T1 = TT + 1.0D0/FN
               RXSQ = 1.0D0/(XDMY*XDMY)
               TA = 0.5D0*RXSQ
               T = FNP*TA
               S = T*B(3)
               IF (ABS(S).GE.TST) THEN
                  TK = 2.0D0
                  DO 120 K = 4, 22
                     T = T*((TK+FN+1.0D0)/(TK+1.0D0))*((TK+FN)
     *                   /(TK+2.0D0))*RXSQ
                     TRM(K) = T*B(K)
                     IF (ABS(TRM(K)).LT.TST) THEN
                        GO TO 140
                     ELSE
                        S = S + TRM(K)
                        TK = TK + 2.0D0
                     END IF
  120             CONTINUE
               END IF
  140          S = (S+T1)*TSS
               IF (XINC.NE.0.0D0) THEN
C                 Backward recur from XDMY to X.
                  NX = INT(XINC)
                  NP = NN + 1
                  IF (NX.GT.NMAX) THEN
                     IER = 7
                     GO TO 320
                  ELSE IF (NN.EQ.0) THEN
                     GO TO 260
                  ELSE
                     XM = XINC - 1.0D0
                     FX = X + XM
C                    This loop should not be changed. FX is accurate
C                    when X is small.
                     DO 160 I = 1, NX
                        TRMR(I) = FX**(-NP)
                        S = S + TRMR(I)
                        XM = XM - 1.0D0
                        FX = X + XM
  160                CONTINUE
                  END IF
               END IF
               ANS(M) = S
               IF (FN.EQ.0.0D0) THEN
                  GO TO 300
               ELSE IF (M.EQ.1) THEN
                  GO TO 320
               ELSE
C                 Generate lower derivatives, J.lt.N+M-1 .
                  DO 240 J = 2, M
                     FNP = FN
                     FN = FN - 1.0D0
                     TSS = TSS*XDMY
                     T1 = TT
                     IF (FN.NE.0.0D0) T1 = TT + 1.0D0/FN
                     T = FNP*TA
                     S = T*B(3)
                     IF (ABS(S).GE.TST) THEN
                        TK = 3.0D0 + FNP
                        DO 180 K = 4, 22
                           TRM(K) = TRM(K)*FNP/TK
                           IF (ABS(TRM(K)).LT.TST) THEN
                              GO TO 200
                           ELSE
                              S = S + TRM(K)
                              TK = TK + 2.0D0
                           END IF
  180                   CONTINUE
                     END IF
  200                S = (S+T1)*TSS
                     IF (XINC.NE.0.0D0) THEN
                        IF (FN.EQ.0.0D0) THEN
                           GO TO 260
                        ELSE
                           XM = XINC - 1.0D0
                           FX = X + XM
                           DO 220 I = 1, NX
                              TRMR(I) = TRMR(I)*FX
                              S = S + TRMR(I)
                              XM = XM - 1.0D0
                              FX = X + XM
  220                      CONTINUE
                        END IF
                     END IF
                     MX = M - J + 1
                     ANS(MX) = S
                     IF (FN.EQ.0.0D0) GO TO 300
  240             CONTINUE
                  GO TO 320
               END IF
C              Recursion for N = 0.
  260          CONTINUE
               DO 280 I = 1, NX
                  S = S + 1.0D0/(X+(NX-I))
  280          CONTINUE
  300          IF (KODE.NE.2) THEN
                  ANS(1) = S - XDMLN
               ELSE IF (XDMY.NE.X) THEN
                  XQ = XDMY/X
                  ANS(1) = S - LOG(XQ)
               END IF
            END IF
         END IF
      END IF
  320 RETURN
C
      END
