      SUBROUTINE C06EBQ(X,PTS,M,P)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 9 REVISED. IER-303 (SEP 1981).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX PRIME HERMITE FOURIER TRANSFORM KERNEL
C     .. Scalar Arguments ..
      INTEGER           M, P, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X(PTS)
C     .. Local Scalars ..
      DOUBLE PRECISION  ANGLE, IS, IU, RS, RU, T, TWOPI, XT, YT
      INTEGER           J, J0, J1, JJ, K, K0, K1, KS1, KS2, MOVER2, MP,
     *                  P2, PM, U, V
C     .. Local Arrays ..
      DOUBLE PRECISION  AA(9,9), BB(9,9), C(18), IA(9), IB(9), RA(9),
     *                  RB(9), S(18)
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         COS, DBLE, SIN
C     .. Executable Statements ..
      TWOPI = 2.0D0*X01AAF(0.0D0)
      MOVER2 = (M-1)/2
      MP = M*P
      P2 = P/2
      PM = P - 1
      DO 20 U = 1, P2
         JJ = P - U
         ANGLE = TWOPI*DBLE(U)/DBLE(P)
         C(U) = COS(ANGLE)
         S(U) = SIN(ANGLE)
         C(JJ) = C(U)
         S(JJ) = -S(U)
   20 CONTINUE
      DO 60 U = 1, P2
         DO 40 V = 1, P2
            JJ = U*V - ((U*V)/P)*P
            AA(V,U) = C(JJ)
            BB(V,U) = S(JJ)
   40    CONTINUE
   60 CONTINUE
C
      DO 160 K = 1, PTS, MP
         XT = X(K)
         KS1 = M + K
         RS = 2.0D0*X(KS1)
         KS1 = (P-1)*M + K
         IU = 2.0D0*X(KS1)
         DO 80 U = 1, P2
            RA(U) = XT + RS*AA(U,1)
            IB(U) = IU*BB(U,1)
   80    CONTINUE
         XT = XT + RS
         DO 120 U = 2, P2
            JJ = P - U
            KS1 = U*M + K
            RS = 2.0D0*X(KS1)
            KS1 = JJ*M + K
            IU = 2.0D0*X(KS1)
            XT = XT + RS
            DO 100 V = 1, P2
               RA(V) = RA(V) + RS*AA(V,U)
               IB(V) = IB(V) + IU*BB(V,U)
  100       CONTINUE
  120    CONTINUE
         X(K) = XT
         DO 140 U = 1, P2
            JJ = P - U
            KS1 = U*M + K
            X(KS1) = RA(U) + IB(U)
            KS1 = JJ*M + K
            X(KS1) = RA(U) - IB(U)
  140    CONTINUE
  160 CONTINUE
C
      IF (MOVER2.LT.1) GO TO 320
      DO 300 J = 1, MOVER2
         J0 = J + 1
         J1 = M - 2*J
         ANGLE = TWOPI*DBLE(J)/DBLE(MP)
         C(1) = COS(ANGLE)
         S(1) = SIN(ANGLE)
         DO 180 U = 2, PM
            C(U) = C(U-1)*C(1) - S(U-1)*S(1)
            S(U) = S(U-1)*C(1) + C(U-1)*S(1)
  180    CONTINUE
C
         DO 280 K0 = J0, PTS, MP
            K1 = K0 + J1
            XT = X(K0)
            KS1 = (P-1)*M + K1
            YT = X(KS1)
            KS1 = M + K0
            RS = X(KS1) + X(K1)
            RU = X(KS1) - X(K1)
            KS1 = (P-2)*M + K1
            KS2 = (P-1)*M + K0
            IS = X(KS1) - X(KS2)
            IU = X(KS1) + X(KS2)
            DO 200 U = 1, P2
               RA(U) = XT + RS*AA(U,1)
               IA(U) = YT + IS*AA(U,1)
               RB(U) = RU*BB(U,1)
               IB(U) = IU*BB(U,1)
  200       CONTINUE
            XT = XT + RS
            YT = YT + IS
            DO 240 U = 2, P2
               JJ = P - U
               KS1 = U*M + K0
               KS2 = (U-1)*M + K1
               RS = X(KS1) + X(KS2)
               RU = X(KS1) - X(KS2)
               KS1 = (JJ-1)*M + K1
               KS2 = JJ*M + K0
               IS = X(KS1) - X(KS2)
               IU = X(KS1) + X(KS2)
               XT = XT + RS
               YT = YT + IS
               DO 220 V = 1, P2
                  RA(V) = RA(V) + RS*AA(V,U)
                  IA(V) = IA(V) + IS*AA(V,U)
                  RB(V) = RB(V) + RU*BB(V,U)
                  IB(V) = IB(V) + IU*BB(V,U)
  220          CONTINUE
  240       CONTINUE
            X(K0) = XT
            X(K1) = YT
            DO 260 U = 1, P2
               JJ = P - U
               XT = RA(U) + IB(U)
               YT = IA(U) - RB(U)
               KS1 = U*M + K0
               X(KS1) = XT*C(U) + YT*S(U)
               KS1 = U*M + K1
               X(KS1) = YT*C(U) - XT*S(U)
               XT = RA(U) - IB(U)
               YT = IA(U) + RB(U)
               KS1 = JJ*M + K0
               X(KS1) = XT*C(JJ) + YT*S(JJ)
               KS1 = JJ*M + K1
               X(KS1) = YT*C(JJ) - XT*S(JJ)
  260       CONTINUE
  280    CONTINUE
  300 CONTINUE
C
  320 CONTINUE
      IF (M.NE.(M/2)*2) GO TO 500
      J = M/2 + 1
      DO 360 V = 1, P2
         ANGLE = TWOPI*DBLE(V)/DBLE(2*P)
         C(V) = COS(ANGLE)
         S(V) = SIN(ANGLE)
         DO 340 U = 1, P2
            T = AA(V,U)*C(V) + BB(V,U)*S(V)
            BB(V,U) = BB(V,U)*C(V) - AA(V,U)*S(V)
            AA(V,U) = T
  340    CONTINUE
  360 CONTINUE
      DO 480 K = J, PTS, MP
         KS1 = P2*M + K
         XT = X(KS1)
         RS = 2.0D0*X(K)
         KS1 = (P-1)*M + K
         IU = -2.0D0*X(KS1)
         DO 380 U = 1, P2, 2
            RA(U) = -XT + RS*AA(U,1)
            IB(U) = IU*BB(U,1)
  380    CONTINUE
         DO 400 U = 2, P2, 2
            RA(U) = XT + RS*AA(U,1)
            IB(U) = IU*BB(U,1)
  400    CONTINUE
         XT = XT + RS
         DO 440 U = 2, P2
            JJ = P - U
            KS1 = (U-1)*M + K
            RS = 2.0D0*X(KS1)
            KS1 = JJ*M + K
            IU = -2.0D0*X(KS1)
            XT = XT + RS
            DO 420 V = 1, P2
               RA(V) = RA(V) + RS*AA(V,U)
               IB(V) = IB(V) + IU*BB(V,U)
  420       CONTINUE
  440    CONTINUE
         X(K) = XT
         DO 460 U = 1, P2
            JJ = P - U
            KS1 = U*M + K
            X(KS1) = RA(U) - IB(U)
            KS1 = JJ*M + K
            X(KS1) = -RA(U) - IB(U)
  460    CONTINUE
  480 CONTINUE
C
  500 CONTINUE
      RETURN
      END
