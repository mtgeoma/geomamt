      SUBROUTINE C06ECQ(X,Y,PTS,M,P)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     RADIX PRIME COMPLEX FOURIER TRANSFORM KERNEL
C     .. Scalar Arguments ..
      INTEGER           M, P, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X(PTS), Y(PTS)
C     .. Local Scalars ..
      DOUBLE PRECISION  ANGLE, IS, IU, RS, RU, T, TWOPI, XT, YT
      INTEGER           J, JJ, K, K0, KS1, KS2, MOVER2, MP, PM, PP, U, V
      LOGICAL           FOLD, ZERO
C     .. Local Arrays ..
      DOUBLE PRECISION  A(18), AA(9,9), B(18), BB(9,9), C(18), IA(9),
     *                  IB(9), RA(9), RB(9), S(18)
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         COS, DBLE, SIN
C     .. Executable Statements ..
      TWOPI = 2.0D0*X01AAF(0.0D0)
      MOVER2 = M/2 + 1
      MP = M*P
      PP = P/2
      PM = P - 1
      DO 20 U = 1, PP
         JJ = P - U
         ANGLE = TWOPI*DBLE(U)/DBLE(P)
         A(U) = COS(ANGLE)
         B(U) = SIN(ANGLE)
         A(JJ) = A(U)
         B(JJ) = -B(U)
   20 CONTINUE
      DO 60 U = 1, PP
         DO 40 V = 1, PP
            JJ = U*V - ((U*V)/P)*P
            AA(V,U) = A(JJ)
            BB(V,U) = B(JJ)
   40    CONTINUE
   60 CONTINUE
C
      DO 300 J = 1, MOVER2
         FOLD = J .GT. 1 .AND. 2*J .LT. M + 2
         K0 = J
         ANGLE = TWOPI*DBLE(J-1)/DBLE(MP)
         ZERO = ANGLE .EQ. 0.0D0
         C(1) = COS(ANGLE)
         S(1) = SIN(ANGLE)
         DO 80 U = 2, PM
            C(U) = C(U-1)*C(1) - S(U-1)*S(1)
            S(U) = S(U-1)*C(1) + C(U-1)*S(1)
   80    CONTINUE
         GO TO 140
  100    CONTINUE
         FOLD = .FALSE.
         K0 = M + 2 - J
         DO 120 U = 1, PM
            T = C(U)*A(U) + S(U)*B(U)
            S(U) = -S(U)*A(U) + C(U)*B(U)
            C(U) = T
  120    CONTINUE
  140    CONTINUE
C
         DO 280 K = K0, PTS, MP
            XT = X(K)
            YT = Y(K)
            KS1 = M + K
            KS2 = (P-1)*M + K
            RS = X(KS1) + X(KS2)
            IS = Y(KS1) + Y(KS2)
            RU = X(KS1) - X(KS2)
            IU = Y(KS1) - Y(KS2)
            DO 160 U = 1, PP
               RA(U) = XT + RS*AA(U,1)
               IA(U) = YT + IS*AA(U,1)
               RB(U) = RU*BB(U,1)
               IB(U) = IU*BB(U,1)
  160       CONTINUE
            XT = XT + RS
            YT = YT + IS
            DO 200 U = 2, PP
               JJ = P - U
               KS1 = U*M + K
               KS2 = JJ*M + K
               RS = X(KS1) + X(KS2)
               IS = Y(KS1) + Y(KS2)
               RU = X(KS1) - X(KS2)
               IU = Y(KS1) - Y(KS2)
               XT = XT + RS
               YT = YT + IS
               DO 180 V = 1, PP
                  RA(V) = RA(V) + RS*AA(V,U)
                  IA(V) = IA(V) + IS*AA(V,U)
                  RB(V) = RB(V) + RU*BB(V,U)
                  IB(V) = IB(V) + IU*BB(V,U)
  180          CONTINUE
  200       CONTINUE
            X(K) = XT
            Y(K) = YT
            DO 260 U = 1, PP
               JJ = P - U
               IF (ZERO) GO TO 220
               XT = RA(U) + IB(U)
               YT = IA(U) - RB(U)
               KS1 = U*M + K
               X(KS1) = XT*C(U) + YT*S(U)
               Y(KS1) = YT*C(U) - XT*S(U)
               XT = RA(U) - IB(U)
               YT = IA(U) + RB(U)
               KS1 = JJ*M + K
               X(KS1) = XT*C(JJ) + YT*S(JJ)
               Y(KS1) = YT*C(JJ) - XT*S(JJ)
               GO TO 240
  220          CONTINUE
               KS1 = U*M + K
               X(KS1) = RA(U) + IB(U)
               Y(KS1) = IA(U) - RB(U)
               KS1 = JJ*M + K
               X(KS1) = RA(U) - IB(U)
               Y(KS1) = IA(U) + RB(U)
  240          CONTINUE
  260       CONTINUE
  280    CONTINUE
         IF (FOLD) GO TO 100
  300 CONTINUE
C
      RETURN
      END
