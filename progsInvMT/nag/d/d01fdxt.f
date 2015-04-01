      SUBROUTINE D01FDX(YLS,MS,JLS,B,NUM,X,VAL,F,REGION)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     CALCULATE ALL PERMUTATIONS WITH NUM(I) ALIKE VALUES OF
C     X(I),I=1,..,MS
C     BASED ON CACM ALGORITHM 242
C     REGION
C     .. Scalar Arguments ..
      DOUBLE PRECISION  VAL, YLS
      INTEGER           JLS, MS
C     .. Array Arguments ..
      DOUBLE PRECISION  B(30), X(30)
      INTEGER           NUM(30)
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
      DOUBLE PRECISION  A, XL
      INTEGER           I, IP, IP2, ISUB, J1, JL, KK, L, M, N, N1, N2,
     *                  NM, NS
C     .. Local Arrays ..
      DOUBLE PRECISION  STACK(400), Y(30), YM(30)
      INTEGER           ISTACK(200), J(30)
C     .. External Functions ..
      DOUBLE PRECISION  D01FDW
      EXTERNAL          D01FDW
C     .. Common blocks ..
      COMMON            /AD01FD/YY, RMOD, TMAX, TOT, QF, FACT, ISG,
     *                  ICOUNT
C     .. Executable Statements ..
      NS = NUM(JLS)
      VAL = 0.0D0
      L = 0
      IP = 1
      IP2 = 1
C
C     INITIAL PSEUDO ARGUMENTS
C
      XL = YLS
      M = MS
      JL = JLS
      N = NS
C
C     PSEUDO ENTRY
C
   20 L = L + 1
      N2 = N + M
      IF (M.EQ.0) GO TO 60
      NM = N + 1
      DO 40 I = NM, N2
         Y(I) = XL
   40 CONTINUE
   60 DO 80 I = 1, N
         J(I) = I
   80 CONTINUE
      J(N+1) = N2 + 1
      J1 = JL - 1
      KK = N
  100 DO 120 I = 1, KK
         ISUB = J(I)
         Y(ISUB) = B(I)
  120 CONTINUE
      IF (J1.GT.1) GO TO 220
      IF (ISG.EQ.0) GO TO 160
C
C     SPHERE
C
      VAL = VAL + F(N2,Y)
      DO 140 I = 1, N2
         YM(I) = -Y(I)
  140 CONTINUE
      VAL = VAL + F(N2,YM)
      GO TO 200
C
C     PRODUCT REGION
C
  160 VAL = VAL + D01FDW(F,Y,N2,REGION)
      DO 180 I = 1, N2
         YM(I) = -Y(I)
  180 CONTINUE
      VAL = VAL + D01FDW(F,YM,N2,REGION)
  200 ICOUNT = ICOUNT + 2
      GO TO 400
  220 A = X(J1-1)
      N1 = NUM(J1-1)
      IP = IP - 1
      NM = N + 1
      DO 240 I = 1, NM
         ISUB = IP + I
         ISTACK(ISUB) = J(I)
  240 CONTINUE
      IP = IP + N + 2
      ISTACK(IP) = N
      IP = IP + 1
      ISTACK(IP) = N2
      IP = IP + 1
      ISTACK(IP) = J1
      IP = IP + 1
      IP2 = IP2 - 1
      DO 260 I = 1, N2
         ISUB = IP2 + I
         STACK(ISUB) = Y(I)
  260 CONTINUE
      IP2 = IP2 + N2
      DO 280 I = 1, N
         ISUB = IP2 + I
         STACK(ISUB) = B(I)
  280 CONTINUE
      IP2 = IP2 + N + 1
      STACK(IP2) = XL
      IP2 = IP2 + 1
      XL = A
      M = N1
      N = N2
      JL = J1
      DO 300 I = 1, N
         B(I) = Y(I)
  300 CONTINUE
      GO TO 20
C
C     PSEUDO RETURN
C
  320 L = L - 1
      IF (L.EQ.0) RETURN
      IP = IP - 1
      J1 = ISTACK(IP)
      IP = IP - 1
      N2 = ISTACK(IP)
      IP = IP - 1
      N = ISTACK(IP)
      IP = IP - N - 2
      NM = N + 1
      DO 340 I = 1, NM
         ISUB = IP + I
         J(I) = ISTACK(ISUB)
  340 CONTINUE
      IP = IP + 1
      IP2 = IP2 - 1
      XL = STACK(IP2)
      IP2 = IP2 - N - 1
      DO 360 I = 1, N
         ISUB = IP2 + I
         B(I) = STACK(ISUB)
  360 CONTINUE
      IP2 = IP2 - N2
      DO 380 I = 1, N2
         ISUB = IP2 + I
         Y(I) = STACK(ISUB)
  380 CONTINUE
      IP2 = IP2 + 1
  400 DO 440 I = 1, N
         ISUB = J(I)
         Y(ISUB) = XL
         J(I) = J(I) + 1
         IF (J(I).GE.J(I+1)) GO TO 420
         KK = I
         GO TO 100
  420    J(I) = I
  440 CONTINUE
      GO TO 320
      END
