      SUBROUTINE E02GAF(M,A,LA,B,NPLUS2,TOL,X,RESID,IRANK,ITER,IWORK,
     *                  IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 10A REVISED. IER-396 (OCT 1982).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C
C     THIS SUBROUTINE USES A MODIFICATION OF THE SIMPLEX METHOD
C     OF LINEAR PROGRAMMING TO CALCULATE AN L1 SOLUTION TO AN
C     OVER-DETERMINED SYSTEM OF LINEAR EQUATIONS.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E02GAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RESID, TOL
      INTEGER           IFAIL, IRANK, ITER, LA, M, NPLUS2
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LA,NPLUS2), B(M), X(NPLUS2)
      INTEGER           IWORK(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  DZERO, SUM
      DOUBLE PRECISION  BIG, CSING, D, MAX, MIN, ONE, ONENEG, PIVOT,
     *                  RMLT, SSING, TOLER, TPIVOT, TSING, TWO, TWONEG,
     *                  YSING, ZERO
      INTEGER           I, IN, IONE, ISAVE, IZERO, J, K, KL, KOUNT, KR,
     *                  L, M1, M2, N, N1, N2, OUT
      LOGICAL           STAGE, TEST
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02ALF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, X02ALF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          E02GAZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE
C     .. Data statements ..
      DATA              IZERO, IONE, ZERO, ONE, TWO, DZERO, ONENEG,
     *                  TWONEG/0, 1, 0.0D0, 1.0D0, 2.0D0, 0.D0, -1.0D0,
     *                  -2.0D0/
C     .. Executable Statements ..
      TOLER = TOL
      IF (TOLER.LE.0.0D0) TOLER = X02AJF()**(2.0D0/3.0D0)
      BIG = X02ALF()
      ISAVE = IFAIL
      N = NPLUS2 - 2
C     TEST INPUT VALUES
      IF ((M.GE.N) .AND. (LA.GE.M+2) .AND. (N.GE.1)) GO TO 20
      SUM = BIG
      IFAIL = 3
      GO TO 800
   20 M1 = M + IONE
      N1 = N + IONE
      M2 = M1 + IONE
      N2 = N1 + IONE
      DO 40 J = 1, N
         A(M2,J) = J
         X(J) = ZERO
   40 CONTINUE
      DO 80 I = 1, M
         A(I,N2) = N + I
         A(I,N1) = B(I)
         IF (B(I).GE.ZERO) GO TO 80
         DO 60 J = 1, N2
            A(I,J) = -A(I,J)
   60    CONTINUE
   80 CONTINUE
C     COMPUTE THE MARGINAL COSTS IN HIGHER PRECISION.
      DO 120 J = 1, N1
         SSING = 0.0D0
         CSING = 0.0D0
         DO 100 I = 1, M
            YSING = CSING + A(I,J)
            TSING = SSING + YSING
            CSING = (SSING-TSING) + YSING
            SSING = TSING
  100    CONTINUE
         A(M1,J) = SSING + CSING
  120 CONTINUE
C     STAGE I.
C     DETERMINE THE VECTOR TO ENTER THE BASIS.
      STAGE = .TRUE.
      KOUNT = IZERO
      KR = IONE
      KL = IONE
  140 MAX = ONENEG
      DO 160 J = KR, N
         IF (ABS(A(M2,J)).GT.DBLE(N)) GO TO 160
         D = ABS(A(M1,J))
         IF (D.LE.MAX) GO TO 160
         MAX = D
         IN = J
  160 CONTINUE
      IF (A(M1,IN).GE.ZERO) GO TO 200
      DO 180 I = 1, M2
         A(I,IN) = -A(I,IN)
  180 CONTINUE
C     DETERMINE THE VECTOR TO LEAVE THE BASIS.
  200 K = IZERO
      DO 220 I = KL, M
         D = A(I,IN)
         IF (D.LE.TOLER) GO TO 220
         K = K + IONE
         B(K) = A(I,N1)/D
         IWORK(K) = I
         TEST = .TRUE.
  220 CONTINUE
  240 IF (K.GT.IZERO) GO TO 260
      TEST = .FALSE.
      GO TO 300
  260 MIN = BIG
      DO 280 I = 1, K
         IF (B(I).GE.MIN) GO TO 280
         J = I
         MIN = B(I)
         OUT = IWORK(I)
  280 CONTINUE
      B(J) = B(K)
      IWORK(J) = IWORK(K)
      K = K - IONE
C     CHECK FOR LINEAR DEPENDENCE IN STAGE I.
  300 IF (TEST .OR. .NOT. STAGE) GO TO 340
      DO 320 I = 1, M2
         D = A(I,KR)
         A(I,KR) = A(I,IN)
         A(I,IN) = D
  320 CONTINUE
      KR = KR + IONE
      GO TO 500
  340 IF (TEST) GO TO 360
      IFAIL = 2
      GO TO 680
  360 PIVOT = A(OUT,IN)
      IF (A(M1,IN)-PIVOT-PIVOT.LE.TOLER) GO TO 400
      DO 380 J = KR, N1
         D = A(OUT,J)
         A(M1,J) = A(M1,J) - D - D
         A(OUT,J) = -D
  380 CONTINUE
      A(OUT,N2) = -A(OUT,N2)
      GO TO 240
C     PIVOT ON A(OUT,IN).
  400 DO 420 J = KR, N1
         IF (J.EQ.IN) GO TO 420
         A(OUT,J) = A(OUT,J)/PIVOT
  420 CONTINUE
      DO 440 J = KR, N1
         IF (J.EQ.IN) GO TO 440
         RMLT = -A(OUT,J)
         CALL E02GAZ(A(1,J),A(1,IN),RMLT,M1,OUT)
  440 CONTINUE
      TPIVOT = -PIVOT
      DO 460 I = 1, M1
         IF (I.EQ.OUT) GO TO 460
         A(I,IN) = A(I,IN)/TPIVOT
  460 CONTINUE
      A(OUT,IN) = ONE/PIVOT
      D = A(OUT,N2)
      A(OUT,N2) = A(M2,IN)
      A(M2,IN) = D
      KOUNT = KOUNT + IONE
      IF ( .NOT. STAGE) GO TO 520
C     INTERCHANGE ROWS IN STAGE I.
      KL = KL + IONE
      DO 480 J = KR, N2
         D = A(OUT,J)
         A(OUT,J) = A(KOUNT,J)
         A(KOUNT,J) = D
  480 CONTINUE
  500 IF (KOUNT+KR.NE.N1) GO TO 140
C     STAGE II.
      STAGE = .FALSE.
C     DETERMINE THE VECTOR TO ENTER THE BASIS.
  520 MAX = -BIG
      DO 560 J = KR, N
         D = A(M1,J)
         IF (D.GE.ZERO) GO TO 540
         IF (D.GT.TWONEG) GO TO 560
         D = TWONEG - D
  540    IF (D.LE.MAX) GO TO 560
         MAX = D
         IN = J
  560 CONTINUE
      IF (MAX.LE.TOLER) GO TO 600
      IF (A(M1,IN).GT.ZERO) GO TO 200
      DO 580 I = 1, M2
         A(I,IN) = -A(I,IN)
  580 CONTINUE
      A(M1,IN) = A(M1,IN) + TWONEG
      GO TO 200
C     PREPARE OUTPUT.
  600 L = KL - IONE
      DO 640 I = 1, L
         IF (A(I,N1).GE.ZERO) GO TO 640
         DO 620 J = KR, N2
            A(I,J) = -A(I,J)
  620    CONTINUE
  640 CONTINUE
      IFAIL = 1
      IF (KR.NE.IONE) GO TO 680
      DO 660 J = 1, N
         D = ABS(A(M1,J))
         IF (D.LE.TOLER .OR. TWO-D.LE.TOLER) GO TO 680
  660 CONTINUE
      IFAIL = 0
  680 DO 700 I = 1, M
         B(I) = ZERO
  700 CONTINUE
      DO 760 I = 1, M
         K = A(I,N2)
         D = A(I,N1)
         IF (K.GT.IZERO) GO TO 720
         K = -K
         D = -D
  720    IF (I.GE.KL) GO TO 740
         X(K) = D
         GO TO 760
  740    K = K - N
         B(K) = D
  760 CONTINUE
      ITER = KOUNT
      IRANK = N1 - KR
      SUM = DZERO
      IF (KL.GT.M) GO TO 800
      DO 780 I = KL, M
         SUM = SUM + A(I,N1)
  780 CONTINUE
  800 RESID = SUM
      IF (IFAIL.NE.0) IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
      END
