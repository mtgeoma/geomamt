      SUBROUTINE E02GCF(M,N,MDIM,NDIM,A,B,TOL1,RELER,X,RESMAX,IRANK,
     *                  ITER,IFAIL)
C     MARK 8 RELEASE.  NAG COPYRIGHT 1980.
C     MARK 9 REVISED. IER-316 (SEP 1981).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C
C     THIS SUBROUTINE USES A MODIFICATION OF THE SIMPLEX METHOD
C     OF LINEAR PROGRAMMING TO CALCULATE A CHEBYSHEV SOLUTION TO
C     AN OVER-DETERMINED SYSTEM OF LINEAR EQUATIONS.
C
C     DERIVED FROM ACM ALGORITHM 495 BY I. BARRODALE AND C. PHILLIPS
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E02GCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RELER, RESMAX, TOL1
      INTEGER           IFAIL, IRANK, ITER, M, MDIM, N, NDIM
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NDIM,MDIM), B(M), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIG, D, DD, ONE, PIVOT, RELERR, RELTMP, TOL,
     *                  TPIVOT, TWO, VAL, ZERO
      INTEGER           I, IER, IRNKP1, J, K, LEV, MM1, MODE, MP1, NP1,
     *                  NP1MK, NP1MR, NP2, NP3, PCOL, PROW
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02ALF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, X02ALF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          E02GCZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SIGN, DBLE
C     .. Data statements ..
      DATA              ONE, TWO, ZERO/1.D0, 2.D0, 0.D0/
C     .. Executable Statements ..
      TOL = MAX(TOL1,10.0D0*X02AJF())
      RELERR = MAX(RELER,0.0D0)
      BIG = X02ALF()
      MP1 = M + 1
      NP1 = N + 1
      NP2 = N + 2
      NP3 = N + 3
      IER = 0
      IF (NDIM.GE.NP3 .AND. MDIM.GE.MP1 .AND. M.GE.N .AND. N.GE.1)
     *    GO TO 20
      IER = 3
      GO TO 960
   20 CONTINUE
      NP1MR = 1
      IRANK = N
      RELTMP = RELERR
      RELERR = ZERO
      DO 40 J = 1, M
         A(NP1,J) = ONE
         A(NP2,J) = -B(J)
         A(NP3,J) = N + J
   40 CONTINUE
      A(NP1,MP1) = ZERO
      ITER = 0
      DO 60 I = 1, N
         X(I) = ZERO
         A(I,MP1) = I
   60 CONTINUE
C     LEVEL 1.
      LEV = 1
      K = 0
   80 K = K + 1
      NP1MK = NP1 - K
      MODE = 0
      DO 100 J = K, M
         B(J) = ONE
  100 CONTINUE
C     DETERMINE THE VECTOR TO ENTER THE BASIS.
  120 D = -BIG
      DO 140 J = K, M
         IF (B(J).EQ.ZERO) GO TO 140
         DD = ABS(A(NP2,J))
         IF (DD.LE.D) GO TO 140
         PCOL = J
         D = DD
  140 CONTINUE
      IF (K.GT.1) GO TO 160
C     TEST FOR ZERO RIGHT-HAND SIDE.
      IF (D.GT.TOL) GO TO 160
      RESMAX = ZERO
      MODE = 2
      GO TO 780
C     DETERMINE THE VECTOR TO LEAVE THE BASIS.
  160 D = TOL
      DO 180 I = 1, NP1MK
         DD = ABS(A(I,PCOL))
         IF (DD.LE.D) GO TO 180
         PROW = I
         D = DD
  180 CONTINUE
      IF (D.GT.TOL) GO TO 700
C     CHECK FOR LINEAR DEPENDENCE IN LEVEL 1.
      B(PCOL) = ZERO
      IF (MODE.EQ.1) GO TO 120
      DO 220 J = K, M
         IF (B(J).EQ.ZERO) GO TO 220
         DO 200 I = 1, NP1MK
            IF (ABS(A(I,J)).LE.TOL) GO TO 200
            MODE = 1
            GO TO 120
  200    CONTINUE
  220 CONTINUE
      IRANK = K - 1
      NP1MR = NP1 - IRANK
      IER = 1
      GO TO 340
  240 IF (PCOL.EQ.K) GO TO 280
C     INTERCHANGE COLUMNS IN LEVEL 1.
      DO 260 I = 1, NP3
         D = A(I,PCOL)
         A(I,PCOL) = A(I,K)
         A(I,K) = D
  260 CONTINUE
  280 IF (PROW.EQ.NP1MK) GO TO 320
C     INTERCHANGE ROWS IN LEVEL 1.
      DO 300 J = 1, MP1
         D = A(PROW,J)
         A(PROW,J) = A(NP1MK,J)
         A(NP1MK,J) = D
  300 CONTINUE
  320 IF (K.LT.N) GO TO 80
  340 IF (IRANK.EQ.M) GO TO 780
      IRNKP1 = IRANK + 1
C     LEVEL 2.
      LEV = 2
C     DETERMINE THE VECTOR TO ENTER THE BASIS.
      D = TOL
      DO 360 J = IRNKP1, M
         DD = ABS(A(NP2,J))
         IF (DD.LE.D) GO TO 360
         PCOL = J
         D = DD
  360 CONTINUE
C     COMPARE CHEBYSHEV ERROR WITH TOL.
      IF (D.GT.TOL) GO TO 380
      RESMAX = ZERO
      MODE = 3
      GO TO 780
  380 IF (A(NP2,PCOL).LT.-TOL) GO TO 420
      A(NP1,PCOL) = TWO - A(NP1,PCOL)
      DO 400 I = NP1MR, NP3
         IF (I.EQ.NP1) GO TO 400
         A(I,PCOL) = -A(I,PCOL)
  400 CONTINUE
C     ARRANGE FOR ALL ENTRIES IN PIVOT COLUMN
C     (EXCEPT PIVOT) TO BE NEGATIVE.
  420 IF (NP1MR.GT.N) GO TO 480
      DO 460 I = NP1MR, N
         IF (A(I,PCOL).LT.TOL) GO TO 460
         DO 440 J = 1, M
            A(NP1,J) = A(NP1,J) + TWO*A(I,J)
            A(I,J) = -A(I,J)
  440    CONTINUE
         A(I,MP1) = -A(I,MP1)
  460 CONTINUE
  480 PROW = NP1
      GO TO 700
  500 IF (IRNKP1.EQ.M) GO TO 780
      IF (PCOL.EQ.M) GO TO 540
C     INTERCHANGE COLUMNS IN LEVEL 2.
      DO 520 I = NP1MR, NP3
         D = A(I,PCOL)
         A(I,PCOL) = A(I,M)
         A(I,M) = D
  520 CONTINUE
  540 MM1 = M - 1
C     LEVEL 3.
      LEV = 3
C     DETERMINE THE VECTOR TO ENTER THE BASIS.
  560 D = -TOL
      VAL = TWO*A(NP2,M)
      DO 600 J = IRNKP1, MM1
         IF (A(NP2,J).GE.D) GO TO 580
         PCOL = J
         D = A(NP2,J)
         MODE = 0
         GO TO 600
  580    DD = VAL - A(NP2,J)
         IF (DD.GE.D) GO TO 600
         MODE = 1
         PCOL = J
         D = DD
  600 CONTINUE
      IF (D.GE.-TOL) GO TO 780
      DD = -D/A(NP2,M)
      IF (DD.GE.RELTMP) GO TO 620
      RELERR = DD
      MODE = 4
      GO TO 780
  620 IF (MODE.EQ.0) GO TO 660
      DO 640 I = NP1MR, NP1
         A(I,PCOL) = TWO*A(I,M) - A(I,PCOL)
  640 CONTINUE
      A(NP2,PCOL) = D
      A(NP3,PCOL) = -A(NP3,PCOL)
C     DETERMINE THE VECTOR TO LEAVE THE BASIS.
  660 D = BIG
      DO 680 I = NP1MR, NP1
         IF (A(I,PCOL).LE.TOL) GO TO 680
         DD = A(I,M)/A(I,PCOL)
         IF (DD.GE.D) GO TO 680
         PROW = I
         D = DD
  680 CONTINUE
      IF (D.LT.BIG) GO TO 700
      IER = 2
      GO TO 780
C     PIVOT ON A(PROW,PCOL).
  700 PIVOT = A(PROW,PCOL)
      DO 720 J = 1, M
         A(PROW,J) = A(PROW,J)/PIVOT
  720 CONTINUE
      DO 740 J = 1, M
         IF (J.EQ.PCOL) GO TO 740
         CALL E02GCZ(A(1,J),A(1,PCOL),A(PROW,J),PROW,NP1MR,NP2)
  740 CONTINUE
      TPIVOT = -PIVOT
      DO 760 I = NP1MR, NP2
         A(I,PCOL) = A(I,PCOL)/TPIVOT
  760 CONTINUE
      A(PROW,PCOL) = ONE/PIVOT
      D = A(PROW,MP1)
      A(PROW,MP1) = A(NP3,PCOL)
      A(NP3,PCOL) = D
      ITER = ITER + 1
      GO TO (240,500,560) LEV
C     PREPARE OUTPUT.
  780 DO 800 J = 1, M
         B(J) = ZERO
  800 CONTINUE
      IF (MODE.EQ.2) GO TO 940
      IF (IRANK.EQ.0) GO TO 840
      DO 820 J = 1, IRANK
         K = A(NP3,J)
         X(K) = A(NP2,J)
  820 CONTINUE
  840 IF (MODE.EQ.3 .OR. IRANK.EQ.M) GO TO 940
      DO 860 I = NP1MR, NP1
         K = ABS(A(I,MP1)) - DBLE(N)
         B(K) = A(NP2,M)*SIGN(ONE,A(I,MP1))
  860 CONTINUE
      IF (IRNKP1.EQ.M) GO TO 900
      DO 880 J = IRNKP1, MM1
         K = ABS(A(NP3,J)) - DBLE(N)
         B(K) = (A(NP2,M)-A(NP2,J))*SIGN(ONE,A(NP3,J))
  880 CONTINUE
C     TEST FOR NON-UNIQUE SOLUTION.
  900 DO 920 I = NP1MR, NP1
         IF (ABS(A(I,M)).GT.TOL) GO TO 920
         IER = 1
         GO TO 940
  920 CONTINUE
  940 IF (MODE.NE.2 .AND. MODE.NE.3) RESMAX = A(NP2,M)
      IF (IRANK.EQ.M) RESMAX = ZERO
      IF (MODE.EQ.4) RESMAX = RESMAX - D
  960 CONTINUE
      RELER = RELERR
      IFAIL = P01ABF(IFAIL,IER,SRNAME,0,P01REC)
      RETURN
      END
