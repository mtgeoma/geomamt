      SUBROUTINE D02SAU(P,M,A,E,N,X,NPOINT,H,HMAX,HMIN,FCN,FCN1,FCN2,BC,
     *                  BC1,YMAX,W,IW,PRSOL,C,M1,IFLAG)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 8A REVISED. IER-260 (AUG 1980)
C     MARK 9 REVISED. IER-310 (SEP 1981).
C     MARK 9A REVISED. IER-355 (NOV 1981)
C     MARK 11 REVISED. IER-419 (FEB 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     CALCULATES SOLUTION
C     BC1, BC, FCN1, FCN2, FCN, PRSOL
C     .. Scalar Arguments ..
      DOUBLE PRECISION  YMAX
      INTEGER           IFLAG, IW, M, M1, N, NPOINT
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IW,2), C(IW,M1), E(N), H(NPOINT),
     *                  HMAX(NPOINT), HMIN(NPOINT), P(M), W(IW,11),
     *                  X(NPOINT)
C     .. Subroutine Arguments ..
      EXTERNAL          BC, BC1, FCN, FCN1, FCN2, PRSOL
C     .. Scalars in Common ..
      DOUBLE PRECISION  COUT1, COUT2, SQEPS, W2, W3
      INTEGER           ICASE, II, IW2
C     .. Arrays in Common ..
      DOUBLE PRECISION  W1(4)
      INTEGER           IW1(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  S, STEP, T, T1, T2
      INTEGER           I, IFAIL, J, N2, NP
C     .. Local Arrays ..
      DOUBLE PRECISION  CIN(6), COMM(5), CON(3), COUT(14)
C     .. External Subroutines ..
      EXTERNAL          D02HAW, D02KDY, D02SAY, D02XAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SIGN, DBLE
C     .. Common blocks ..
      COMMON            /AD02SA/W1, SQEPS, W2, W3, IW1, IW2, ICASE
      COMMON            /BD02SA/COUT1, COUT2, II
C     .. Executable Statements ..
      IFLAG = 0
      IF (IW.GE.N .AND. IW.GE.M .AND. NPOINT.GE.2 .AND. M1.GE.1)
     *    GO TO 40
C     INPUT ERROR
   20 IFLAG = 1
      RETURN
   40 IF (X(1).EQ.X(NPOINT)) GO TO 20
      IF (ICASE.GT.1 .AND. M1.EQ.1) RETURN
      NP = NPOINT - 1
      S = SIGN(1.D0,X(NPOINT)-X(1))
      DO 60 J = 1, NP
         IF (S*(X(J+1)-X(J)).LE.0.D0) GO TO 20
   60 CONTINUE
      T = X(1)
      GO TO (80,100,120) ICASE
   80 CALL BC1(W(1,1),W(1,2),P,M,N)
      GO TO 160
  100 CALL BC(W(1,1),W(1,2),P)
      GO TO 160
  120 DO 140 J = 1, IW
         W(J,1) = A(J,1)
  140 CONTINUE
  160 IF (YMAX.EQ.0.D0) GO TO 200
      DO 180 J = 1, N
         IF (ABS(W(J,1)).GE.YMAX .OR. ABS(W(J,2)).GE.YMAX) GO TO 20
  180 CONTINUE
  200 IF (ICASE.GT.1) GO TO 220
      CALL PRSOL(T,W(1,1),N)
      GO TO 260
  220 DO 240 J = 1, N
         C(J,1) = W(J,1)
  240 CONTINUE
      STEP = (X(NPOINT)-X(1))/DBLE(M1-1)
      T = T + STEP
      N2 = 2
  260 IF (T.EQ.X(1)) RETURN
      IF (T.EQ.X(NPOINT)) GO TO 280
      S = SIGN(1.D0,X(NPOINT)-X(1))
      IF (SIGN(1.D0,T-X(NPOINT)).EQ.S) RETURN
      IF (SIGN(1.D0,T-X(1)).NE.S) RETURN
  280 COMM(1) = 0.D0
      COMM(2) = YMAX
      COMM(3) = 0.D0
      DO 300 J = 1, 3
         CON(J) = 0.D0
  300 CONTINUE
      CIN(2) = 3.D0
      DO 320 J = 1, N
         W(J,9) = E(J)
         W(J,10) = SQEPS
  320 CONTINUE
      DO 760 I = 1, NP
         COMM(4) = 0.D0
         IF (T.EQ.X(I+1)) GO TO 340
         IF (SIGN(1.D0,T-X(I)).EQ.SIGN(1.D0,T-X(I+1))) GO TO 340
         COMM(4) = -1.D0
         COMM(5) = T
  340    CIN(1) = 1.D0
         CIN(3) = HMIN(I)
         CIN(4) = HMAX(I)
         CIN(5) = H(I)
         T1 = X(I)
  360    IFAIL = 1
         IF (ICASE.EQ.3) GO TO 380
         II = I
         CALL D02KDY(T1,X(I+1),N,W(1,1),CIN,1.D0,D02SAY,COMM,CON,COUT,
     *               W(1,3),IW,9,FCN1,FCN2,P,M,IFAIL)
         GO TO 400
  380    CALL D02KDY(T1,X(I+1),N,W(1,1),CIN,1.D0,D02HAW,COMM,CON,COUT,
     *               W(1,3),IW,9,FCN,FCN,P,M,IFAIL)
  400    IF (IFAIL.EQ.0) GO TO 480
         W(1,1) = T1
         GO TO (420,440,440,460,420,420,420) IFAIL
  420    IFLAG = 6
         RETURN
  440    IFLAG = 3
         W(1,2) = T1
         RETURN
  460    IFLAG = 4
         W(1,2) = T1
         RETURN
  480    IF (CIN(1).EQ.2.D0 .AND. T1.EQ.X(I+1)) GO TO 520
         IF (CIN(1).EQ.6.D0 .AND. T1.NE.X(I+1)) GO TO 640
         IF (CIN(1).NE.4.D0) GO TO 420
C        SOLUTION TOO LARGE
  500    W(1,2) = T1
         IFLAG = 5
         RETURN
  520    IF (YMAX.EQ.0.0D0) GO TO 560
         DO 540 J = 1, N
            IF (ABS(W(J,1)).GE.YMAX) GO TO 500
  540    CONTINUE
  560    IF (T.NE.T1) GO TO 620
  580    DO 600 J = 1, N
            W(J,11) = W(J,1)
  600    CONTINUE
         GO TO 660
  620    IF (SIGN(1.D0,T-T1).EQ.SIGN(1.D0,T-X(I))) GO TO 760
  640    IFAIL = 1
         CALL D02XAF(T,T1,COUT,N,W(1,1),W(1,3),IW,W(1,11),IFAIL)
         IF (IFAIL.EQ.0) GO TO 660
C        IMPOSS ERROR IN D02XAF
         IFLAG = 2
         RETURN
  660    T2 = T
         IF (ICASE.GT.1) GO TO 680
         CALL PRSOL(T,W(1,11),N)
         GO TO 720
  680    DO 700 J = 1, N
            C(J,N2) = W(J,11)
  700    CONTINUE
         IF (M1.EQ.2) RETURN
         N2 = N2 + 1
         T = T + STEP
         IF (N2.EQ.M1) T = X(NPOINT)
  720    IF (T.EQ.T2 .OR. T2.EQ.X(NPOINT)) RETURN
         IF (T.EQ.X(NPOINT)) GO TO 740
         IF (SIGN(1.D0,T-X(NPOINT)).EQ.S) RETURN
         IF (SIGN(1.D0,T-T2).NE.S) RETURN
  740    IF (T.EQ.T1) GO TO 580
         IF (SIGN(1.D0,T-T1).NE.SIGN(1.D0,T-COUT(5))) GO TO 640
         IF (T.EQ.X(I+1)) GO TO 360
         IF (SIGN(1.D0,T-X(I)).EQ.SIGN(1.D0,T-X(I+1))) GO TO 360
         COMM(4) = -1.D0
         COMM(5) = T
         IF (CIN(1).EQ.6.D0) GO TO 360
  760 CONTINUE
      RETURN
      END
