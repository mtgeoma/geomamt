C     MARK 11.5(F77) REVISED. (SEPT 1985.)
      SUBROUTINE D02SAZ(F,P,M,A,E,N,N1,X,NPOINT,H,HMAX,HMIN,MARK,FCN,
     *                  FCN1,FCN2,EQN,BC,BC1,YMAX,W,IW,IFLAG)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 8A REVISED. IER-261 (AUG 1980).
C     MARK 8F REVISED. IER-299 (APR 1981).
C     MARK 11 REVISED. IER-419 (FEB 1984).
C     CALCULATES RESIDUAL
C     BC1, BC, EQN, FCN1, FCN2, FCN
C     .. Scalar Arguments ..
      DOUBLE PRECISION  YMAX
      INTEGER           IFLAG, IW, M, MARK, N, N1, NPOINT
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IW,2), E(N), F(M), H(NPOINT), HMAX(NPOINT),
     *                  HMIN(NPOINT), P(M), W(IW,11), X(NPOINT)
C     .. Subroutine Arguments ..
      EXTERNAL          BC, BC1, EQN, FCN, FCN1, FCN2
C     .. Scalars in Common ..
      DOUBLE PRECISION  COUT12, COUT13, SQEPS, W2, W3
      INTEGER           COUNT, ICASE, IFAIL1, II, IW2
C     .. Arrays in Common ..
      DOUBLE PRECISION  W1(4)
      INTEGER           IW1(2)
C     .. Local Scalars ..
      DOUBLE PRECISION  T
      INTEGER           I, I1, IFAIL, J, J1, NP
      LOGICAL           CH
C     .. Local Arrays ..
      DOUBLE PRECISION  CIN(6), COMM(5), CON(3), COUT(14)
C     .. External Subroutines ..
      EXTERNAL          D02HAW, D02KDY, D02SAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE
C     .. Common blocks ..
      COMMON            /AD02SA/W1, SQEPS, W2, W3, IW1, COUNT, IFAIL1,
     *                  IW2, ICASE
      COMMON            /BD02SA/COUT12, COUT13, II
C     .. Executable Statements ..
      CH = .FALSE.
      IF (IFLAG.EQ.1) CH = .TRUE.
      IFLAG = 0
      COUNT = COUNT + 1
      IF (IW.GE.N .AND. IW.GE.M .AND. M.GE.N1 .AND. N1.GT.0 .AND.
     *    NPOINT.GE.2 .AND. N.GE.N1) GO TO 20
C     INPUT  ERROR
      IFLAG = 1
      RETURN
C     SET UP FOR INTEGRATION
   20 HMAX(NPOINT) = 0.D0
      H(NPOINT) = 0.D0
      GO TO (40,60,80) ICASE
   40 CALL BC1(W(1,1),W(1,2),P,M,N)
      GO TO 140
   60 CALL BC(W(1,1),W(1,2),P)
      GO TO 140
   80 DO 120 I1 = 1, N
         DO 100 J = 1, 2
            W(I1,J) = A(I1,J)
  100    CONTINUE
  120 CONTINUE
  140 IF (YMAX.EQ.0.D0) GO TO 180
      DO 160 I1 = 1, N
         IF (ABS(W(I1,1)).LT.YMAX .AND. ABS(W(I1,2)).LT.YMAX)
     *       GO TO 160
C        SOLUTION TOO LARGE
         IFLAG = 4
         H(NPOINT) = DBLE(NPOINT)
         RETURN
  160 CONTINUE
  180 DO 200 J = 1, N
         W(J,9) = E(J)
         W(J,10) = SQEPS
  200 CONTINUE
      NP = NPOINT - 1
      DO 340 I = 1, NP
         DO 220 J = 1, 3
            CON(J) = 0.D0
  220    CONTINUE
         COMM(1) = 0.D0
         COMM(2) = YMAX
         COMM(3) = 0.D0
         COMM(4) = 0.D0
         CIN(1) = 1.D0
         IF (MARK.EQ.0) GO TO 240
         CIN(1) = 7.D0
         COUT(12) = COUT12
  240    CIN(2) = 3.D0
         CIN(3) = HMIN(I)
         CIN(4) = HMAX(I)
         CIN(5) = H(I)
         T = X(I)
         IFAIL = 1
C        INTEGRATE
         IF (ICASE.EQ.3) GO TO 260
         II = I
         CALL D02KDY(T,X(I+1),N,W(1,1),CIN,1.D0,D02SAY,COMM,CON,COUT,
     *               W(1,3),IW,9,FCN1,FCN2,P,M,IFAIL)
         GO TO 280
  260    CALL D02KDY(T,X(I+1),N,W(1,1),CIN,1.D0,D02HAW,COMM,CON,COUT,
     *               W(1,3),IW,9,FCN,FCN,P,M,IFAIL)
  280    IF (IFAIL.GT.0) GO TO 380
         IF (CIN(1).EQ.4.D0) GO TO 460
         IF (CIN(1).NE.2.D0) GO TO 400
         IF (YMAX.EQ.0.0D0) GO TO 320
         DO 300 I1 = 1, N
            IF (ABS(W(I1,1)).GE.YMAX) GO TO 460
  300    CONTINUE
  320    IF (MARK.NE.0) GO TO 340
         COUT12 = COUT(12)
         IF (CH) H(I) = CIN(5)
         IF (COUT(3).GT.0.3D0*COUT(8)) HMAX(NPOINT) = DBLE(I)
  340 CONTINUE
C     CALCULATE RESIDUAL FROM INTEGRATED VALUES
      J = N - N1
      DO 360 I1 = 1, N1
         J1 = J + I1
         F(I1) = W(J1,1) - W(J1,2)
  360 CONTINUE
      J = M - N1
      I1 = J
      IF (I1.EQ.0) I1 = 1
      IF (J.EQ.0) RETURN
      CALL EQN(F(N1+1),I1,P,M)
      RETURN
C     ERROR EXITS
  380 H(NPOINT) = DBLE(I)
      W(1,2) = T
      GO TO (400,420,420,440,400,400,400) IFAIL
C     IMPOSSIBLE EXIT FROM D02KDY OR D02PAF
  400 IFLAG = 5
      RETURN
C     STEP TOO SMALL
  420 IFLAG = 2
      RETURN
C     INITIAL STEP TOO SMALL
  440 IFLAG = 3
      RETURN
C     SOLUTION TOO LARGE
  460 IFLAG = 4
      W(1,2) = T
      H(NPOINT) = DBLE(I)
      RETURN
      END
