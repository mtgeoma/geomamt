      SUBROUTINE D02SAW(D,U,V,F,P,PF,M,A,B,DP,E,N,N1,W1,W2,NPOINT,X,X1,
     *                  H,HMAX,HMIN,FCN,FCN1,FCN2,EQN,BC,BC1,RANGE,
     *                  RANGE1,CONSTR,YMAX,W,IW,IW2,IFLAG)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 8F REVISED. IER-299 (APR 1981).
C     MARK 11 REVISED. IER-419 (FEB 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 15 REVISED. IER-902 (APR 1991).
C     CALCULATES JACOBIAN
C     BC1, BC, EQN, FCN1, FCN2, FCN, RANGE1, RANGE
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X, X1, YMAX
      INTEGER           IFLAG, IW, IW2, M, N, N1, NPOINT
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IW,2), B(IW,2), D(M), DP(M), E(N), F(M),
     *                  H(NPOINT), HMAX(NPOINT), HMIN(NPOINT), P(M),
     *                  PF(M), U(IW,M), V(IW,M), W(IW,IW2), W1(NPOINT),
     *                  W2(NPOINT)
C     .. Function Arguments ..
      LOGICAL           CONSTR
      EXTERNAL          CONSTR
C     .. Subroutine Arguments ..
      EXTERNAL          BC, BC1, EQN, FCN, FCN1, FCN2, RANGE, RANGE1
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPS, SQEPS, W3, W5
      INTEGER           ICASE, IW3
C     .. Arrays in Common ..
      DOUBLE PRECISION  W4(3)
      INTEGER           IW1(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  S, SQEPS1, T
      INTEGER           I, J, J1, K, MARK, NP
C     .. External Subroutines ..
      EXTERNAL          D02SAT, D02SAZ, F02WEF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SIGN, DBLE
C     .. Common blocks ..
      COMMON            /AD02SA/W4, EPS, SQEPS, W3, W5, IW1, IW3, ICASE
C     .. Executable Statements ..
      SQEPS1 = MAX(MIN(SQEPS,1.D-3),1.D-5)
      I = M + 15
      J = 2*M + 4
      IF (N.GE.N1 .AND. M.GE.N1 .AND. N1.GT.0 .AND. NPOINT.GE.2 .AND.
     *    IW.GE.N .AND. IW.GE.M .AND. IW2.GE.J .AND. IW2.GE.I)
     *    GO TO 20
C     INPUT ERROR
      IFLAG = 1
      RETURN
C     SET UP INITIAL VALUES
   20 T = MIN(1000.D0*EPS,0.1D0*SQEPS)
      NP = NPOINT
      IF (ICASE.EQ.3) CALL D02SAT(P,A,B,IW,0)
      DO 40 I = 1, M
         W(I,1) = P(I)
         W(I,3) = P(I)
         IF (DP(I).EQ.0.D0) DP(I) = SQEPS1*MAX(PF(I),ABS(P(I)))
   40 CONTINUE
C     SET UP JACOBIAN COLUMN-BY-COLUMN
C     FORWARD DIFFERENCES
      DO 680 I = 1, M
         W(I,1) = P(I) + DP(I)
         MARK = 0
C        CHECK CONSTRAINTS
         IF (CONSTR(W(1,1),M)) GO TO 80
C        BACKWARD DIFFERENCES
   60    W(I,1) = P(I) - DP(I)
         MARK = 1
         IF (CONSTR(W(1,1),M)) GO TO 80
C        CONSTRAINTS VIOLATED
         IFLAG = 2
         GO TO 640
C        CALCULATE RANGE
   80    GO TO (100,120,140) ICASE
  100    CALL RANGE1(W1,NP,W(1,1),M)
         GO TO 160
  120    CALL RANGE(W1(1),W1(2),W(1,1))
         NP = 2
         GO TO 160
  140    W1(1) = X
         W1(2) = X1
         NP = 2
         CALL D02SAT(W(1,1),A,B,IW,1)
  160    IF (W1(NP).EQ.W1(1)) GO TO 200
         IF (NP.EQ.2) GO TO 220
         S = SIGN(1.D0,W1(NP)-W1(1))
         DO 180 J = 2, NP
            IF (S*(W1(J)-W1(J-1)).LE.0.D0) GO TO 200
  180    CONTINUE
         GO TO 220
  200    IF (MARK.EQ.0) GO TO 60
C        RANGE VIOLATION
         IFLAG = 3
         GO TO 640
C        INTEGRATE
  220    IFLAG = 0
         CALL D02SAZ(W(1,2),W(1,1),M,A,E,N,N1,W1,NP,H,HMAX,HMIN,1,FCN,
     *               FCN1,FCN2,EQN,BC,BC1,YMAX,W(1,5),IW,IFLAG)
         IF (IFLAG.EQ.0) GO TO 240
C        ERROR IN INTEGRATION
         IFLAG = IFLAG + 3
         GO TO 640
  240    K = 0
C        CALCULATE RESIDUAL AND CHECK FOR SIGNIFICANCE
         DO 260 J = 1, M
            W(J,I+15) = W(J,2) - F(J)
            IF (ABS(W(J,I+15)).LT.MAX(ABS(W(J,2)),ABS(F(J)))*T) K = K +
     *          1
            W(J,I+15) = W(J,I+15)/DP(I)
            IF (MARK.EQ.1) W(J,I+15) = -W(J,I+15)
  260    CONTINUE
         IF (K.EQ.M) GO TO 280
         W(I,1) = P(I)
         GO TO 680
C        CENTRAL DIFFERENCES
  280    MARK = 2
         DP(I) = 10.D0*DP(I)
         W(I,1) = P(I) + DP(I)
         W(I,3) = P(I) - DP(I)
C        CHECK CONSTRAINTS
         IF (CONSTR(W(1,1),M) .AND. CONSTR(W(1,3),M)) GO TO 300
C        CONSTRAINT VIOLATION
         IFLAG = 9
         GO TO 640
C        SET UP RANGE
  300    GO TO (320,340,360) ICASE
  320    CALL RANGE1(W1,NP,W(1,1),M)
         GO TO 380
  340    CALL RANGE(W1(1),W1(2),W(1,1))
         NP = 2
         GO TO 380
  360    W1(1) = X
         W1(2) = X1
         NP = 2
         CALL D02SAT(W(1,1),A,B,IW,1)
  380    IF (W1(NP).EQ.W1(1)) GO TO 420
         IF (NP.EQ.2) GO TO 440
         S = SIGN(1.D0,W1(NP)-W1(1))
         DO 400 J = 2, NP
            IF (S*(W1(J)-W1(J-1)).LE.0.D0) GO TO 420
  400    CONTINUE
         GO TO 440
C        RANGE VIOLATION
  420    IFLAG = 10
         GO TO 640
  440    GO TO (460,480,500) ICASE
  460    CALL RANGE1(W2,NP,W(1,3),M)
         GO TO 520
  480    CALL RANGE(W2(1),W2(2),W(1,3))
         NP = 2
         GO TO 520
  500    W2(1) = X
         W2(2) = X1
         NP = 2
  520    IF (W2(NP).EQ.W2(1)) GO TO 420
         IF (NP.EQ.2) GO TO 560
         S = SIGN(1.D0,W2(NP)-W2(1))
         DO 540 J = 2, NP
            IF (S*(W2(J)-W2(J-1)).LE.0.D0) GO TO 420
  540    CONTINUE
C        INTEGRATE
  560    IFLAG = 0
         CALL D02SAZ(W(1,2),W(1,1),M,A,E,N,N1,W1,NP,H,HMAX,HMIN,1,FCN,
     *               FCN1,FCN2,EQN,BC,BC1,YMAX,W(1,5),IW,IFLAG)
         IF (IFLAG.EQ.0) GO TO 580
C        ERROR IN INTEGRATION
         IFLAG = IFLAG + 3
         GO TO 640
  580    IFLAG = 0
         IF (ICASE.EQ.3) CALL D02SAT(W(1,3),A,B,IW,1)
         CALL D02SAZ(W(1,4),W(1,3),M,A,E,N,N1,W2,NP,H,HMAX,HMIN,1,FCN,
     *               FCN1,FCN2,EQN,BC,BC1,YMAX,W(1,5),IW,IFLAG)
         IF (IFLAG.EQ.0) GO TO 600
C        ERROR IN INTEGRATION
         IFLAG = IFLAG + 3
         GO TO 640
C        CALCULATE RESIDUAL AND CHECK FOR SIGNIFICANCE
  600    K = 0
         DO 620 J = 1, M
            W(J,I+15) = W(J,2) - W(J,4)
            IF (ABS(W(J,I+15)).LT.MAX(ABS(W(J,2)),ABS(W(J,4)))*T)
     *          K = K + 1
            W(J,I+15) = W(J,I+15)/(2.D0*DP(I))
  620    CONTINUE
         IF (K.NE.M) GO TO 660
C        COLUMN OF JACOBIAN INSIGNIFICANT
         IFLAG = 11
  640    H(NP) = DBLE(I)
         RETURN
  660    W(I,1) = P(I)
         W(I,3) = P(I)
  680 CONTINUE
C     SHIFT JACOBIAN
      DO 720 I = 1, M
         DO 700 J = 1, M
C Changed because of call to more up to date F02 routine
C      R.W.Brankin, NAG, 20th Feb 1991
C           W(J,I+4) = W(J,I+15)
            U(J,I) = W(J,I+15)
  700    CONTINUE
  720 CONTINUE
C The following section has been has been completely replaced as the
C linear algebra routines called were scheduled for withdrawal at
C Mark 15. The functionality of the whole section can be duplicated by
C a call to F02WEF.
C      R.W.Brankin, NAG, 20th Feb 1991
C
C     CALCULATE SVD (CODE EXTRACTED FROM F02WCF)
C
C      IFLAG = 1
C      CALL F01QAF(M,M,W(1,5),IW,W(1,5),IW,W(1,3),IFLAG)
C      DO 760 J = 1, M
C         MJ = M + J
C         DO 740 I = 1, M
C            W(I,MJ+4) = W(I,J+4)
C  740    CONTINUE
C  760 CONTINUE
C      CALL F02WCZ(M,M,W(1,5),IW,W(1,3),W(1,5),IW)
C      IFLAG = 1
C      CALL F01LZF(M,W(1,M+5),IW,W(1,M+5),IW,.FALSE.,W(1,2)
C     *            ,.FALSE.,.TRUE.,W(1,5),IW,M,.FALSE.,W(1,2),1,1,W(1,1)
C     *            ,W(1,2),W(1,2),W(1,2),IFLAG)
C      CALL F02WAY(M,W(1,M+5),IW,W(1,M+5),IW)
C      IFLAG = 1
C      CALL F02SZF(M,W(1,1),W(1,2),W(1,1),.FALSE.,W(1,2),.TRUE.,W(1,5)
C     *            ,IW,M,.TRUE.,W(1,M+5),IW,M,W(1,2),W(1,3),W(1,4),IFLAG)
C      IF (IFLAG.EQ.0) GO TO 780
CC     F02SZF FAILURE
C      IFLAG = 12
C      RETURN
CC
CC     SET UP SVD AFTER SUCCESSFUL COMPUTATION
CC
C  780 DO 820 J = 1, M
C         J1 = J + M + 4
C         D(J) = W(J,1)
C         DO 800 I = 1, M
C            U(I,J) = W(I,J+4)
C            V(J,I) = W(I,J1)
C  800    CONTINUE
C  820 CONTINUE
      IFLAG = 1
      CALL F02WEF(M,M,U,IW,0,W(1,1),IW,.TRUE.,W(1,1),IW,D,.TRUE.,V,IW,
     *            W(1,1),IFLAG)
      IF (IFLAG.NE.0) IFLAG = 12
      RETURN
      END
