      SUBROUTINE D02BDF(X,XEND,N,Y,TOL,IRELAB,FCN,STIFF,YNORM,W,IW,M,
     *                  OUTPUT,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     CALCULATES A GLOBAL ERROR ESTIMATE FOR A SYSTEM
C     OF ODES, AND MAKES A STIFFNESS CHECK IF REQUIRED.
C     FCN, OUTPUT
C
C     PARAMETER TESTS
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02BDF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  STIFF, TOL, X, XEND, YNORM
      INTEGER           IFAIL, IRELAB, IW, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  W(N,IW), Y(N)
C     .. Subroutine Arguments ..
      EXTERNAL          FCN, OUTPUT
C     .. Local Scalars ..
      DOUBLE PRECISION  C, C2, H, OLDX, S, S1, ST2, ST3, STIFFM, TOL1,
     *                  U, Z
      INTEGER           I, IND, J, K1, N3, Q
      LOGICAL           STEP1, STIFF1, STIFF2
C     .. Local Arrays ..
      DOUBLE PRECISION  CIN(6), COMM(5), CONST(3), COUT(14)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02PAF, D02YAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SIGN
C     .. Executable Statements ..
      IF (N.GT.0 .AND. TOL.GT.0.D0 .AND. IRELAB.GE.0 .AND. IRELAB.LE.
     *    2 .AND. IW.GE.2 .AND. YNORM.GE.0.D0) GO TO 20
C     INPUT ERROR
      IND = 1
      GO TO 920
   20 IF (X.NE.XEND) GO TO 80
      DO 60 I = 1, N
         DO 40 J = 1, 3
            W(I,J) = 0.D0
   40    CONTINUE
   60 CONTINUE
      IND = 0
      GO TO 920
   80 IF (IW.GE.14) GO TO 100
      IF (STIFF.EQ.0.D0 .AND. IW.GE.12) GO TO 100
      IF (STIFF.LT.0.D0 .AND. IW.GE.13) GO TO 100
      IND = 1
      GO TO 920
  100 IF (YNORM.EQ.0.D0) GO TO 140
      DO 120 I = 1, N
         IF (ABS(Y(I)).LT.YNORM) GO TO 120
C        NORM OF Y TOO LARGE INITIALLY
         IND = 10
         GO TO 920
  120 CONTINUE
C
C     INITIALISATION
C
  140 TOL1 = 16.D0*TOL
      Q = 0
      STIFF2 = STIFF .EQ. 0.D0
      STIFF1 = STIFF .LT. 0.D0
      STIFF = 0.D0
      IF (STIFF2) GO TO 180
      ST2 = 0.D0
      ST3 = 0.D0
      STIFFM = 0.D0
      OLDX = X
      DO 160 I = 1, N
         W(I,IW) = Y(I)
  160 CONTINUE
  180 CIN(1) = 1.D0
      CIN(2) = 0.D0
      IF (IRELAB.EQ.1) CIN(2) = 1.D0
      IF (IRELAB.EQ.2) CIN(2) = 2.D0
      CIN(3) = 0.D0
      CIN(4) = 0.D0
      CIN(5) = 0.D0
      CONST(1) = 0.D0
      CONST(2) = 0.D0
      CONST(3) = 0.D0
      COMM(1) = 0.D0
      COMM(2) = YNORM
      COMM(3) = 0.D0
      COMM(4) = 1.D0
      DO 200 I = 1, N
         W(I,1) = 0.D0
         W(I,2) = 0.D0
         W(I,3) = -1.D0
         W(I,4) = Y(I)
         W(I,5) = 0.D0
  200 CONTINUE
      STEP1 = .TRUE.
C
C     COMPUTE SOLUTION INTERRUPTING AT EACH STEP
C
  220 IND = 1
      CALL D02PAF(X,XEND,N,W(1,4),CIN,TOL1,FCN,COMM,CONST,COUT,W(1,6)
     *            ,N,IW-5,IND)
      IF (IND.EQ.0) GO TO 260
      IF (IND.GE.2 .AND. IND.LE.4) GO TO 920
  240 IND = 8
      GO TO 920
  260 IF (CIN(1).NE.4.D0) GO TO 300
C     NORM OF SOLUTION TOO LARGE
  280 IND = 7
      GO TO 920
  300 IF (X.EQ.XEND .AND. CIN(1).NE.2.D0) GO TO 240
      IF (CIN(1).NE.5.D0 .AND. X.NE.XEND) GO TO 240
      IF ( .NOT. STEP1) GO TO 320
      N3 = 7
      Z = COUT(5)
      H = 0.5D0*CIN(5)
      C2 = COUT(12)/(COUT(11)*COUT(11))
      GO TO 340
C
C     COMPUTE TWO STEP SOLUTION
C
  320 Z = COUT(4)
      H = 0.5D0*(X-COUT(4))
  340 CALL FCN(Z,Y,W(1,9))
      CALL D02YAF(Z,H,N,Y,FCN,W(1,9),N,4)
      IF (YNORM.EQ.0.D0) GO TO 380
      DO 360 I = 1, N
         IF (ABS(Y(I)).LT.YNORM) GO TO 360
C        NORM OF SOLUTION TOO LARGE
         IND = 9
         GO TO 920
  360 CONTINUE
  380 Z = Z + H
      CALL FCN(Z,Y,W(1,9))
      CALL D02YAF(Z,H,N,Y,FCN,W(1,9),N,4)
      IF (YNORM.EQ.0.D0) GO TO 420
      DO 400 I = 1, N
         IF (ABS(Y(I)).LT.YNORM) GO TO 400
         IND = 9
         GO TO 920
  400 CONTINUE
  420 IF (STIFF2 .OR. STEP1) GO TO 620
C
C     COMPUTE STIFFNESS FACTORS
C
      S = 1.D0/(COUT(8)+100.D0)
      IF (S.GT.COUT(11)) GO TO 440
C     COUT(8) TOO LARGE
      IND = 5
      GO TO 920
  440 Z = COUT(4)
      CALL D02YAF(Z,H,N,W(1,7),FCN,W(1,8),N,6)
      IF (YNORM.EQ.0.D0) GO TO 480
      DO 460 I = 1, N
         IF (ABS(W(I,7)).GE.YNORM) GO TO 280
  460 CONTINUE
  480 S1 = 0.D0
      DO 500 I = 1, N
         S1 = MAX(S1,ABS(W(I,11)))
  500 CONTINUE
      IF (IRELAB.EQ.1) GO TO 540
      U = 0.D0
      DO 520 I = 1, N
         U = MAX(U,ABS(W(I,4)))
  520 CONTINUE
  540 K1 = IRELAB + 1
      ST2 = ST2 + 1.D0
      GO TO (560,580,600) K1
  560 IF (S1.LE.TOL1*MAX(1.D0,U)) ST3 = ST3 + 1.D0
      GO TO 620
  580 IF (S1.LE.TOL1) ST3 = ST3 + 1.D0
      GO TO 620
  600 IF (S1.LE.TOL1*MAX(C2,U)) ST3 = ST3 + 1.D0
C
C     COMPUTE ERROR ESTIMATE
C
  620 DO 680 I = 1, N
         W(I,1) = (Y(I)-W(I,N3))/15.D0
         IF (W(I,1).EQ.0.D0) GO TO 640
         IF (ABS(W(I,1)).GT.ABS(W(I,2))) W(I,2) = W(I,1)
         IF (W(I,5).EQ.0.D0) GO TO 660
         IF (SIGN(1.D0,W(I,5)).NE.SIGN(1.D0,W(I,1))) GO TO 640
         IF (ABS(W(I,1)).LT.ABS(W(I,5)) .AND. W(I,3).EQ.-1.D0) W(I,3)
     *       = 0.D0
         GO TO 660
  640    IF (W(I,3).GE.0.D0) W(I,3) = W(I,3) + 1.D0
         IF (W(I,3).EQ.-1.D0) W(I,3) = 1.D0
  660    W(I,5) = W(I,1)
  680 CONTINUE
      N3 = 4
      IF (STIFF2) GO TO 700
      STIFF = 0.D0
      IF (STEP1) GO TO 700
      STIFF = ST3/ST2
  700 IF (M.LE.0) GO TO 720
      Q = Q + 1
      IF (Q.LT.M) GO TO 720
      Z = X
      IF (STEP1) Z = COUT(4)
      CALL OUTPUT(Z,Y,W,STIFF)
      Q = 0
  720 IF (X.EQ.XEND .AND. .NOT. STEP1) GO TO 760
      IF ( .NOT. STIFF1) GO TO 740
      IF (STIFF.GT.0.9D0 .AND. COUT(8).GT.100.D0) GO TO 940
      IF (COUT(8).LT.1000.D0) GO TO 740
      COUT(8) = 0.D0
      COUT(9) = 0.D0
      COUT(3) = 0.D0
      ST2 = 0.D0
      ST3 = 0.D0
      STIFFM = MAX(STIFFM,STIFF)
      GO TO 220
  740 IF ( .NOT. STEP1) GO TO 220
      STEP1 = .FALSE.
      GO TO 320
  760 IF (3.D0*COUT(3).LT.COUT(8)) GO TO 780
      IND = 6
      GO TO 920
  780 IF (STIFF2) GO TO 940
      IF (COUT(8).LE.100.D0) GO TO 920
C
C     COMPUTE SECOND STIFFNESS FACTOR
C
      STIFF = STIFFM
      IF (STIFF1) GO TO 940
      C = COUT(8)
      X = OLDX
      CIN(1) = 1.D0
      CIN(5) = 0.5D0*CIN(5)
      COMM(4) = 0.D0
      IND = 1
      CALL D02PAF(X,XEND,N,W(1,IW),CIN,TOL,FCN,COMM,CONST,COUT,W(1,6)
     *            ,N,IW-5,IND)
      IF (IND.EQ.0) GO TO 820
      IND = IND + 10
      IF (IND.GE.12 .AND. IND.LE.14) GO TO 920
  800 IND = 18
      GO TO 920
  820 IF (CIN(1).NE.4.D0) GO TO 840
      IND = 17
      GO TO 920
  840 IF (X.NE.XEND .OR. CIN(1).NE.2.D0) GO TO 800
      IF (3.D0*COUT(3).LT.COUT(8)) GO TO 860
      IND = 16
      GO TO 920
  860 IF (COUT(8).GT.100.D0) GO TO 880
      GO TO 920
  880 S = 1.D0/(COUT(8)+100.D0)
      IF (S.GT.COUT(11)) GO TO 900
      IND = 15
      GO TO 920
  900 STIFF = MAX(0.D0,(2.D0*C-COUT(8))/COUT(8))
      STIFF = STIFF + ST3/ST2
      STIFF = STIFF*0.5D0
      GO TO 940
  920 STIFF = 0.D0
  940 IFAIL = P01ABF(IFAIL,IND,SRNAME,0,P01REC)
      RETURN
      END
