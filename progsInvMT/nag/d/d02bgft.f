      SUBROUTINE D02BGF(X,XEND,N,Y,TOL,HMAX,M,VAL,FCN,W,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     TO FIND POSITION AT WHICH A COMPONENT OF THE SOLUTION ATTAINS
C     A GIVEN VALUE
C     FCN
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02BGF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HMAX, TOL, VAL, X, XEND
      INTEGER           IFAIL, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  W(N,10), Y(N)
C     .. Subroutine Arguments ..
      EXTERNAL          FCN
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B, FX, NEWF, OLDF, TOLX, X1, X2, XS
      INTEGER           I, IND, IND1, IR, J
C     .. Local Arrays ..
      DOUBLE PRECISION  C(6), COMM(5), CON(3), COUT(14), D(17)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          C05AZF, D02PAF, D02XBF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SIGN
C     .. Executable Statements ..
      IF (TOL.GT.0.D0 .AND. N.GT.0 .AND. M.GT.0 .AND. M.LE.N)
     *    GO TO 20
C     INPUT ERROR
      IND = 1
      GO TO 380
   20 IF (Y(M).NE.VAL) GO TO 40
C     SOLUTION IS X
      IND = 0
      GO TO 380
   40 C(1) = 1.D0
      C(2) = 0.D0
      C(3) = 0.D0
      C(4) = HMAX
      C(5) = 0.D0
      COMM(1) = 0.D0
      COMM(2) = 0.D0
      COMM(3) = 1.D0
      COMM(4) = 0.D0
      CON(1) = 0.D0
      CON(2) = 0.D0
      CON(3) = 0.D0
      DO 60 I = 1, N
         W(I,9) = MAX(1.D0,ABS(Y(I)))
         IF (Y(I).GT.0.D0) W(I,9) = -W(I,9)
   60 CONTINUE
      W(M,9) = VAL
   80 IND = 1
      CALL D02PAF(X,XEND,N,Y,C,TOL,FCN,COMM,CON,COUT,W,N,10,IND)
      XS = X
      IF (IND.NE.0) GO TO 100
      IF (C(1).EQ.3.D0 .AND. X.NE.XEND) GO TO 120
      IF (X.EQ.XEND .AND. C(1).EQ.2.D0) GO TO 120
C     ERROR EXITS FROM D02PAF
  100 IF (IND.EQ.3 .OR. IND.EQ.4) IND = IND - 1
      IF (IND.NE.2 .AND. IND.NE.3) IND = 6
      IF (IND.EQ.2) GO TO 360
      GO TO 380
  120 IF (Y(M).NE.VAL) GO TO 140
      IND = 0
      GO TO 360
  140 IF (W(M,2).NE.VAL) GO TO 160
      XS = COUT(4)
      IND = 0
      GO TO 360
  160 IF (SIGN(1.D0,W(M,2)-VAL).EQ.SIGN(1.D0,W(M,4)-VAL)) GO TO 180
      X1 = COUT(4)
      X2 = COUT(5)
      OLDF = W(M,4)
      NEWF = W(M,2)
      GO TO 200
  180 IF (SIGN(1.D0,Y(M)-VAL).EQ.SIGN(1.D0,W(M,2)-VAL)) GO TO 280
      X1 = X
      X2 = COUT(4)
      OLDF = W(M,2)
      NEWF = Y(M)
  200 IND = -1
      I = 1
      IR = 2
      A = X2 - X1
      B = 2.D0*X1 - X2
      X1 = 1.D0
      X2 = 2.D0
      FX = NEWF - VAL
      D(1) = OLDF - VAL
      TOLX = 5.D0*COUT(11)
  220 CALL C05AZF(X1,X2,FX,TOLX,IR,D,IND,I)
      XS = A*X1 + B
      IF (IND.EQ.0) GO TO 240
      IF (IND.NE.4) GO TO 260
      IND1 = 1
      FX = 0.D0
      CALL D02XBF(XS,X,COUT,N,Y,W,N,M,FX,IND1)
      FX = FX - VAL
      IF (IND1.EQ.0) GO TO 220
C     IMPOSSIBLE EXIT FROM D02XBF
      IND = 7
      GO TO 360
  240 IF (I.EQ.0) GO TO 360
C     IMPOSSIBLE EXIT FROM C05AZF
  260 IND = 5
      GO TO 360
  280 IF (C(1).EQ.3.D0) GO TO 300
C     NO ROOT IN RANGE
      IND = 4
      GO TO 360
  300 J = 0
      DO 340 I = 1, N
         IF (I.EQ.M) GO TO 320
         IF ((Y(I).GE.W(I,9) .AND. W(I,2).GE.W(I,9) .AND. W(I,4)
     *       .GE.W(I,9)) .OR. (Y(I).LT.W(I,9) .AND. W(I,2).LT.W(I,9)
     *        .AND. W(I,4).LT.W(I,9))) GO TO 320
         W(I,9) = -2.D0*MAX(ABS(W(I,9)),ABS(Y(I)),ABS(W(I,2)),ABS(W(I,4)
     *            ),1.D0)
         IF (Y(I).LT.0.0D0) W(I,9) = -W(I,9)
         GO TO 340
  320    J = J + 1
  340 CONTINUE
      IF (J.EQ.N) GO TO 100
      COMM(3) = 1.D0
      GO TO 80
  360 IF (3.D0*COUT(3).GE.COUT(8)) TOL = -TOL
      X = XS
  380 IFAIL = P01ABF(IFAIL,IND,SRNAME,0,P01REC)
      RETURN
      END
