      SUBROUTINE D02BHF(X,XEND,N,Y,TOL,IRELAB,HMAX,FCN,G,W,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     CALCULATES THE POINT X,WHERE G(X,Y)=0
C     FCN
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02BHF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HMAX, TOL, X, XEND
      INTEGER           IFAIL, IRELAB, N
C     .. Array Arguments ..
      DOUBLE PRECISION  W(N,7), Y(N)
C     .. Function Arguments ..
      DOUBLE PRECISION  G
      EXTERNAL          G
C     .. Subroutine Arguments ..
      EXTERNAL          FCN
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B, FX, NEWF, OLDF, TOLX, X1, X2, XS
      INTEGER           I, IND, IND1, IR
      LOGICAL           START
C     .. Local Arrays ..
      DOUBLE PRECISION  C(6), COMM(5), CON(3), COUT(14), D(17)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          C05AZF, D02PAF, D02XAF
C     .. Intrinsic Functions ..
      INTRINSIC         SIGN
C     .. Executable Statements ..
      IF (TOL.GT.0.D0 .AND. N.GT.0 .AND. IRELAB.GE.0 .AND. IRELAB.LE.2)
     *    GO TO 20
C     INPUT ERROR
      IND = 1
      GO TO 440
   20 OLDF = G(X,Y)
      IF (OLDF.NE.0.D0) GO TO 40
      IND = 0
      GO TO 440
   40 START = .TRUE.
      C(1) = 1.D0
      C(2) = 0.D0
      IF (IRELAB.EQ.1) C(2) = 1.D0
      IF (IRELAB.EQ.2) C(2) = 2.D0
      C(3) = 0.D0
      C(4) = HMAX
      C(5) = 0.D0
      COMM(1) = 0.D0
      COMM(2) = 0.D0
      COMM(3) = 0.D0
      COMM(4) = 1.D0
      CON(1) = 0.D0
      CON(2) = 0.D0
      CON(3) = 0.D0
   60 IND = 1
      CALL D02PAF(X,XEND,N,Y,C,TOL,FCN,COMM,CON,COUT,W,N,7,IND)
      XS = X
      IF (IND.NE.0) GO TO 80
      IF (C(1).EQ.2.D0 .OR. C(1).EQ.5.D0) GO TO 100
C     ERROR EXITS FROM D02PAF
   80 IF (IND.EQ.3) IND = 2
      IF (IND.EQ.4) IND = 3
      IF (IND.NE.2 .AND. IND.NE.3) IND = 6
      IF (IND.EQ.2) GO TO 420
      GO TO 440
  100 IF (C(1).EQ.2.D0 .AND. X.EQ.XEND) GO TO 120
      IF (C(1).EQ.5.D0 .AND. X.NE.XEND) GO TO 120
      IND = 6
      GO TO 420
  120 IF ( .NOT. START) GO TO 200
      NEWF = G(COUT(4),W(1,2))
      IF (NEWF.NE.0.D0) GO TO 160
      XS = COUT(4)
      DO 140 I = 1, N
         Y(I) = W(I,2)
  140 CONTINUE
      GO TO 420
  160 IF (SIGN(1.D0,NEWF).NE.SIGN(1.D0,OLDF)) GO TO 180
      START = .FALSE.
      OLDF = NEWF
      GO TO 200
  180 X1 = COUT(5)
      X2 = COUT(4)
      GO TO 260
  200 NEWF = G(X,Y)
      IF (NEWF.EQ.0.D0) GO TO 420
      IF (SIGN(1.D0,NEWF).NE.SIGN(1.D0,OLDF)) GO TO 240
      IF (C(1).EQ.5.D0) GO TO 220
C     NO ROOT IN RANGE
      IND = 4
      GO TO 420
  220 OLDF = NEWF
      GO TO 60
  240 X1 = COUT(4)
      X2 = X
  260 TOLX = 5.D0*COUT(11)
      IR = 2
      A = X2 - X1
      B = 2.D0*X1 - X2
      X1 = 1.D0
      X2 = 2.D0
      IND = -1
      FX = OLDF
      D(1) = NEWF
      I = 1
  280 CALL C05AZF(X1,X2,FX,TOLX,IR,D,IND,I)
      XS = X1*A + B
      IF (IND.EQ.0) GO TO 340
      IF (IND.NE.4) GO TO 360
      IND1 = 1
      CALL D02XAF(XS,X,COUT,N,Y,W,N,W(1,7),IND1)
      IF (IND1.EQ.0) GO TO 320
C     IMPOSSIBLE EXIT FROM D02XAF
  300 IND = 7
      GO TO 420
  320 FX = G(XS,W(1,7))
      GO TO 280
  340 IF (I.EQ.0) GO TO 380
C     IMPOSSIBLE EXIT FROM C05AZF
  360 IND = 5
      GO TO 420
  380 IND1 = 1
      CALL D02XAF(XS,X,COUT,N,Y,W,N,W(1,7),IND1)
      IF (IND1.NE.0) GO TO 300
      DO 400 I = 1, N
         Y(I) = W(I,7)
  400 CONTINUE
  420 IF (3.D0*COUT(3).GE.COUT(8)) TOL = -TOL
      X = XS
  440 IFAIL = P01ABF(IFAIL,IND,SRNAME,0,P01REC)
      RETURN
      END
