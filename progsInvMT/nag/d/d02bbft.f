      SUBROUTINE D02BBF(X,XEND,N,Y,TOL,IRELAB,FCN,OUTPUT,W,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     FCN, OUTPUT
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02BBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL, X, XEND
      INTEGER           IFAIL, IRELAB, N
C     .. Array Arguments ..
      DOUBLE PRECISION  W(N,7), Y(N)
C     .. Subroutine Arguments ..
      EXTERNAL          FCN, OUTPUT
C     .. Local Scalars ..
      DOUBLE PRECISION  XS, XSOL
      INTEGER           IND
      LOGICAL           MARK, MARK1
C     .. Local Arrays ..
      DOUBLE PRECISION  CIN(6), COMM(5), CONST(3), COUT(14)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02PAF, D02XAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SIGN
C     .. Executable Statements ..
      IF (N.GT.0 .AND. TOL.GT.0.D0 .AND. IRELAB.GE.0 .AND. IRELAB.LE.2)
     *    GO TO 20
C     INPUT ERROR
      IND = 1
      GO TO 240
   20 XSOL = X
      XS = X
      IND = 0
      CALL OUTPUT(XSOL,Y)
      IF (X.NE.XEND) GO TO 40
      IF (XSOL.EQ.XEND) GO TO 240
C     X.EQ.XEND.NE.XSOL
      IND = 4
      GO TO 240
   40 MARK = .FALSE.
      MARK1 = .FALSE.
      CIN(1) = 1.D0
      CIN(2) = 0.D0
      IF (IRELAB.EQ.1) CIN(2) = 1.D0
      IF (IRELAB.EQ.2) CIN(2) = 2.D0
      CIN(3) = 0.D0
      CIN(4) = 0.D0
      CIN(5) = 0.D0
      COMM(1) = 0.D0
      COMM(2) = 0.D0
      COMM(3) = 0.D0
      COMM(4) = 0.D0
      CONST(1) = 0.D0
      CONST(2) = 0.D0
      CONST(3) = 0.D0
   60 IF (XSOL.EQ.XS) GO TO 200
      IF (XSOL.EQ.XEND .AND. X.EQ.XEND) GO TO 160
      IF (XSOL.EQ.XEND) GO TO 100
      IF (SIGN(1.D0,XSOL-XS).NE.SIGN(1.D0,XSOL-XEND)) GO TO 80
      IF (ABS(XSOL-X).LT.ABS(XSOL-XEND)) GO TO 200
      MARK = .TRUE.
      GO TO 100
   80 IF (XSOL.EQ.X) GO TO 160
      IF (SIGN(1.D0,XSOL-X).NE.SIGN(1.D0,XSOL-XS)) GO TO 160
      COMM(4) = -1.D0
      COMM(5) = XSOL
  100 IND = 1
      MARK1 = .TRUE.
      CALL D02PAF(X,XEND,N,Y,CIN,TOL,FCN,COMM,CONST,COUT,W,N,7,IND)
      IF (IND.EQ.0) GO TO 120
      IF (IND.EQ.1 .OR. IND.EQ.5 .OR. IND.EQ.7) IND = 6
      IF (IND.EQ.3) IND = 2
      IF (IND.EQ.4) IND = 3
      IF (IND.EQ.2) GO TO 220
      GO TO 240
  120 IF (CIN(1).EQ.2.D0 .AND. X.EQ.XEND) GO TO 140
      IF (CIN(1).EQ.6.D0 .AND. X.NE.XEND) GO TO 140
      IND = 6
      GO TO 220
  140 IF (MARK) GO TO 220
  160 IND = 1
      CALL D02XAF(XSOL,X,COUT,N,Y,W,N,W(1,7),IND)
      IF (IND.EQ.0) GO TO 180
      IND = 7
      GO TO 220
  180 XS = XSOL
      CALL OUTPUT(XSOL,W(1,7))
      IF (XS.EQ.XEND .AND. X.EQ.XEND) GO TO 220
      GO TO 60
C     XSOL OUTSIDE RANGE
  200 IND = 5
  220 IF ( .NOT. MARK1) GO TO 240
      IF (3.D0*COUT(3).GE.COUT(8)) TOL = -TOL
  240 IFAIL = P01ABF(IFAIL,IND,SRNAME,0,P01REC)
      RETURN
      END
