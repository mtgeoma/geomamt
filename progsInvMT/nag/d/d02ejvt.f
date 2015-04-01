      SUBROUTINE D02EJV(N,X,Y,HU,XLAST,H,RNQU,YSAV,ACOR,CHK,IMON,NSTPS,
     *                  GLAST,G,IFLAG,D,IFIN,DIR,ROOT,OUTPUT,PATH,XOUT)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DIR, GLAST, H, HU, RNQU, ROOT, X, XLAST, XOUT
      INTEGER           IFIN, IFLAG, IMON, N, NSTPS, PATH
C     .. Array Arguments ..
      DOUBLE PRECISION  ACOR(N), CHK(N), D(17), Y(N), YSAV(N,6)
C     .. Function Arguments ..
      DOUBLE PRECISION  G
      EXTERNAL          G
C     .. Subroutine Arguments ..
      EXTERNAL          OUTPUT
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B, FX, GNEW, TOLX, X1, X2, XNEW, XPREV, XS
      INTEGER           I, IND, IR, J, NQU
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          C05AZF, D02EJU, D02XKF
C     .. Intrinsic Functions ..
      INTRINSIC         SIGN, INT
C     .. Executable Statements ..
      IFIN = 0
      IF (IMON.NE.1) RETURN
      NSTPS = NSTPS + 1
      NQU = INT(RNQU)
      XPREV = XLAST
      IF (PATH.EQ.2) THEN
         CALL D02EJU(XOUT,OUTPUT,CHK,N,YSAV,ACOR,X,NQU,HU,H,IFLAG,XPREV,
     *               X,DIR)
         IF (IFLAG.NE.0) GO TO 100
         RETURN
      END IF
      IFLAG = 1
      CALL D02XKF(XLAST,CHK,N,YSAV,N,6,ACOR,N,X,NQU,HU,H,IFLAG)
      IF (IFLAG.NE.0) GO TO 100
      GNEW = G(XLAST,CHK)
      IF (GLAST*SIGN(1.0D0,GNEW).LE.0.0D0) THEN
         IFIN = 1
         ROOT = XLAST
         DO 20 I = 1, N
            Y(I) = CHK(I)
   20    CONTINUE
         IFLAG = 0
         RETURN
      END IF
      GLAST = GNEW
      GNEW = G(X,Y)
      IF (GLAST*SIGN(1.0D0,GNEW).GT.0.0D0) THEN
         XLAST = X
         GLAST = GNEW
         IF (PATH.EQ.4) THEN
            CALL D02EJU(XOUT,OUTPUT,CHK,N,YSAV,ACOR,X,NQU,HU,H,IFLAG,
     *                  XPREV,X,DIR)
            IF (IFLAG.NE.0) GO TO 100
         END IF
         RETURN
      END IF
      XNEW = X
   40 CONTINUE
      X1 = XLAST
      X2 = XNEW
      TOLX = 5.0D0*X02AJF()
      IR = 2
      A = X2 - X1
      B = 2.0D0*X1 - X2
      X1 = 1.0D0
      X2 = 2.0D0
      IND = -1
      FX = GLAST
      D(1) = GNEW
      I = 1
   60 CALL C05AZF(X1,X2,FX,TOLX,IR,D,IND,I)
      XS = X1*A + B
      IF (IND.EQ.4) THEN
         IFLAG = 1
         CALL D02XKF(XS,CHK,N,YSAV,N,6,ACOR,N,X,NQU,HU,H,IFLAG)
         IF (IFLAG.NE.0) GO TO 100
         FX = G(XS,CHK)
         GO TO 60
      END IF
      IMON = -2
      IFLAG = 7
      IF (IND.NE.0) RETURN
      IF (I.NE.0) RETURN
      IFLAG = 1
      CALL D02XKF(XS,CHK,N,YSAV,N,6,ACOR,N,X,NQU,HU,H,IFLAG)
      IF (IFLAG.NE.0) GO TO 100
      DO 80 J = 1, N
         Y(J) = CHK(J)
   80 CONTINUE
      ROOT = XS
      IFLAG = 0
      IFIN = 1
      IF (PATH.EQ.4) THEN
         CALL D02EJU(XOUT,OUTPUT,CHK,N,YSAV,ACOR,X,NQU,HU,H,IFLAG,XPREV,
     *               ROOT,DIR)
         IF (IFLAG.NE.0) IFLAG = 10
      END IF
      RETURN
  100 IMON = -2
      IF (IFLAG.NE.5) IFLAG = 8
      RETURN
      END
