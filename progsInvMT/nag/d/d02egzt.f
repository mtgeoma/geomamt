      SUBROUTINE D02EGZ(N,X,Y,HU,XLAST,H,RNQU,YSAV,ACOR,CHK,IMON,NSTPS,
     *                  GLAST,USRHMX,HMAX,M,VAL,IFLAG,D,IFIN,DIR,ROOT)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DIR, GLAST, H, HMAX, HU, RNQU, ROOT, VAL, X,
     *                  XLAST
      INTEGER           IFIN, IFLAG, IMON, M, N, NSTPS
      LOGICAL           USRHMX
C     .. Array Arguments ..
      DOUBLE PRECISION  ACOR(N), CHK(N), D(17), Y(N), YSAV(N,6)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B, FX, GNEW, TOLX, X1, X2, XNEW, XS
      INTEGER           I, IND, IR, J, NQU
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          C05AZF, D02XKF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SIGN, INT
C     .. Executable Statements ..
      IFIN = 0
      IF (IMON.NE.1) RETURN
      NSTPS = NSTPS + 1
      NQU = INT(RNQU)
      IFLAG = 1
      CALL D02XKF(XLAST,CHK,N,YSAV,N,6,ACOR,N,X,NQU,HU,H,IFLAG)
      IF (IFLAG.NE.0) GO TO 120
      GNEW = CHK(M) - VAL
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
      J = 0
      IF (USRHMX) J = INT(ABS(X-XLAST)/HMAX)
      DO 40 I = 1, J
         XNEW = XLAST + HMAX*DIR
         IFLAG = 1
         CALL D02XKF(XNEW,CHK,N,YSAV,N,6,ACOR,N,X,NQU,HU,H,IFLAG)
         IF (IFLAG.NE.0) GO TO 120
         GNEW = CHK(M) - VAL
         IF (GLAST*SIGN(1.0D0,GNEW).LE.0.0D0) GO TO 60
         GLAST = GNEW
         XLAST = XNEW
   40 CONTINUE
      GNEW = Y(M) - VAL
      IF (GLAST*SIGN(1.0D0,GNEW).GT.0.0D0) THEN
         XLAST = X
         GLAST = GNEW
         RETURN
      END IF
      XNEW = X
   60 CONTINUE
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
   80 CALL C05AZF(X1,X2,FX,TOLX,IR,D,IND,I)
      XS = X1*A + B
      IF (IND.EQ.4) THEN
         IFLAG = 1
         CALL D02XKF(XS,CHK,N,YSAV,N,6,ACOR,N,X,NQU,HU,H,IFLAG)
         IF (IFLAG.NE.0) GO TO 120
         FX = CHK(M) - VAL
         GO TO 80
      END IF
      IMON = -2
      IFLAG = 5
      IF (IND.NE.0) RETURN
      IF (I.NE.0) RETURN
      IFLAG = 1
      CALL D02XKF(XS,CHK,N,YSAV,N,6,ACOR,N,X,NQU,HU,H,IFLAG)
      IF (IFLAG.NE.0) GO TO 120
      DO 100 J = 1, N
         Y(J) = CHK(J)
  100 CONTINUE
      ROOT = XS
      IFLAG = 0
      IFIN = 1
      RETURN
  120 IMON = -2
      IFLAG = 7
      RETURN
      END
