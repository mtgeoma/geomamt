      SUBROUTINE A02AAF(XXR,XXI,YR,YI)
C     MARK 2A RELEASE.  NAG COPYRIGHT 1973
C     MARK 4.5 REVISED
C     MARK 5C REVISED
C     MARK 11C REVISED. IER-467 (MAR 1985)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     COMPUTES THE SQUARE ROOT OF A COMPLEX NUMBER
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  XXI, XXR, YI, YR
C     .. Local Scalars ..
      DOUBLE PRECISION  H, HALF, ONE, XI, XR, ZERO
C     .. External Functions ..
      DOUBLE PRECISION  A02ABF
      EXTERNAL          A02ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Data statements ..
      DATA              ZERO/0.0D0/, HALF/0.5D0/, ONE/1.0D0/
C     .. Executable Statements ..
C
      XR = ABS(XXR)
      XI = XXI
      IF (XR.GT.ONE) H = SQRT(XR*HALF+A02ABF(XR*HALF,XI*HALF))
      IF (XR.LE.ONE) H = SQRT(XR+A02ABF(XR,XI))*SQRT(HALF)
      IF (XI.NE.ZERO) XI = XI/(H+H)
      IF (XXR.LT.ZERO) GO TO 20
      YR = H
      YI = XI
      RETURN
   20 IF (XI.LT.ZERO) GO TO 40
      YR = XI
      YI = H
      RETURN
   40 YR = -XI
      YI = -H
      RETURN
      END
