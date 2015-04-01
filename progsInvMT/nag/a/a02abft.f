      DOUBLE PRECISION FUNCTION A02ABF(XXR,XXI)
C     NAG COPYRIGHT 1975
C     MARK 4.5 REVISED
C     MARK 5C REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     RETURNS THE ABSOLUTE VALUE OF A COMPLEX NUMBER VIA ROUTINE
C     NAME
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 XXI, XXR
C     .. Local Scalars ..
      DOUBLE PRECISION                 H, ONE, XI, XR, ZERO
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, SQRT
C     .. Data statements ..
      DATA                             ZERO/0.0D0/, ONE/1.0D0/
C     .. Executable Statements ..
C
      XR = ABS(XXR)
      XI = ABS(XXI)
      IF (XI.LE.XR) GO TO 20
      H = XR
      XR = XI
      XI = H
   20 IF (XI.NE.ZERO) GO TO 40
      A02ABF = XR
      RETURN
   40 H = XR*SQRT(ONE+(XI/XR)**2)
      A02ABF = H
      RETURN
      END
