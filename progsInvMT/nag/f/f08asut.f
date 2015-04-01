      DOUBLE PRECISION FUNCTION F08ASU(X,Y,Z)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     DOUBLE PRECISION                 DLAPY3
C     ENTRY                            DLAPY3(X,Y,Z)
C
C  Purpose
C  =======
C
C  DLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause
C  unnecessary overflow.
C
C  Arguments
C  =========
C
C  X       (input) DOUBLE PRECISION
C  Y       (input) DOUBLE PRECISION
C  Z       (input) DOUBLE PRECISION
C          X, Y and Z specify the values x, y and z.
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C     .. Parameters ..
      DOUBLE PRECISION                 ZERO
      PARAMETER                        (ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X, Y, Z
C     .. Local Scalars ..
      DOUBLE PRECISION                 W, XABS, YABS, ZABS
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, MAX, SQRT
C     .. Executable Statements ..
C
      XABS = ABS(X)
      YABS = ABS(Y)
      ZABS = ABS(Z)
      W = MAX(XABS,YABS,ZABS)
      IF (W.EQ.ZERO) THEN
         F08ASU = ZERO
      ELSE
         F08ASU = W*SQRT((XABS/W)**2+(YABS/W)**2+(ZABS/W)**2)
      END IF
      RETURN
C
C     End of F08ASU (DLAPY3)
C
      END
