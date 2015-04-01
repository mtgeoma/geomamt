      SUBROUTINE A02ACF(XXR,XXI,YYR,YYI,ZR,ZI)
C     MARK 2A RELEASE.  NAG COPYRIGHT 1973
C     MARK 4.5 REVISED
C     MARK 5C REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     DIVIDES ONE COMPLEX NUMBER BY A SECOND
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  XXI, XXR, YYI, YYR, ZI, ZR
C     .. Local Scalars ..
      DOUBLE PRECISION  A, H, ONE
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Data statements ..
      DATA              ONE/1.0D0/
C     .. Executable Statements ..
C
      IF (ABS(YYR).LE.ABS(YYI)) GO TO 20
      H = YYI/YYR
      A = ONE/(H*YYI+YYR)
      ZR = (XXR+H*XXI)*A
      ZI = (XXI-H*XXR)*A
      RETURN
   20 H = YYR/YYI
      A = ONE/(H*YYR+YYI)
      ZR = (H*XXR+XXI)*A
      ZI = (H*XXI-XXR)*A
      RETURN
      END
