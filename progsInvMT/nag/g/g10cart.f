      SUBROUTINE G10CAR(Y,N)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     estimate smoothed values for both end points
C     using the end point extrapolation rule
C     all the values except the end points have been smoothed
C
C     .. Scalar Arguments ..
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  Y0, YMED
      LOGICAL           CHANGE
C     .. External Subroutines ..
      EXTERNAL          G10CAP
C     .. Executable Statements ..
      CHANGE = .FALSE.
C
C     left end
C
      Y0 = 3.0D0*Y(2) - 2.0D0*Y(3)
      CALL G10CAP(Y0,Y(1),Y(2),YMED,CHANGE)
      Y(1) = YMED
C
C     right end
C
      Y0 = 3.0D0*Y(N-1) - 2.0D0*Y(N-2)
      CALL G10CAP(Y0,Y(N),Y(N-1),YMED,CHANGE)
      Y(N) = YMED
      RETURN
      END
