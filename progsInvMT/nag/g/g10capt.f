      SUBROUTINE G10CAP(X1,X2,X3,XMED,CHANGE)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     PUT THE MEDIAN OF X1, X2, X3 IN XMED AND
C     SET  CHANGE  .TRUE. IF THE MEDIAN ISNT X2.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X1, X2, X3, XMED
      LOGICAL           CHANGE
C     .. Local Scalars ..
      DOUBLE PRECISION  Y1, Y2, Y3
C     .. Executable Statements ..
C
      Y1 = X1
      Y2 = X2
      Y3 = X3
      XMED = Y2
      IF ((Y2-Y1)*(Y3-Y2).GE.0.0D0) GO TO 20
      CHANGE = .TRUE.
      XMED = Y1
      IF ((Y3-Y1)*(Y3-Y2).GT.0.0D0) GO TO 20
      XMED = Y3
   20 RETURN
      END
