      SUBROUTINE D05BAR(WKY1,WKY2,THRESH,NDIM,INDEX)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     ------------------------------------------------------------
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     This subroutine finds the index of the maximum value of
C     the relative error.
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     ------------------------------------------------------------
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  THRESH
      INTEGER           INDEX, NDIM
C     .. Array Arguments ..
      DOUBLE PRECISION  WKY1(NDIM), WKY2(2*NDIM)
C     .. Local Scalars ..
      DOUBLE PRECISION  DENMAX
      INTEGER           I, I2
C     .. External Functions ..
      INTEGER           IDAMAX
      EXTERNAL          IDAMAX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Executable Statements ..
C
      DO 20 I = 1, NDIM
         I2 = 2*I
         DENMAX = MAX(ABS(WKY1(I)),ABS(WKY2(I2)),ABS(THRESH))
         WKY1(I) = (WKY2(I2)-WKY1(I))/DENMAX
   20 CONTINUE
      INDEX = IDAMAX(NDIM,WKY1,1)
      RETURN
      END
