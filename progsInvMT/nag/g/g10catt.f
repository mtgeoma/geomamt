      SUBROUTINE G10CAT(Y,N)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     3-point smooth by moving averages weighted 1/4,  1/2,  1/4
C     this is called hanning
C
C     .. Scalar Arguments ..
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  Y1, Y2, Y3
      INTEGER           I
C     .. Executable Statements ..
C
      Y2 = Y(1)
      Y3 = Y(2)
C
      DO 20 I = 2, N - 1
         Y1 = Y2
         Y2 = Y3
         Y3 = Y(I+1)
         Y(I) = (Y1+Y2+Y2+Y3)/4.0D0
   20 CONTINUE
      RETURN
      END
