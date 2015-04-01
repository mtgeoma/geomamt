      SUBROUTINE G04BBZ(N,Y,IFAC,LEVELS,B,RSS)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     General sweep routine
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RSS
      INTEGER           LEVELS, N
C     .. Array Arguments ..
      DOUBLE PRECISION  B(LEVELS), Y(N)
      INTEGER           IFAC(N)
C     .. Local Scalars ..
      INTEGER           I
C     .. Executable Statements ..
      RSS = 0.0D0
      DO 20 I = 1, N
         Y(I) = Y(I) - B(IFAC(I))
         RSS = RSS + Y(I)**2
   20 CONTINUE
      RETURN
      END
