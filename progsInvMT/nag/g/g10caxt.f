      SUBROUTINE G10CAX(Y,N,ENDSAV)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     smooth by running medians (means) of 2
C     used to recenter results of running medians of 4
C     ENDSAV holds the original Y(N)
C
C     .. Parameters ..
      DOUBLE PRECISION  TWO
      PARAMETER         (TWO=2.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ENDSAV
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  Y(N)
C     .. Local Scalars ..
      INTEGER           I, NM1
C     .. Executable Statements ..
C
      NM1 = N - 1
      DO 20 I = 2, NM1
         Y(I) = (Y(I+1)+Y(I))/TWO
   20 CONTINUE
      Y(N) = ENDSAV
      RETURN
      END
