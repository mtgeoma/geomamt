      SUBROUTINE F06FLF( N, X, INCX, XMAX, XMIN )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   XMAX, XMIN
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06FLF returns the values xmax and xmin given by
C
C     xmax = max( abs( x( i ) ) ),   xmin = min( abs( x( i ) ) ).
C             i                              i
C
C  If n is less than unity then xmax and xmin are returned as zero.
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 27-February-1986.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      INTEGER            IX
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
C     ..
C     .. Executable Statements ..
      IF( N.LT.1 )THEN
         XMAX = ZERO
         XMIN = ZERO
      ELSE
         XMAX = ABS( X( 1 ) )
         XMIN = XMAX
         DO 10 IX = 1 + INCX, 1 + ( N - 1 )*INCX, INCX
            XMAX = MAX( XMAX, ABS( X( IX ) ) )
            XMIN = MIN( XMIN, ABS( X( IX ) ) )
   10    CONTINUE
      END IF
C
      RETURN
C
C     End of F06FLF. ( SCOND )
C
      END
