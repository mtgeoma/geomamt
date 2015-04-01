      DOUBLE PRECISION FUNCTION F06BNF( A, B )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                  A, B
C     ..
C
C  F06BNF returns the value
C
C     p = sqrt( a*a + b*b )
C
C  via the function name.
C
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 17-January-1985.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ZERO
      PARAMETER           ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      P
C     .. Intrinsic Functions ..
      INTRINSIC             ABS, SQRT
C     ..
C     .. Executable Statements ..
      IF( A.EQ.ZERO )THEN
         P = ABS( B )
      ELSE IF( B.EQ.ZERO )THEN
         P = ABS( A )
      ELSE IF( ABS( A ).GE.ABS( B ) )THEN
         P = ABS( A )*SQRT( 1 + ( B/A )**2 )
      ELSE
         P = ABS( B )*SQRT( 1 + ( A/B )**2 )
      END IF
C
      F06BNF = P
      RETURN
C
C     End of F06BNF. ( SPYTH )
C
      END
