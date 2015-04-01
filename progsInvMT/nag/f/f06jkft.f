      DOUBLE PRECISION FUNCTION F06JKF( N, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      DOUBLE PRECISION          DZASUM
      ENTRY                     DZASUM( N, X, INCX )
C     .. Scalar Arguments ..
      INTEGER                           INCX, N
C     .. Array Arguments ..
      COMPLEX*16                        X( * )
C     ..
C
C  F06JKF returns the value, asum, given by
C
C     asum = abs( real( x( 1 ) ) ) + abs( aimag( x( 1 ) ) ) + ... +
C            abs( real( x( n ) ) ) + abs( aimag( x( n ) ) ),
C
C  for the vector x. asum is returned via the function name F06JKF.
C
C
C  Nag Fortran 77 version of the Blas routine DZASUM.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 12-February-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ZERO
      PARAMETER           ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      SUM
      INTEGER               IX
C     .. Intrinsic Functions ..
      INTRINSIC             ABS, DIMAG, DBLE
C     ..
C     .. Executable Statements ..
      SUM = ZERO
      IF( N.GT.0 )THEN
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            SUM = SUM + ABS( DBLE ( X( IX ) ) )
     $                + ABS( DIMAG( X( IX ) ) )
   10    CONTINUE
      END IF
C
      F06JKF = SUM
      RETURN
C
C     End of F06JKF. ( DZASUM )
C
      END
