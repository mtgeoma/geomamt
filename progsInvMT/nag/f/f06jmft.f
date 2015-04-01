      INTEGER FUNCTION F06JMF( N, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      INTEGER          IZAMAX
      ENTRY            IZAMAX( N, X, INCX )
C     .. Scalar Arguments ..
      INTEGER                  INCX, N
C     .. Array Arguments ..
      COMPLEX*16               X( * )
C     ..
C
C  F06JMF returns the smallest value of i such that
C
C     alpha( i ) = max( abs( real( x( j ) ) ) + abs( imag( x( j ) ) ) )
C                   j
C
C  Nag Fortran 77 version of the Blas routine IZAMAX.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 31-May-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      DOUBLE PRECISION         TEMP, XMAX
      INTEGER                  I, IMAX, IX
C     .. Intrinsic Functions ..
      INTRINSIC                ABS, DIMAG, DBLE
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IMAX = 1
         IF( N.GT.1 )THEN
            XMAX = ABS( DBLE( X( 1 ) ) ) + ABS( DIMAG( X( 1 ) ) )
            IX   = 1
            DO 10, I = 2, N
               IX   = IX                     + INCX
               TEMP = ABS( DBLE( X( IX ) ) ) + ABS( DIMAG( X( IX ) ) )
               IF( XMAX.LT.TEMP )THEN
                  XMAX = TEMP
                  IMAX = I
               END IF
   10       CONTINUE
         END IF
      ELSE
         IMAX = 0
      END IF
C
      F06JMF = IMAX
      RETURN
C
C     End of F06JMF. ( IZAMAX )
C
      END
