      DOUBLE PRECISION FUNCTION F06EKF( N, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      DOUBLE PRECISION          DASUM
      ENTRY                     DASUM ( N, X, INCX )
C     .. Scalar Arguments ..
      INTEGER                           INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION                  X( * )
C     ..
C
C  F06EKF returns the rectangular length, l1, of the vector x given by
C
C     l1 = abs( x( 1 ) ) + abs( x( 2 ) ) + ... + abs( x( n ) ).
C
C  l1 is returned via the function name F06EKF.
C
C
C  Nag Fortran 77 version of the Blas routine DASUM.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 27-November-1984.
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
      INTRINSIC             ABS
C     ..
C     .. Executable Statements ..
      SUM = ZERO
      IF( N.GT.0 )THEN
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            SUM = SUM + ABS( X( IX ) )
   10    CONTINUE
      END IF
C
      F06EKF = SUM
      RETURN
C
C     End of F06EKF. ( DASUM )
C
      END
