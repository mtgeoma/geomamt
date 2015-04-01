      SUBROUTINE F06FUF( N, Z, INCZ, Z1, ALPHA, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, Z1
      INTEGER            INCX, INCZ, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Z( * )
C     ..
C
C  F06FUF performs a Householder reflection given by
C
C     ( alpha ) = P*( alpha ) ,
C     (   x   )     (   x   )
C
C  where the orthogonal matrix p is given in the form
C
C     P = I - ( 1/z( 1 ) )*z*z'.
C
C  z( 1 ) must be supplied in Z1 and the remaining n elements in Z.
C  If Z1 is zero then P is assumed to be the unit matrix and the
C  transformation is skipped, otherwise Z1 must be in the range
C  ( 1.0, 2.0 ). Z1 and Z will usually be supplied by routine F06FSF.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 2-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   BETA
C     .. External Functions ..
      DOUBLE PRECISION   DDOT
      EXTERNAL           DDOT
C     .. External Subroutines ..
      EXTERNAL           DAXPY
C     ..
C     .. Executable Statements ..
      IF( Z1.NE.ZERO )THEN
         BETA  = ALPHA*Z1 + DDOT( N, Z, INCZ, X, INCX )
         ALPHA = ALPHA    - BETA
         CALL DAXPY( N, -BETA/Z1, Z, INCZ, X, INCX )
      END IF
C
      RETURN
C
C     End of F06FUF. ( SREF )
C
      END
