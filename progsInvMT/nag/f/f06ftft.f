      SUBROUTINE F06FTF( N, DELTA, Y, INCY, ZETA, Z, INCZ )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   DELTA, ZETA
      INTEGER            INCY, INCZ, N
C     .. Array Arguments ..
      DOUBLE PRECISION   Y( * ), Z( * )
C     ..
C
C  F06FTF performs a Householder reflection given by
C
C     ( delta ) := P*( delta ),
C     (   y   )      (   y   )
C
C  where the orthogonal matrix P is given in the form
C
C     P = I - ( zeta )*( zeta  z' ),
C             (   z  )
C
C  z being an n element vector and zeta a scalar.
C
C  If  zeta = 0.0  then P is assumed to be the unit matrix and the
C  transformation is skipped otherwise zeta must be in the range
C  ( 1.0, sqrt( 2.0 ) ).
C
C  zeta and z should normally be as supplied by routine F06FRF.
C  If n = 0 then F06FTF must be called with ZETA = 0.0.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 30-August-1984.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
C     .. External Functions ..
      DOUBLE PRECISION   DDOT
      EXTERNAL           DDOT
C     .. External Subroutines ..
      EXTERNAL           DAXPY
C     ..
C     .. Executable Statements ..
      IF( ZETA.NE.ZERO )THEN
         TEMP  = ZETA*DELTA + DDOT( N, Z, INCZ, Y, INCY )
         CALL DAXPY( N, -TEMP, Z, INCZ, Y, INCY )
         DELTA = DELTA      - TEMP*ZETA
      END IF
C
      RETURN
C
C     End of F06FTF. ( SGRF )
C
      END
