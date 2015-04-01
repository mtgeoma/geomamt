      SUBROUTINE F06HTF( N, DELTA, Y, INCY, THETA, Z, INCZ )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      COMPLEX*16         DELTA, THETA
      INTEGER            INCY, INCZ, N
C     .. Array Arguments ..
      COMPLEX*16         Y( * ), Z( * )
C     ..
C
C  F06HTF performs a generalized Householder reflection given by
C
C     ( delta ) := P*( delta ),
C     (   y   )      (   y   )
C
C  where the unitary matrix P is given in the form
C
C     P = I - gamma*( zeta )*( zeta  conjg( z' ) ),
C                   (   z  )
C
C  z is an n element vector, gamma is a scalar such that
C
C     real( gamma ) = 1.0
C
C  and zeta is a real scalar that satisfies
C
C     1.0 .le. zeta .le. sqrt( 2.0 ).
C
C  gamma and zeta must be supplied in THETA such that
C
C     THETA = ( zeta, aimag( gamma ) ).
C
C  If  THETA = 0.0  then P is assumed to be the unit matrix and the
C  transformation is skipped, if  THETA .ne. 0.0   but
C  real( THETA ) .le. 0.0  then P is assumed to be
C
C     P = ( THETA  0 ).
C         (   0    I )
C
C  To perform the transformation
C
C     ( delta ) := conjg( P' )*( delta )
C     (   y   )                (   y   )
C
C  call F06HTF with DCONJG( THETA ) in place of THETA.
C
C  THETA and Z should normally be as supplied by routine F06HRF.
C  If  n = 0  then F06HTF must be called with  real( THETA ) .le. 0.0.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- written on 30-August-1984.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER        ( ONE  = 1.0D+0 )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      COMPLEX*16         TEMP
      DOUBLE PRECISION   ZETA
C     .. External Functions ..
      COMPLEX*16         ZDOTC
      EXTERNAL           ZDOTC
C     .. External Subroutines ..
      EXTERNAL           ZAXPY
C     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, DIMAG, DBLE
C     ..
C     .. Executable Statements ..
      IF( THETA.NE.ZERO )THEN
         IF( DBLE( THETA ).LE.DBLE( ZERO ) )THEN
            DELTA = THETA*DELTA
         ELSE
            ZETA  = DBLE( THETA )
            TEMP  = ZETA*DELTA    + ZDOTC( N, Z, INCZ, Y, INCY )
            IF( DIMAG( THETA ).NE.DBLE( ZERO ) )
     $         TEMP = DCMPLX( ONE, DIMAG( THETA ) )*TEMP
            CALL ZAXPY( N, -TEMP, Z, INCZ, Y, INCY )
            DELTA = DELTA         - TEMP*ZETA
         END IF
      END IF
C
      RETURN
C
C     End of F06HTF. ( CGRF )
C
      END
