      SUBROUTINE F06CBF( A, B, C, S )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      COMPLEX*16         A, B, C
      DOUBLE PRECISION   S
C     ..
C
C  F06CBF returns c and s such that
C
C     (  conjg( c )  s )*( a ) = ( alpha ),
C     (        -s    c ) ( b ) = (   0   )
C
C  where
C
C     c = ( sign( real( t ) )*abs( t ) )/( t*sqrt( 1 + abs( t )**2 ) ),
C
C     s = ct,   t = b/a.
C
C  alpha is overwritten on a and t is overwritten on b. Note that when b
C  is real then alpha is also real.
C
C  When  abs( b ) .le. eps*abs( a ),  where  eps is the relative machine
C  precision  as returned  by routine  X02AJF,  but  aimag( a ) .ne. 0.0
C  then c and s are returned as
C
C     c = ( sign( real( a ) )*a )/abs( a )   and   s = 0.0.
C
C  If in addition  b .eq. 0.0  then t is returned as
C
C     t = tol*c,
C
C  where  tol  is the  small positive value given by
C
C     tol = flmin/eps,
C
C  flmin being the value returned by routine X02AMF.
C
C  When   abs( b ) .le. eps*abs( a ),  but  aimag( a ) .ne. 0.0  then  c
C  and s are returned as
C
C     c = 1.0   and   s = 0.0.
C
C  If in addition  b .eq. 0.0  then t is returned as
C
C     t = 0.0.
C
C  c and s  can be reconstructed from the  tangent, t,  by a call to the
C  Nag basic linear algebra routine F06CDF.
C
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 3-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      COMPLEX*16         T
      DOUBLE PRECISION   TEMP, TOL
      LOGICAL            FAIL, FIRST
C     .. External Functions ..
      COMPLEX*16         F06CLF
      DOUBLE PRECISION   X02AJF, X02AMF
      EXTERNAL           F06CLF, X02AJF, X02AMF
C     .. External Subroutines ..
      EXTERNAL           F06CDF
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, DCONJG, DIMAG, DBLE
C     .. Save statement ..
      SAVE               FIRST, TOL
C     .. Data statements ..
      DATA               FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( FIRST )THEN
         FIRST = .FALSE.
         TOL   =  X02AMF( )/X02AJF( )
      END IF
C
      IF( B.EQ.ZERO )THEN
         S = ZERO
         IF( DIMAG( A ).NE.DBLE( ZERO ) )THEN
            TEMP = SIGN( ABS( A ), DBLE( A ) )
            C    = A/TEMP
            A    = TEMP
            B    = TOL*C
         ELSE
            C    = ONE
         END IF
      ELSE
         T = F06CLF( B, A, FAIL )
         CALL F06CDF( T, C, S )
         IF( DIMAG( B ).EQ.DBLE( ZERO ) )THEN
            A = DBLE( DCONJG( C )*A ) + S*DBLE( B )
         ELSE
            A =       DCONJG( C )*A +   S*      B
         END IF
         B = T
      END IF
      RETURN
C
C     End of F06CBF. ( CROTGS )
C
      END
