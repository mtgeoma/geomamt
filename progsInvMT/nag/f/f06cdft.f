      SUBROUTINE F06CDF( T, C, S )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-604 (MAR 1988).
C     .. Scalar Arguments ..
      COMPLEX*16         C, T
      DOUBLE PRECISION   S
C     ..
C
C  F06CDF returns values c and s such that
C
C     c = cos( theta ),   s = sin( theta )
C
C  for a given value of
C
C     t = tan( theta ).
C
C  s is always real and c and s are given by
C
C     c = ( sign( real( t ) )*abs( t ) )/( t*sqrt( 1 + abs( t )**2 ) ),
C
C     s = c*t.
C
C  If  abs( t ) .le. eps, where eps is the relative machine precision as
C  returned by routine X02AJF then c and s are returned as
C
C     c = ( sign( real( t ) )*abs( t ) )/t   and   s = 0.0.
C
C  If in addition  abs( real( t ) ).le.tol and  abs( aimag( t ) ).le.tol
C  where tol is the small positive value
C
C     tol = flmin/eps,
C
C  flmin being the values returned by routine  X02AMF, then  c and s are
C  returned as
C
C     c = t/tol,  s = 0.0.
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 3-January-1986.
C     Sven Hammarling, Nag Central Office.
C  -- Modified on 4-December-1987.
C     Sven Hammarling and Jeremy Du Croz, Nag Central office.
C        No longer sets s to zero when t is less than eps.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER        ( ONE  = 1.0D+0 )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      DOUBLE PRECISION   ABST, EPS, RRTEPS, RTEPS, TOL
      LOGICAL            FIRST
C     .. External Functions ..
      DOUBLE PRECISION   X02AJF, X02AMF
      EXTERNAL           X02AJF, X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT, DIMAG, DBLE
C     .. Save statement ..
      SAVE               FIRST, EPS, RTEPS, RRTEPS, TOL
C     .. Data statements ..
      DATA               FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( T.EQ.ZERO )THEN
         S = ZERO
         C = ONE
      ELSE
C
         IF( FIRST )THEN
            FIRST  = .FALSE.
            EPS    =  X02AJF( )
            RTEPS  =  SQRT( EPS )
            RRTEPS =  1/RTEPS
            TOL    =  X02AMF( )/EPS
         END IF
C
         IF( ( ABS( DBLE ( T ) ).LE.TOL ).AND.
     $       ( ABS( DIMAG( T ) ).LE.TOL )      )THEN
            S = ZERO
            C = T/TOL
         ELSE
            ABST = ABS( T )
            IF( ABST.LT.RTEPS )THEN
               S = SIGN( ABST, DBLE( T ) )
               C = S/T
            ELSE IF( ABST.GT.RRTEPS )THEN
               S = SIGN( ONE, DBLE( T ) )
               C = S/T
            ELSE
               S = SIGN( ABST, DBLE( T ) )/SQRT( 1 + ABST**2 )
               C = S/T
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06CDF. ( CCSGS )
C
      END
