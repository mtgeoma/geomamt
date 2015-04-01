      SUBROUTINE F06CCF( T, C, S )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-603 (MAR 1988).
C     MARK 15 REVISED. IER-943 (APR 1991).
C     .. Scalar Arguments ..
      COMPLEX*16         S, T
      DOUBLE PRECISION   C
C     ..
C
C  F06CCF returns values c and s such that
C
C     c = cos( theta ),   s = sin( theta )
C
C  for a given value of
C
C     t = tan( theta ).
C
C  c is always real and non-negative and c and s are given by
C
C     c = 1.0/sqrt( 1.0 + abs( t )**2 ),   s = c*t.
C
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 28-February-1986.
C     Sven Hammarling, Nag Central Office.
C  -- Modified 4-December-1987.
C     Sven Hammarling and Jeremy Du Croz, Nag Central Office.
C        No longer sets s to zero when t is less than eps.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
      DOUBLE PRECISION   ONE
      PARAMETER        ( ONE  = 1.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   ABST, EPS, FLMAX, RRTEPS, RTEPS
      LOGICAL            FIRST
C     .. External Functions ..
      DOUBLE PRECISION   X02AJF, X02AMF
      EXTERNAL           X02AJF, X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DIMAG, SQRT
C     .. Save statement ..
      SAVE               FIRST, EPS, FLMAX, RTEPS, RRTEPS
C     .. Data statements ..
      DATA               FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( FIRST )THEN
         FIRST  = .FALSE.
         EPS    =  X02AJF( )
         FLMAX  =  1/( 2*X02AMF( ) )
         RTEPS  =  SQRT  ( EPS )
         RRTEPS =  1/RTEPS
      END IF
C
      ABST = ABS( T )
      IF( ABST.GT.RRTEPS )THEN
         IF( ABST.LT.FLMAX )THEN
            C = 1/ABST
            S = T/ABST
         ELSE
            C = ZERO
            S = DCMPLX( DBLE( T )/ABST,DIMAG( T )/ABST )
         END IF
      ELSE IF( ABST.GE.RTEPS )THEN
         C = 1/SQRT( 1 + ABST**2 )
         S = C*T
      ELSE
         C = ONE
         S = T
      END IF
C
      RETURN
C
C     End of F06CCF. ( CCSG )
C
      END
