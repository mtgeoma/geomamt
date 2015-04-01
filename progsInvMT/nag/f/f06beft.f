      SUBROUTINE F06BEF( JOB, X, Y, Z, C, S )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   C, S, X, Y, Z
      CHARACTER*1        JOB
C     ..
C
C  F06BEF  returns  details  of  the  Jacobi  plane  rotation  given  by
C
C     ( a  0 ) = (  c  s )*( x  y )*( c  -s ).
C     ( 0  b )   ( -s  c ) ( y  z ) ( s   c )
C
C  If  JOB = 'B' or 'b'  then the rotation is chosen so that
C
C     c .ge. 1/sqrt( 2 ).
C
C  If  JOB = 'S' or 's'  then the rotation is chosen so that
C
C     0 .le. c .le. 1/sqrt( 2 ).
C
C  If  JOB = 'M' or 'm'  then the rotation is chosen so that
C
C     abs( a ) .ge. abs( b ).
C
C  a is overwritten on X,  b is overwritten on Z and Y is overwritten by
C  the tangent of the angle defining the rotation t = s/c.
C
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 27-May-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      DOUBLE PRECISION   ABST2, FLMIN, RRTEPS, RTEPS, T, TEMP, T2
      LOGICAL            FAIL, FIRST
C     .. External Functions ..
      DOUBLE PRECISION   F06BLF, X02AJF, X02AMF
      EXTERNAL           F06BLF, X02AJF, X02AMF
C     .. External Subroutines ..
      EXTERNAL           F06BCF
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
C     .. Save statement ..
      SAVE               FIRST, FLMIN, RTEPS, RRTEPS
C     .. Data statements ..
      DATA               FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( FIRST )THEN
         FIRST  = .FALSE.
         FLMIN  =  X02AMF( )
         RTEPS  =  SQRT( X02AJF( ) )
         RRTEPS =  1/RTEPS
      END IF
C
      T2    = F06BLF( 2*Y, X - Z, FAIL )
      ABST2 = ABS( T2 )
      IF( ABST2.LT.RTEPS )THEN
         T = T2/2
      ELSE IF( ABST2.GT.RRTEPS )THEN
         T = T2/( 1 + ABST2 )
      ELSE
         T = T2/( 1 + SQRT( 1 + T2**2 ) )
      END IF
      TEMP = T*Y
      X    = X + TEMP
      Z    = Z - TEMP
      IF( ( ( JOB.EQ.'S' ).OR.( JOB.EQ.'s' )          ).OR.
     $    ( ( ( JOB.EQ.'M' ).OR.( JOB.EQ.'m' ) ).AND.
     $      ( ABS( X ).LT.ABS( Z )             )      )     )THEN
         IF( ABS( T ).LT.FLMIN )THEN
            T =  SIGN( 1/FLMIN, T2 )
         ELSE
            T = -1/T
         END IF
         TEMP = X
         X    = Z
         Z    = TEMP
      END IF
      CALL F06BCF( T, C, S )
      Y = T
C
      RETURN
C
C     End of F06BEF. ( SROTJ )
C
      END
