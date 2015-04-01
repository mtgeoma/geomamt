      SUBROUTINE F06HQF( PIVOT, DIRECT, N, ALPHA, X, INCX, C, S )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      COMPLEX*16         ALPHA
      INTEGER            INCX, N
      CHARACTER*1        DIRECT, PIVOT
C     .. Array Arguments ..
      COMPLEX*16         S( * ), X( * )
      DOUBLE PRECISION   C( * )
C     ..
C
C  F06HQF generates the parameters of a unitary matrix P such that
C
C     when   PIVOT = 'F' or 'f'   and   DIRECT = 'F' or 'f'
C     or     PIVOT = 'V' or 'v'   and   DIRECT = 'B' or 'b'
C
C        P*( alpha ) = ( beta ),
C          (   x   )   (   0  )
C
C     when   PIVOT = 'F' or 'f'   and   DIRECT = 'B' or 'b'
C     or     PIVOT = 'V' or 'v'   and   DIRECT = 'F' or 'f'
C
C        P*(   x   ) = (   0  ),
C          ( alpha ) = ( beta )
C
C  where alpha is a scalar and x is an n element vector.
C
C  When  PIVOT = 'F' or 'f'  ( fixed pivot )
C  and  DIRECT = 'F' or 'f'  ( forward sequence ) then
C
C     P is given as the sequence of plane rotation matrices
C
C        P = P( n )*P( n - 1 )*...*P( 1 )
C
C     where P( k ) is a plane rotation matrix for the ( 1, k + 1 ) plane
C     designed to annihilate the kth element of x.
C
C  When  PIVOT = 'V' or 'v'  ( variable pivot )
C  and  DIRECT = 'B' or 'b'  ( backward sequence ) then
C
C     P is given as the sequence of plane rotation matrices
C
C        P = P( 1 )*P( 2 )*...*P( n )
C
C     where P( k ) is a plane rotation matrix for the ( k, k + 1 ) plane
C     designed to annihilate the kth element of x.
C
C  When  PIVOT = 'F' or 'f'  ( fixed pivot )
C  and  DIRECT = 'B' or 'b'  ( backward sequence ) then
C
C     P is given as the sequence of plane rotation matrices
C
C        P = P( 1 )*P( 2 )*...*P( n )
C
C     where P( k ) is a plane rotation matrix for the ( k, n + 1 ) plane
C     designed to annihilate the kth element of x.
C
C  When  PIVOT = 'V' or 'v'  ( variable pivot )
C  and  DIRECT = 'F' or 'f'  ( forward sequence ) then
C
C     P is given as the sequence of plane rotation matrices
C
C        P = P( n )*P( n - 1 )*...*P( 1 )
C
C     where P( k ) is a plane rotation matrix for the ( k, k + 1 ) plane
C     designed to annihilate the kth element of x.
C
C  The routine returns the cosine, c( k ), and sine, s( k ) that define
C  the matrix P( k ), such that the two by two rotation part of P( k ),
C  R( k ), has the form
C
C     R( k ) = (  c( k )  conjg( s( k ) ) )  with  c( k ) real.
C              ( -s( k )         c( k )   )
C
C  On entry,  ALPHA must contain the scalar  alpha and on exit, ALPHA is
C  overwritten by  beta. If  alpha is real then  beta will also be real.
C
C  The vector  x  is overwritten by the  tangents of the plane rotations
C  ( t( k ) = s( k )/c( k ) ).
C
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 17-April-1985.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, IX
C     .. External Subroutines ..
      EXTERNAL           F06CAF
C     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
            IX = 1 + ( N - 1 )*INCX
            IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
               DO 10, I = N, 2, -1
                  CALL F06CAF( X( IX - INCX ), X( IX ), C( I ), S( I ) )
                  IX = IX - INCX
   10          CONTINUE
               CALL F06CAF( ALPHA, X( IX ), C( 1 ), S( 1 ) )
            ELSE IF( ( PIVOT.EQ.'F' ).OR.( PIVOT.EQ.'f' ) )THEN
C
C              Here we choose c and s so that
C
C                 ( alpha ) := (  c  conjg( s ) )*( alpha  )
C                 (   0   )    ( -s         c   ) ( x( i ) )
C
C              which is equivalent to
C
C                 (   0   ) := (        c    -s )*( x( i ) )
C                 ( alpha )    ( conjg( s )   c ) ( alpha  )
C
C              and so we need to return  s( i ) = -conjg( s )  in order
C              to make R( i ) look like
C
C                 R( i ) = (  c( i )  s( i ) ).
C                          ( -s( i )  c( i ) )
C
               DO 20, I = N, 1, -1
                  CALL F06CAF( ALPHA, X( IX ), C( I ), S( I ) )
                  S( I )  = -DCONJG( S( I ) )
                  X( IX ) = -DCONJG( X( IX ) )
                  IX      =  IX                - INCX
   20          CONTINUE
            END IF
         ELSE IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
            IX = 1
            IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
C
C              Here we choose c and s so that
C
C                 ( x( i + 1 ) ) := (  c  conjg( s ) )*( x( i + 1 ) )
C                 (    0       )    ( -s         c   ) ( x( i )     )
C
C              which is equivalent to
C
C                 (    0       ) := (        c    -s )*( x( i )     )
C                 ( x( i + 1 ) )    ( conjg( s )   c ) ( x( i + 1 ) )
C
C              and so we need to return  s( i ) = -conjg( s )  in order
C              to make R( i ) look like
C
C                 R( i ) = (  c( i )  conjg( s( i ) ) ).
C                          ( -s( i )         c( i )   )
C
               DO 30, I = 1, N - 1
                  CALL F06CAF( X( IX + INCX ), X( IX ), C( I ), S( I ) )
                  S( I )  = -DCONJG( S( I ) )
                  X( IX ) = -DCONJG( X( IX ) )
                  IX      =  IX                + INCX
   30          CONTINUE
               CALL F06CAF( ALPHA, X( IX ), C( N ), S( N ) )
               S( N )  = -DCONJG( S( N ) )
               X( IX ) = -DCONJG( X( IX ) )
            ELSE IF( ( PIVOT.EQ.'F' ).OR.( PIVOT.EQ.'f' ) )THEN
               DO 40, I = 1, N
                  CALL F06CAF( ALPHA, X( IX ), C( I ), S( I ) )
                  IX = IX + INCX
   40          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06HQF. ( CSROTG )
C
      END
