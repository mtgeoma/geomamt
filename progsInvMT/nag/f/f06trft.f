      SUBROUTINE F06TRF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( * )
      DOUBLE PRECISION   S( * )
C     ..
C
C  F06TRF restores an upper Hessenberg matrix H to upper triangular form
C  by  applying a sequence of  plane rotations  from either the left, or
C  the right. The matrix H is assumed to have real non-zero sub-diagonal
C  elements  in  positions  h( k + 1, k ),  k = k1, k1 + 1, ..., k2 - 1,
C  only  and  h( k + 1, k )  must  be  supplied  in  s( k ).  The  upper
C  triangular  matrix  will  be  returned with  real  diagonal elements.
C
C  H is restored to the upper triangular matrix R either as
C
C     R = P*H,             when   SIDE = 'L' or 'l'  (  Left-hand side )
C
C  where P is a unitary matrix of the form
C
C     P = D( k2 )*P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  or as
C
C     R = H*conjg( P' ),   when   SIDE = 'R' or 'r'  ( Right-hand side )
C
C  where P is a unitary matrix of the form
C
C     P = D( k1 )*P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C  in both cases  P( k )  being a  plane rotation  for the  ( k, k + 1 )
C  plane and D( k ) is a unitary diagonal matrix with a non-unit element
C  only in the kth position.  The cosine and sine that define P( k ) are
C  returned in c( k ) and s( k ) respectively and  d( k ) is returned in
C  c( k2 ).  The two by two rotation part of  P( k ), Q( k ),  is of the
C  form
C
C     Q( k ) = ( conjg( c( k ) )  s( k ) ),  with s( k ) real.
C              (       -s( k )    c( k ) )
C
C  The upper triangular part of the matrix H must be supplied in the n
C  by n leading upper triangular part of A, and this is overwritten by
C  the upper triangular matrix R.
C
C  If n or k1 are less than unity, or k1 is not less than k2, or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 18-April-1985.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      COMPLEX*16         AIJ, CTEMP, CTEMPC, SUBH, TEMP
      DOUBLE PRECISION   STEMP
      INTEGER            I, J
C     .. External Subroutines ..
      EXTERNAL           ZSCAL, F06CBF
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DIMAG, DCONJG, MIN, DBLE
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Restore   H  to  upper  triangular  form  by  annihilating  the
C        sub-diagonal elements of  H. The jth rotation is chosen so that
C
C        ( h( j, j ) ) := (  conjg( c )  s )*( h( j, j )     ).
C        (     0     )    (        -s    c ) ( h( j + 1, j ) )
C
C        Apply the rotations in columns k1 up to n.
C
         DO 20 J = K1, N
            AIJ = A( K1, J )
            DO 10 I = K1, MIN( J, K2 ) - 1
               TEMP = A( I + 1, J )
               A( I, J ) = S( I )*TEMP + DCONJG( C( I ) )*AIJ
               AIJ = C( I )*TEMP - S( I )*AIJ
   10       CONTINUE
            IF( J.LT.K2 )THEN
C
C              Set up the rotation.
C
               SUBH = S( J )
               CALL F06CBF( AIJ, SUBH, C( J ), S( J ) )
               A( J, J ) = DBLE( AIJ )
            ELSE IF( J.EQ.K2 )THEN
C
C              Form D( k2 ), chosen to make a( k2, k2 ) real.
C
               A( K2, K2 ) = AIJ
               IF( DIMAG( A( K2, K2 ) ).NE.ZERO )THEN
                  STEMP = ABS( A( K2, K2 ) )
                  C( K2 ) = DCONJG( A( K2, K2 ) )/STEMP
                  A( K2, K2 ) = STEMP
               ELSE
                  C( K2 ) = ONE
               END IF
            ELSE
C
C              Scale by D( k2 ).
C
               A( K2, J ) = AIJ*C( K2 )
            END IF
   20    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        Restore   H  to  upper  triangular  form  by  annihilating  the
C        sub-diagonal elements of  H. The jth rotation is chosen so that
C
C           ( h( j + 1, j + 1 ) ) :=
C           (         0         )
C
C              ( conjg( c )  s )*( h( j + 1, j + 1 ) ),
C              (       -s    c ) ( h( j + 1, j )     )
C
C        which can be expressed as
C
C           ( 0  h( j + 1, j + 1 ) ) :=
C
C               ( h( j + 1, j )  h( j + 1, j + 1 ) )*(  c         s   ).
C                                                    ( -s  conjg( c ) )
C
C        Thus we return  c( j ) = c  and  s( j ) = -s  to make the plane
C        rotation matrix look like
C
C           Q( j ) = ( conjg( c( j ) )  s( j ) ).
C                    (       -s( j )    c( j ) )
C
         DO 40 J = K2 - 1, K1, -1
            SUBH = S( J )
            CALL F06CBF( A( J + 1, J + 1 ), SUBH, CTEMP, STEMP )
            STEMP = -STEMP
            S( J ) = STEMP
            C( J ) = CTEMP
            IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
               CTEMPC = DCONJG( CTEMP )
               DO 30 I = J, 1, -1
                  TEMP = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMPC*TEMP - STEMP*A( I, J )
                  A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
   30          CONTINUE
            END IF
   40    CONTINUE
C
C        Form and apply D( k1 ).
C
         IF( DIMAG( A( K1, K1 ) ).NE.ZERO )THEN
            STEMP = ABS( A( K1, K1 ) )
            C( K2 ) = A( K1, K1 )/STEMP
            A( K1, K1 ) = STEMP
            CALL ZSCAL( K1 - 1, DCONJG( C( K2 ) ), A( 1, K1 ), 1 )
         ELSE
            C( K2 ) = ONE
         END IF
      END IF
C
      RETURN
C
C     End of F06TRF. ( CUHQR )
C
      END
