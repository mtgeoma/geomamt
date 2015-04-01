      SUBROUTINE F06TSF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), S( * )
      DOUBLE PRECISION   C( * )
C     ..
C
C  F06TSF restores an upper spiked matrix  H to upper triangular form by
C  applying a sequence of plane rotations, in planes  k1 up to k2,  from
C  either the left, or the right.
C
C  The matrix  H is assumed to have non-zero elements only in the spiked
C  positions, h( k2, k ) for a row spike and h( k + 1, k1 ) for a column
C  spike, k = k1, k1 + 1, ..., k2 - 1, and these must be supplied in the
C  elements  s( k ).  The matrix  H  is assumed  to have  real  diagonal
C  elements in each position, except where the spike joins the diagonal,
C  that is,  h( k2, k2 )  for a row spike and  h( k1, k1 )  for a column
C  spike.  The  triangular matrix  will be  returned with  real diagonal
C  elements.
C
C  When  SIDE = 'L' or 'l'  ( Left-hand side )
C
C     H  is  assumed  to have a  row spike  and is restored to the upper
C     triangular matrix  R as
C
C        R = P*H,
C
C     where P is a unitary matrix of the form
C
C        P = D( k2 )*P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C     P( k )  being a plane rotation matrix for the  ( k, k2 ) plane and
C     D( k )  being a  unitary diagonal matrix  with a  non-unit element
C     only in the kth position.
C
C  When  SIDE = 'R' or 'r'  ( Right-hand side )
C
C     H  is assumed to have a  column spike and is restored to the upper
C     triangular matrix R as
C
C        R = H*conjg( P' ),
C
C     where P is a unitary matrix of the form
C
C        P = D( k1 )*P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C     P( k )  being a plane rotation matrix for the  ( k1, k + 1 ) plane
C     and D( k ) being a unitary diagonal matrix with a non-zero only in
C     the kth position.
C
C  The two by two rotation part of  P( k ), Q( k ),  is of the form
C
C     Q( k ) = (  c( k )  conjg( s( k ) ) ),  with c( k ) real,
C              ( -s( k )         c( k )   )
C
C  and  c( k ) and s( k ) are returned in the kth elements of the arrays
C  C and S respectively. d( k ) is returned in s( k2 ).
C
C  The upper triangular part of the matrix  H  must be supplied in the n
C  by n  leading upper triangular part of  A, and this is overwritten by
C  the upper triangular matrix R.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
C     .. Local Scalars ..
      COMPLEX*16         AIJ, SPIKE, STEMP, STEMPC, TEMP
      DOUBLE PRECISION   CTEMP
      INTEGER            I, J
C     .. External Subroutines ..
      EXTERNAL           ZSCAL, F06CAF
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DIMAG, DCONJG, MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Restore H to upper triangular form by annihilating the elements
C        in  the  spike  of  H.  The  jth rotation  is  chosen  so  that
C
C           ( h( j, j ) ) := (  c  conjg( s ) )*( h( j , j ) ).
C           (     0     )    ( -s         c   ) ( h( k2, j ) )
C
C        Apply the rotations in columns k1 up to ( k2 - 1 ).
C
         DO 20 J = K1, K2 - 1
            SPIKE = S( J )
            DO 10 I = K1, J - 1
               AIJ = A( I, J )
               A( I, J ) = DCONJG( S( I ) )*SPIKE + C( I )*AIJ
               SPIKE = C( I )*SPIKE - S( I )*AIJ
   10       CONTINUE
C
C           Set up the rotation.
C
            CALL F06CAF( A( J, J ), SPIKE, C( J ), S( J ) )
   20    CONTINUE
C
C        Apply the rotations to column k2.
C
         TEMP = A( K2, K2 )
         DO 30 I = K1, K2 - 1
            AIJ = A( I, K2 )
            A( I, K2 ) = DCONJG( S( I ) )*TEMP + C( I )*AIJ
            TEMP = C( I )*TEMP - S( I )*AIJ
   30    CONTINUE
C
C        Form D( k2 ), chosen to make a( k2, k2 ) real.
C
         IF( DIMAG( TEMP ).NE.ZERO )THEN
            CTEMP = ABS( TEMP )
            S( K2 ) = DCONJG( TEMP )/CTEMP
            A( K2, K2 ) = CTEMP
         ELSE
            S( K2 ) = ONE
            A( K2, K2 ) = TEMP
         END IF
C
C        Apply the rotations to columns ( k2 + 1 ) up to n.
C
         DO 50 J = K2 + 1, N
            TEMP = A( K2, J )
            DO 40 I = K1, K2 - 1
               AIJ = A( I, J )
               A( I, J ) = DCONJG( S( I ) )*TEMP + C( I )*AIJ
               TEMP = C( I )*TEMP - S( I )*AIJ
   40       CONTINUE
C
C           Scale by D( k2 ).
C
            A( K2, J ) = S( K2 )*TEMP
   50    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        Restore  H  to upper triangular form  by annihilating the spike
C        of H. The jth rotation is chosen so that
C
C           ( h( j, j ) ) := (  c  conjg( s ) )*( h( j, j )  ),
C           (     0     )    ( -s         c   ) ( h( j, k1 ) )
C
C        which can be expressed as
C
C           ( 0  h( j, j ) ) :=
C
C               ( h( j, k1 )  h( j, j ) )*(  c  conjg( s ) ).
C                                         ( -s         c   )
C
C        Thus we return  c( j ) = c  and  s( j ) = -s  to make the plane
C        rotation matrix look like
C
C           Q( j ) = (  c( j )  conjg( s( j ) ) ).
C                    ( -s( j )         c( j )   )
C
         DO 80 J = K2, K1 + 1, -1
            CALL F06CAF( A( J, J ), S( J - 1 ), CTEMP, STEMP )
            STEMP = -STEMP
            S( J - 1 ) = STEMP
            C( J - 1 ) = CTEMP
            IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
               STEMPC = DCONJG( STEMP )
               DO 60 I = J - 1, K1 + 1, -1
                  SPIKE = S( I - 1 )
                  S( I - 1 ) = STEMP*A( I, J ) + CTEMP*SPIKE
                  A( I, J ) = CTEMP*A( I, J ) - STEMPC*SPIKE
   60          CONTINUE
               DO 70 I = K1, 1, -1
                  TEMP = A( I, K1 )
                  A( I, K1 ) = STEMP*A( I, J ) + CTEMP*TEMP
                  A( I, J ) = CTEMP*A( I, J ) - STEMPC*TEMP
   70          CONTINUE
            END IF
   80    CONTINUE
C
C        Form and apply D( k1 ).
C
         IF( DIMAG( A( K1, K1 ) ).NE.ZERO )THEN
            CTEMP = ABS( A( K1, K1 ) )
            S( K2 ) = A( K1, K1 )/CTEMP
            A( K1, K1 ) = CTEMP
            CALL ZSCAL( K1 - 1, DCONJG( S( K2 ) ), A( 1, K1 ), 1 )
         ELSE
            S( K2 ) = ONE
         END IF
      END IF
C
      RETURN
C
C     End of F06TSF. ( CUSQR )
C
      END
