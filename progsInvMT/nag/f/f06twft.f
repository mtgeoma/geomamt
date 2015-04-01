      SUBROUTINE F06TWF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), S( * )
      DOUBLE PRECISION   C( * )
C     ..
C
C  F06TWF applies a  given sequence  of  plane rotations  to either  the
C  left, or the right, of the n by n upper triangular matrix U with real
C  diagonal elements.  U  is transformed to an upper spiked matrix.  The
C  rotations are applied in planes k1 up to k2.
C
C  The upper spiked matrix, H, is formed as
C
C     H = P*U             when   SIDE = 'L' or 'l',  (  Left-hand side )
C
C  where P is a unitary matrix of the form
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C  P( k ) being a plane rotation matrix for the ( k, k2 ) plane, and is
C  formed as
C
C     H = U*conjg( P' )   when   SIDE = 'R' or 'r',  ( Right-hand side )
C
C  where P is a unitary matrix of the form
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  P( k )  being a  plane rotation matrix for the  ( k1, k + 1 )  plane.
C
C  The cosine and sine that define  P( k ), k = k1, k1 + 1, ..., k2 - 1,
C  must be  supplied  in  c( k ) and s( k ) respectively. The two by two
C  rotation part of P( k ), R( k ), is assumed to have the form
C
C     R( k ) = (  c( k )  conjg( s( k ) ) ),  c( k ) real.
C              ( -s( k )         c( k )   )
C
C  The matrix  U must be supplied in the n by n leading upper triangular
C  part of the array  A, and this is overwritten by the upper triangular
C  part of  H.  The imaginary parts of the diagonal elements must be set
C  to zero.
C
C  When  SIDE = 'L' or 'l'  then a  row spike  is  generated  in  H  and
C  when  SIDE = 'R' or 'r'  then a  column spike is generated. For a row
C  spike the elements  h( k2, k )  and for a  column spike  the elements
C  h( k + 1, k1 ), k = k1, k1 + 1, ..., k2 - 1, are returned in  s( k ).
C  The  diagonal elements of  H  are all real on return,  except for the
C  position where the spike joins the diagonal, that is, h( k2, k2 ) for
C  a row spike and h( k1, k1 ) for a column spike.
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
C     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MIN, DBLE
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Apply the plane rotations to columns n back to k2.
C
         DO 20 J = N, K2, -1
            TEMP = A( K2, J )
            DO 10 I = K2 - 1, K1, -1
               AIJ = A( I, J )
               A( I, J ) = DCONJG( S( I ) )*TEMP + C( I )*AIJ
               TEMP = C( I )*TEMP - S( I )*AIJ
   10       CONTINUE
            A( K2, J ) = TEMP
   20    CONTINUE
C
C        Form  the spike  and apply the rotations in columns  ( k2 - 1 )
C        back to k1.
C
         DO 40 J = K2 - 1, K1, -1
            SPIKE = -S( J )*DBLE( A( J, J ) )
            A( J, J ) = C( J )*DBLE( A( J, J ) )
            DO 30 I = J - 1, K1, -1
               AIJ = A( I, J )
               A( I, J ) = DCONJG( S( I ) )*SPIKE + C( I )*AIJ
               SPIKE = C( I )*SPIKE - S( I )*AIJ
   30       CONTINUE
            S( J ) = SPIKE
   40    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        Apply  the plane rotations to columns  ( k1 + 1 ) up to k2  and
C        form the spike.
C
         DO 70 J = K1 + 1, K2
            CTEMP = C( J - 1 )
            STEMP = S( J - 1 )
            IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
               STEMPC = DCONJG( STEMP )
               DO 50 I = 1, K1
                  TEMP = A( I, K1 )
                  A( I, K1 ) = STEMP*A( I, J ) + CTEMP*TEMP
                  A( I, J ) = CTEMP*A( I, J ) - STEMPC*TEMP
   50          CONTINUE
               DO 60 I = K1 + 1, J - 1
                  SPIKE = S( I - 1 )
                  S( I - 1 ) = STEMP*A( I, J ) + CTEMP*SPIKE
                  A( I, J ) = CTEMP*A( I, J ) - STEMPC*SPIKE
   60          CONTINUE
               S( J - 1 ) = STEMP*DBLE( A( J, J ) )
               A( J, J ) = CTEMP*DBLE( A( J, J ) )
            END IF
   70    CONTINUE
      END IF
C
      RETURN
C
C     End of F06TWF. ( CUTSRS )
C
      END
