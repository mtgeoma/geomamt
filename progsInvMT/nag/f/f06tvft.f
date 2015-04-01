      SUBROUTINE F06TVF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( * )
      DOUBLE PRECISION   S( * )
C     ..
C
C  F06TVF applies a  given sequence  of  plane rotations  to either  the
C  left, or the right, of the n by n upper triangular matrix U with real
C  diagonal elements.  U  is transformed to an  upper Hessenberg matrix.
C  The rotations are applied in planes k1 up to k2.
C
C  The upper Hessenberg matrix, H, is formed as
C
C     H = P*U             when   SIDE = 'L' or 'l',  (  Left-hand side )
C
C  where P is a unitary matrix of the form
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 )
C
C  and is formed as
C
C     H = U*conjg( P' )   when   SIDE = 'R' or 'r',  ( Right-hand side )
C
C  where P is a unitary matrix of the form
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  P( k ) being a plane rotation matrix for the  ( k, k + 1 ) plane. The
C  cosine and sine that define P( k ), k = k1, k1 + 1, ..., k2 - 1, must
C  be  supplied  in  c( k )  and  s( k )  respectively.  The  two by two
C  rotation part of P( k ), R( k ), is assumed to have the form
C
C     R( k ) = ( conjg( c( k ) )  s( k ) ),  s( k ) real.
C              (       -s( k )    c( k ) )
C
C  The matrix  U must be supplied in the n by n leading upper triangular
C  part of the array  A, and this is overwritten by the upper triangular
C  part of  H.  The imaginary parts of the diagonal elements must be set
C  to zero.
C
C  The  sub-diagonal elements of  H, h( k + 1, k ),  are returned in the
C  elements s( k ),  k = k1, k1 + 1, ..., k2 - 1. Note that the diagonal
C  elements of H, h( k1, k1 ), h( k1 + 1, k1 + 1 ), ..., h( k2, k2 ) are
C  not real.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 18-January-1985.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      COMPLEX*16         AIJ, CTEMP, CTEMPC, TEMP
      DOUBLE PRECISION   STEMP
      INTEGER            I, J
C     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MIN, DBLE
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Apply the plane rotations to columns n back to k1.
C
         DO 20 J = N, K1, -1
            IF( J.GE.K2 )THEN
               AIJ = A( K2, J )
            ELSE
C
C              Form  the  additional sub-diagonal element  h( j + 1, j )
C              and store it in  s( j ).
C
               AIJ = DCONJG( C( J ) )*DBLE( A( J, J ) )
               S( J ) = -S( J )*DBLE( A( J, J ) )
            END IF
            DO 10 I = MIN( K2, J ) - 1, K1, -1
               TEMP = A( I, J )
               A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
               AIJ = S( I )*AIJ + DCONJG( C( I ) )*TEMP
   10       CONTINUE
            A( K1, J ) = AIJ
   20    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        Apply the  plane rotations to columns  k1 up to ( k2 - 1 )  and
C        form    the    additional   sub-diagonal   elements,    storing
C        h( j + 1, j ) in s( j ).
C
         DO 40 J = K1, K2 - 1
            IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
               STEMP = S( J )
               CTEMP = C( J )
               CTEMPC = DCONJG( CTEMP )
               DO 30 I = 1, J
                  TEMP = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMPC*TEMP - STEMP*A( I, J )
                  A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
   30          CONTINUE
               S( J ) = STEMP*DBLE( A( J + 1, J + 1 ) )
               A( J + 1, J + 1 ) = CTEMPC*DBLE( A( J + 1, J + 1 ) )
            END IF
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F06TVF. ( CUTSRH )
C
      END
