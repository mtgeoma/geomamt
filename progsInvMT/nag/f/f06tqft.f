      SUBROUTINE F06TQF( N, ALPHA, X, INCX, A, LDA, C, S )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      COMPLEX*16         ALPHA
      INTEGER            INCX, LDA, N
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), S( * ), X( * )
      DOUBLE PRECISION   C( * )
C     ..
C
C  F06TQF performs the factorization
C
C     (     U    ) = Q*( R ),
C     ( alpha*x' )     ( 0 )
C
C  where U and R are n by n upper triangular matrices, x is an n element
C  vector,  alpha  is  a  scalar  and  Q  is  an  ( n + 1 ) by ( n + 1 )
C  unitary matrix.  If the diagonal elements of  U  are supplied as real
C  then the diagonal elements of R will also be real.
C
C  U must be supplied in the n by n upper triangular part of the array A
C  and this is  overwritten by  R.  Q is formed as a  sequence of  plane
C  rotations  in  planes  ( 1, n + 1 ), ( 2, n + 1 ), ..., ( n, n + 1 ),
C  the  rotation, Q( j ), in the  ( j, n + 1 )th  plane being  chosen to
C  annihilate the jth element of  x. The cosine and sine that define the
C  jth rotation are returned in c( j ) and s( j )  respectively, and the
C  tangent, t( j ), is overwritten on x( j ). The two by two part of the
C  jth plane rotation, P( j ), is of the form
C
C     P( j ) = (  c( j )  conjg( s( j ) ) ),   with  c( j )  real.
C              ( -s( j )         c( j )   )
C
C  The unitary matrix Q is given by
C
C     Q = conjg( ( Q( n )*...*Q( 2 )*Q( 1 ) )').
C
C  If n is less than unity then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 18-April-1985.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      COMPLEX*16         TEMP, XJX
      INTEGER            I, J, JX
C     .. External Subroutines ..
      EXTERNAL           F06CAF
C     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
C     ..
C     .. Executable Statements ..
      IF( ( N.GT.0 ).AND.( ALPHA.NE.ZERO ) )THEN
         JX = 1
         DO 20 J = 1, N
            XJX = ALPHA*X( JX )
            DO 10 I = 1, J - 1
               TEMP = A( I, J )
               A( I, J ) = DCONJG( S( I ) )*XJX + C( I )*TEMP
               XJX = C( I )*XJX - S( I )*TEMP
   10       CONTINUE
            CALL F06CAF( A( J, J ), XJX, C( J ), S( J ) )
            X( JX ) = XJX
            JX = JX + INCX
   20    CONTINUE
      END IF
C
      RETURN
C
C     End of F06TQF. ( CUTUPD )
C
      END
