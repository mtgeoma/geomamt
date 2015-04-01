      SUBROUTINE F06QPF( N, ALPHA, X, INCX, Y, INCY, A, LDA, C, S )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, N
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * ), X( * ), Y( * )
C     ..
C
C  F06QPF performs the factorization
C
C     alpha*x*y' + U = Q*R,
C
C  where  U and R  are  n by n upper triangular matrices,  x and y are n
C  element vectors,  alpha is a scalar and  Q  is an  n by n  orthogonal
C  matrix.
C
C  U must be supplied in the n by n upper triangular part of the array A
C  and this is overwritten by R.  x and y must be supplied in the arrays
C  X and Y and Y is unaltered on return.
C
C  Q is formed as two sequences of plane rotations,  P1 and P2.  P1 is a
C  sequence in planes  ( n - 1, n ), ( n - 2, n ), ..., ( 1, n ), chosen
C  so that
C
C     P1*x = beta*e( n ),
C
C  where e( n ) is the last column of the unit matrix.  P2 is a sequence
C  in planes ( 1, n ), ( 2, n ), ..., ( n - 1, n ), chosen to reduce the
C  matrix ( alpha*beta*e( n )*y' + P1*U ) back to upper triangular form,
C  so that
C
C     P2*( alpha*beta*e( n )*y' + P1*U ) = R
C
C  Q is given as
C
C     Q' = P2*P1.
C
C  The  cosine and sine  that define the rotation in the  ( j, n ) plane
C  for the sequence  P2 are returned in  c( j ) and s( j ) respectively.
C  The  two by two  part  of the  jth  plane rotation  is  of  the  form
C
C     P( j ) = (  c( j )  s( j ) ).
C              ( -s( j )  c( j ) )
C
C  The tangent, t( j ), that defines the rotation in the  ( j, n ) plane
C  for the  sequence  P1  is returned in  x( j )  and the  corresponding
C  cosine and sine  can  be  recovered  by calling the  Nag basic linear
C  algebra routine F06BCF.
C
C  If n is less than unity then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      DOUBLE PRECISION   BETA
C     .. External Subroutines ..
      EXTERNAL           F06FQF, F06QSF, F06QWF, DAXPY
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
C
C        Form  P1  such that  P1*x = beta*e( n ).
C
         BETA = X( ( N - 1 )*INCX + 1 )
         CALL F06FQF( 'Fixed pivot', 'Backward sequence', N - 1, BETA,
     $                X, INCX, C, S )
C
C        Apply P1 to U, to generate a spiked matrix.
C
         CALL F06QWF( 'Left-hand side', N, 1, N, C, S, A, LDA )
C
C        Form  H = alpha*beta*e( n )*y' + P1*U.
C
         CALL DAXPY( N - 1, ALPHA*BETA, Y, INCY, S, 1 )
         A( N, N ) = ALPHA*BETA*Y( ( N - 1 )*INCY + 1 ) + A( N, N )
C
C        Reduce H to upper triangular form, by applying P2.
C
         CALL F06QSF( 'Left-hand side', N, 1, N, C, S, A, LDA )
      END IF
C
      RETURN
C
C     End of F06QPF. ( SUTR1 )
C
      END
