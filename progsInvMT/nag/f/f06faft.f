      DOUBLE PRECISION FUNCTION F06FAF( N, J, TOLX, X, INCX,
     $                                  TOLY, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                  TOLX, TOLY
      INTEGER                           INCX, INCY, J, N
C     .. Array Arguments ..
      DOUBLE PRECISION                  X( * ), Y( * )
C     ..
C
C  F06FAF returns the cosine, vcos, of the angle between the vectors x
C  and y, given by
C
C     vcos = ( x'*y )/( norm( x )*norm( y ) ),
C
C  where
C
C     norm( z ) = sqrt( z'*z ).
C
C  cos( theta ) is returned via the function name.
C
C  When  1.le.j.le.n  then y is taken to be the unit vector e( j ), in
C  which case the array Y is not referenced. Otherwise y must be
C  supplied in the incremented array Y.
C
C  If  norm( x ).le.tolx  then vcos is returned as  2.0 otherwise
C  if  norm( y ).le.toly  then vcos is returned as -2.0, otherwise vcos
C  is returned in the range ( -1.0, 1.0 ).
C
C  If tol is negative then the value zero is used in place of tol.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ONE         , TWO         , ZERO
      PARAMETER           ( ONE = 1.0D+0, TWO = 2.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      C, XNORM, YNORM
      INTEGER               I, IX, IY
C     .. External Functions ..
      DOUBLE PRECISION      DNRM2
      EXTERNAL              DNRM2
C     .. Intrinsic Functions ..
      INTRINSIC             ABS, MAX, DBLE, SQRT
C     ..
C     .. Executable Statements ..
      IF( N.LT.1 )THEN
         C = TWO
      ELSE
         IF( INCX.EQ.0 )THEN
            XNORM = SQRT( DBLE( N ) )*ABS( X( 1 ) )
         ELSE
            XNORM = DNRM2( N, X, ABS( INCX ) )
         END IF
         IF( XNORM.LE.MAX( TOLX, ZERO ) )THEN
            C = TWO
         ELSE
            IF( INCX.GT.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( ( J.GT.0 ).AND.( J.LE.N ) )THEN
               C = X( IX + ( J - 1 )*INCX )/XNORM
               IF( C.GT.ONE )THEN
                  C = ONE
               ELSE IF( C.LT.( -ONE ) )THEN
                  C = -ONE
               END IF
            ELSE
               IF( INCY.EQ.0 )THEN
                  YNORM = SQRT( DBLE( N ) )*ABS( Y( 1 ) )
               ELSE
                  YNORM = DNRM2( N, Y, ABS( INCY ) )
               END IF
               IF( YNORM.LE.MAX( TOLY, ZERO ) )THEN
                  C = -TWO
               ELSE
                  IF( INCY.GT.0 )THEN
                     IY = 1
                  ELSE
                     IY = 1 - ( N - 1 )*INCY
                  END IF
                  C = ZERO
                  DO 10, I = 1, N
                     C  = C  + ( X( IX )/XNORM )*( Y( IY )/YNORM )
                     IX = IX + INCX
                     IY = IY + INCY
   10             CONTINUE
                  IF( C.GT.ONE )THEN
                     C =  ONE
                  ELSE IF( C.LT.( -ONE ) )THEN
                     C = -ONE
                  END IF
               END IF
            END IF
         END IF
      END IF
C
      F06FAF = C
      RETURN
C
C     End of F06FAF. ( SVCOS )
C
      END
