      DOUBLE PRECISION FUNCTION F06FKF( N, W, INCW, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      INTEGER                           INCW, INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION                  W( * ), X( * )
C     ..
C
C  F06FKF returns the weighted euclidean norm of a vector via the
C  function name, so that
C
C     F06FKF := sqrt( x'*W*x ),   where   W = diag( w ).
C
C  The elements of w are assumed to be non-negative.
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 25-June-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      ABSYI, NORM, SCALE, SSQ
      INTEGER               I, IW, IX
C     .. External Functions ..
      DOUBLE PRECISION      F06BMF
      EXTERNAL              F06BMF
C     .. Intrinsic Functions ..
      INTRINSIC             ABS, SQRT
C     ..
C     .. Executable Statements ..
      IF( N.LT.1 )THEN
         NORM  = ZERO
      ELSE IF( N.EQ.1 )THEN
         NORM  = SQRT( W( 1 ) )*ABS( X( 1 ) )
      ELSE
         IF( INCW.GT.0 )THEN
            IW = 1
         ELSE
            IW = 1 - ( N - 1 )*INCW
         END IF
         IF( INCX.GT.0 )THEN
            IX = 1
         ELSE
            IX = 1 - ( N - 1 )*INCX
         END IF
         SCALE = ZERO
         SSQ   = ONE
         DO 10, I = 1, N
            IF( ( W( IW ).NE.ZERO ).AND.( X( IX ).NE.ZERO ) )THEN
               ABSYI = SQRT( W( IW ) )*ABS( X( IX ) )
               IF( SCALE.LT.ABSYI )THEN
                  SSQ   = 1     + SSQ*( SCALE/ABSYI )**2
                  SCALE = ABSYI
               ELSE
                  SSQ   = SSQ   +     ( ABSYI/SCALE )**2
               END IF
            END IF
            IW = IW + INCW
            IX = IX + INCX
   10    CONTINUE
         NORM = F06BMF( SCALE, SSQ )
      END IF
C
      F06FKF = NORM
      RETURN
C
C     End of F06FKF.
C
      END
