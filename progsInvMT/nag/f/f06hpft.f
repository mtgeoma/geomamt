      SUBROUTINE F06HPF( N, X, INCX, Y, INCY, C, S )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      COMPLEX*16         C, S
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      COMPLEX*16         X( * ), Y( * )
C     ..
C
C  F06HPF performs the plane rotation
C
C     ( x  y ) = ( x  y )*( c  -conjg( s ) ).
C                         ( s   conjg( c ) )
C
C  Advantage is taken of the case where  aimag( c ) = 0.0  and of the
C  case where  aimag( s ) = 0.0.
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 17-April-1985.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      COMPLEX*16         TEMP
      DOUBLE PRECISION   R
      INTEGER            I, IX, IY
C     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, DIMAG, DBLE
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( S.NE.ZERO ).OR.( C.NE.ONE ) )THEN
            IF( DIMAG( C ).EQ.DBLE( ZERO ) )THEN
               R = DBLE( C )
               IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
                  DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     TEMP    = X( IX )
                     X( IX ) = S*Y( IX ) +         R  *TEMP
                     Y( IX ) = R*Y( IX ) - DCONJG( S )*TEMP
   10             CONTINUE
               ELSE
                  IF( INCY.GE.0 )THEN
                     IY = 1
                  ELSE
                     IY = 1 - ( N - 1 )*INCY
                  END IF
                  IF( INCX.GT.0 )THEN
                     DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
                        TEMP    = X( IX )
                        X( IX ) = S*Y( IY ) +         R  *TEMP
                        Y( IY ) = R*Y( IY ) - DCONJG( S )*TEMP
                        IY      = IY        + INCY
   20                CONTINUE
                  ELSE
                     IX = 1 - ( N - 1 )*INCX
                     DO 30, I = 1, N
                        TEMP    = X( IX )
                        X( IX ) = S*Y( IY ) +         R  *TEMP
                        Y( IY ) = R*Y( IY ) - DCONJG( S )*TEMP
                        IX      = IX        + INCX
                        IY      = IY        + INCY
   30                CONTINUE
                  END IF
               END IF
            ELSE IF( DIMAG( S ).EQ.DBLE( ZERO ) )THEN
               R = DBLE( S )
               IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
                  DO 40, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     TEMP    = X( IX )
                     X( IX ) =         R  *Y( IX ) + C*TEMP
                     Y( IX ) = DCONJG( C )*Y( IX ) - R*TEMP
   40             CONTINUE
               ELSE
                  IF( INCY.GE.0 )THEN
                     IY = 1
                  ELSE
                     IY = 1 - ( N - 1 )*INCY
                  END IF
                  IF( INCX.GT.0 )THEN
                     DO 50, IX = 1, 1 + ( N - 1 )*INCX, INCX
                        TEMP    = X( IX )
                        X( IX ) =         R  *Y( IY ) + C*TEMP
                        Y( IY ) = DCONJG( C )*Y( IY ) - R*TEMP
                        IY      = IY                  + INCY
   50                CONTINUE
                  ELSE
                     IX = 1 - ( N - 1 )*INCX
                     DO 60, I = 1, N
                        TEMP    = X( IX )
                        X( IX ) =         R  *Y( IY ) + C*TEMP
                        Y( IY ) = DCONJG( C )*Y( IY ) - R*TEMP
                        IX      = IX                  + INCX
                        IY      = IY                  + INCY
   60                CONTINUE
                  END IF
               END IF
            ELSE
               IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
                  DO 70, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     TEMP    = X( IX )
                     X( IX ) =         S  *Y( IX ) +         C  *TEMP
                     Y( IX ) = DCONJG( C )*Y( IX ) - DCONJG( S )*TEMP
   70             CONTINUE
               ELSE
                  IF( INCY.GE.0 )THEN
                     IY = 1
                  ELSE
                     IY = 1 - ( N - 1 )*INCY
                  END IF
                  IF( INCX.GT.0 )THEN
                     DO 80, IX = 1, 1 + ( N - 1 )*INCX, INCX
                        TEMP    = X( IX )
                        X( IX ) =         S  *Y( IY ) +         C  *TEMP
                        Y( IY ) = DCONJG( C )*Y( IY ) - DCONJG( S )*TEMP
                        IY      = IY                  + INCY
   80                CONTINUE
                  ELSE
                     IX = 1 - ( N - 1 )*INCX
                     DO 90, I = 1, N
                        TEMP    = X( IX )
                        X( IX ) =         S  *Y( IY ) +         C  *TEMP
                        Y( IY ) = DCONJG( C )*Y( IY ) - DCONJG( S )*TEMP
                        IX      = IX                  + INCX
                        IY      = IY                  + INCY
   90                CONTINUE
                  END IF
               END IF
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06HPF. ( CROT )
C
      END
