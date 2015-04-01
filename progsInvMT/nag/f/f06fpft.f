      SUBROUTINE F06FPF( N, X, INCX, Y, INCY, C, S )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   C, S
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * )
C     ..
C
C  F06FPF performs the symmetric plane rotation
C
C     ( x  y ) = ( x  y )*( c   s ),   s .ne. 0.0.
C                         ( s  -c )
C
C  If s is supplied as zero then x and y are unaltered.
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 21-October-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP1
      INTEGER            I, IX, IY
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( S.NE.ZERO )THEN
            IF( ( C.EQ.ZERO ).AND.( S.EQ.ONE ) )THEN
               IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
                  DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     TEMP1   = X( IX )
                     X( IX ) = Y( IX )
                     Y( IX ) = TEMP1
   10             CONTINUE
               ELSE
                  IF( INCY.GE.0 )THEN
                     IY = 1
                  ELSE
                     IY = 1 - ( N - 1 )*INCY
                  END IF
                  IF( INCX.GT.0 )THEN
                     DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
                        TEMP1   = X( IX )
                        X( IX ) = Y( IY )
                        Y( IY ) = TEMP1
                        IY      = IY      + INCY
   20                CONTINUE
                  ELSE
                     IX = 1 - ( N - 1 )*INCX
                     DO 30, I = 1, N
                        TEMP1   = X( IX )
                        X( IX ) = Y( IY )
                        Y( IY ) = TEMP1
                        IX      = IX      + INCX
                        IY      = IY      + INCY
   30                CONTINUE
                  END IF
               END IF
            ELSE IF( ( C.EQ.ZERO ).AND.( S.EQ.( -ONE ) ) )THEN
               IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
                  DO 40, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     TEMP1   = -X( IX )
                     X( IX ) = -Y( IX )
                     Y( IX ) =  TEMP1
   40             CONTINUE
               ELSE
                  IF( INCY.GE.0 )THEN
                     IY = 1
                  ELSE
                     IY = 1 - ( N - 1 )*INCY
                  END IF
                  IF( INCX.GT.0 )THEN
                     DO 50, IX = 1, 1 + ( N - 1 )*INCX, INCX
                        TEMP1   = -X( IX )
                        X( IX ) = -Y( IY )
                        Y( IY ) =  TEMP1
                        IY      =  IY      + INCY
   50                CONTINUE
                  ELSE
                     IX = 1 - ( N - 1 )*INCX
                     DO 60, I = 1, N
                        TEMP1   = -X( IX )
                        X( IX ) = -Y( IY )
                        Y( IY ) =  TEMP1
                        IX      =  IX      + INCX
                        IY      =  IY      + INCY
   60                CONTINUE
                  END IF
               END IF
            ELSE
               IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
                  DO 70, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     TEMP1   = X( IX )
                     X( IX ) = C*TEMP1 + S*Y( IX )
                     Y( IX ) = S*TEMP1 - C*Y( IX )
   70             CONTINUE
               ELSE
                  IF( INCY.GE.0 )THEN
                     IY = 1
                  ELSE
                     IY = 1 - ( N - 1 )*INCY
                  END IF
                  IF( INCX.GT.0 )THEN
                     DO 80, IX = 1, 1 + ( N - 1 )*INCX, INCX
                        TEMP1   = X( IX )
                        X( IX ) = C*TEMP1 + S*Y( IY )
                        Y( IY ) = S*TEMP1 - C*Y( IY )
                        IY      = IY      + INCY
   80                CONTINUE
                  ELSE
                     IX = 1 - ( N - 1 )*INCX
                     DO 90, I = 1, N
                        TEMP1   = X( IX )
                        X( IX ) = C*TEMP1 + S*Y( IY )
                        Y( IY ) = S*TEMP1 - C*Y( IY )
                        IX      = IX      + INCX
                        IY      = IY      + INCY
   90                CONTINUE
                  END IF
               END IF
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06FPF. ( SROTS )
C
      END
