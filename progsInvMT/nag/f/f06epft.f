      SUBROUTINE F06EPF( N, X, INCX, Y, INCY, C, S )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      DROT  ( N, X, INCX, Y, INCY, C, S )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   C, S
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * )
C     ..
C
C  F06EPF performs the plane rotation
C
C     ( x  y ) = ( x  y )*( c  -s ).
C                         ( s   c )
C
C
C  Nag Fortran 77 version of the Blas routine DROT.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 23-January-1984.
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
         IF( ( S.NE.ZERO ).OR.( C.NE.ONE ) )THEN
            IF( ( C.EQ.ZERO ).AND.( S.EQ.ONE ) )THEN
               IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
                  DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     TEMP1   = -X( IX )
                     X( IX ) =  Y( IX )
                     Y( IX ) =  TEMP1
   10             CONTINUE
               ELSE
                  IF( INCY.GE.0 )THEN
                     IY = 1
                  ELSE
                     IY = 1 - ( N - 1 )*INCY
                  END IF
                  IF( INCX.GT.0 )THEN
                     DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
                        TEMP1   = -X( IX )
                        X( IX ) =  Y( IY )
                        Y( IY ) =  TEMP1
                        IY      =  IY       + INCY
   20                CONTINUE
                  ELSE
                     IX = 1 - ( N - 1 )*INCX
                     DO 30, I = 1, N
                        TEMP1   = -X( IX )
                        X( IX ) =  Y( IY )
                        Y( IY ) =  TEMP1
                        IX      =  IX      + INCX
                        IY      =  IY      + INCY
   30                CONTINUE
                  END IF
               END IF
            ELSE IF( ( C.EQ.ZERO ).AND.( S.EQ.( -ONE ) ) )THEN
               IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
                  DO 40, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     TEMP1   =  X( IX )
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
                        TEMP1   =  X( IX )
                        X( IX ) = -Y( IY )
                        Y( IY ) =  TEMP1
                        IY      =  IY       + INCY
   50                CONTINUE
                  ELSE
                     IX = 1 - ( N - 1 )*INCX
                     DO 60, I = 1, N
                        TEMP1   =  X( IX )
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
                     X( IX ) = S*Y( IX ) + C*TEMP1
                     Y( IX ) = C*Y( IX ) - S*TEMP1
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
                        X( IX ) = S*Y( IY ) + C*TEMP1
                        Y( IY ) = C*Y( IY ) - S*TEMP1
                        IY      = IY        + INCY
   80                CONTINUE
                  ELSE
                     IX = 1 - ( N - 1 )*INCX
                     DO 90, I = 1, N
                        TEMP1   = X( IX )
                        X( IX ) = S*Y( IY ) + C*TEMP1
                        Y( IY ) = C*Y( IY ) - S*TEMP1
                        IX      = IX        + INCX
                        IY      = IY        + INCY
   90                CONTINUE
                  END IF
               END IF
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06EPF. ( DROT )
C
      END
