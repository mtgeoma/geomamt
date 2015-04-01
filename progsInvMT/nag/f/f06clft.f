      COMPLEX*16       FUNCTION F06CLF( A, B, FAIL )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-605 (MAR 1988).
C     .. Scalar Arguments ..
      COMPLEX*16                        A, B
      LOGICAL                           FAIL
C     ..
C
C  F06CLF returns the value div given by
C
C     div = ( a/b      if a/b does not overflow,
C           (
C           ( 0.0      if a .eq. 0.0,
C           (
C           ( cflmax   if a .ne. 0.0 and a/b would overflow,
C
C  where
C
C     cflmax = ( flmax*sign( re( a/b ) ), flmax*sign( im( a/b ) ) )
C
C  and flmax is a large value, via the function name. In addition if a/b
C  would  overflow then  fail  is returned as  true, otherwise  fail  is
C  returned as false.
C
C  Note that when a and b are both zero, fail is returned as .true., but
C  div  is returned as  0.0. In all other cases of overflow  div is such
C  that abs( re( div ) ) and abs( im( div ) ) are flmax.
C
C  For  real  x and y,  if  y = 0,  sign( x/y )  is taken as  sign( x ).
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 27-April-1983.
C     Sven Hammarling, Nag Central Office.
C  -- Amended on 4-December-1987.
C     Sven Hammarling and Mick Pont, Nag Central Office.
C        To avoid extremely unlikely division by zero.
C
C
C     .. Parameters ..
      DOUBLE PRECISION         ONE
      PARAMETER              ( ONE  = 1.0D+0 )
      COMPLEX*16               ZERO
      PARAMETER              ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      COMPLEX*16               VALUE
      DOUBLE PRECISION         AI, AR, BI, BIG, BR, DIV, FLMAX, FLMIN,
     $                         NUMI, NUMR, TEMP
      LOGICAL                  FIRST
C     .. External Functions ..
      DOUBLE PRECISION         X02AMF
      EXTERNAL                 X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC                ABS, DBLE, DCMPLX, DIMAG, MAX, SIGN
C     .. Save statement ..
      SAVE                     BIG, FIRST, FLMAX
C     .. Data statements ..
      DATA                     FIRST/ .TRUE. /
C     ..
C     .. Executable Statements ..
      IF( A.EQ.ZERO )THEN
         VALUE = ZERO
         IF( B.EQ.ZERO )THEN
            FAIL = .TRUE.
         ELSE
            FAIL = .FALSE.
         END IF
      ELSE
C
         IF( FIRST )THEN
            FIRST = .FALSE.
            FLMIN =  X02AMF( )
            FLMAX =  1/FLMIN
            BIG   =  FLMAX/2
         END IF
C
         AR    =  DBLE ( A )
         AI    =  DIMAG( A )
         BR    =  DBLE ( B )
         BI    =  DIMAG( B )
         TEMP  =  MAX( ABS( AR ), ABS( AI ), ABS( BR ), ABS( BI ) )
         IF( TEMP.GE.BIG )THEN
            AR = AR/2
            AI = AI/2
            BR = BR/2
            BI = BI/2
         END IF
         IF( DCMPLX( BR, BI ).EQ.ZERO ) THEN
            VALUE =  DCMPLX( SIGN( FLMAX, DBLE ( A ) ),
     $                       SIGN( FLMAX, DIMAG( A ) )  )
            FAIL  = .TRUE.
         ELSE
            IF( ABS( BR ).GE.ABS( BI ) )THEN
               TEMP = BI/BR
               DIV  = BR     + TEMP*BI
               NUMR = AR     + TEMP*AI
               NUMI = AI     - TEMP*AR
            ELSE
               TEMP = BR/BI
               DIV  = BI      + TEMP*BR
               NUMR = AI      + TEMP*AR
               NUMI = TEMP*AI - AR
            END IF
            IF( ABS( DIV ).GE.ONE )THEN
               VALUE =  DCMPLX( NUMR/DIV, NUMI/DIV )
               FAIL  = .FALSE.
            ELSE
               TEMP  =  ABS( DIV )*FLMAX
               IF( ( ABS( NUMR ).LE.TEMP ).AND.
     $             ( ABS( NUMI ).LE.TEMP )      )THEN
                  VALUE =  DCMPLX( NUMR/DIV, NUMI/DIV )
                  FAIL  = .FALSE.
               ELSE
                  IF( DIV.GE.DBLE( ZERO ) )THEN
                     VALUE = DCMPLX( SIGN( FLMAX,  NUMR ),
     $                               SIGN( FLMAX,  NUMI )  )
                  ELSE
                     VALUE = DCMPLX( SIGN( FLMAX, -NUMR ),
     $                               SIGN( FLMAX, -NUMI )  )
                  END IF
                  FAIL = .TRUE.
               END IF
            END IF
         END IF
      END IF
C
      F06CLF = VALUE
      RETURN
C
C     End of F06CLF. ( CDIV )
C
      END
