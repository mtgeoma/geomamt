      SUBROUTINE F02WEX(NX,NY,ALPHA,X,INCX,Y,INCY,TOL,ZETA)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14A REVISED. IER-687 (DEC 1989).
C
C  F02WEX  generates  details  of  a  Householder reflection  such  that
C
C     P*(   x   ) = (   0  ),   P'*P = I.
C       ( alpha )   ( beta ),
C       (   y   )   (   0  )
C
C  P is given in the form
C
C     P = I - (   w  )*( w'  zeta  z' ),
C             ( zeta )
C             (   z  )
C
C  where  w is an nx element vector, z is an ny element vector, and zeta
C  is a scalar that satisfies
C
C     1.0 .le. zeta .le. sqrt( 2.0 ).
C
C  zeta is returned in ZETA unless the vector v given by
C
C     v = ( x )
C         ( y )
C
C  is such that
C
C     max( abs( v( i ) ) ) .le. max( tol, eps*abs( alpha ) ),
C
C  where  eps  is the  relative machine precision  and  tol  is the user
C  supplied tolerance  TOL, in which case  ZETA  is returned as  0.0 and
C  P  can be taken to be the unit matrix.
C
C  beta  is  overwritten on  alpha,  w  is overwritten on  x  and  z  is
C  overwritten on  y.
C
C  The  routine  may be  called  with  either  or  both  nx = 0, ny = 0.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 17-November-1987.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION ONE, ZERO
      PARAMETER        (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA, TOL, ZETA
      INTEGER          INCX, INCY, NX, NY
C     .. Array Arguments ..
      DOUBLE PRECISION X(*), Y(*)
C     .. Local Scalars ..
      DOUBLE PRECISION BETA, EPS, SCALE, SSQ
      LOGICAL          FIRST
C     .. External Functions ..
      DOUBLE PRECISION X02AJF
      EXTERNAL         X02AJF
C     .. External Subroutines ..
      EXTERNAL         F06FJF, F06FRF, DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC        ABS, MAX, SQRT
C     .. Save statement ..
      SAVE             EPS, FIRST
C     .. Data statements ..
      DATA             FIRST/.TRUE./
C     .. Executable Statements ..
      IF (NX.LT.1) THEN
         CALL F06FRF(NY,ALPHA,Y,INCY,TOL,ZETA)
      ELSE IF (NY.LT.1) THEN
         CALL F06FRF(NX,ALPHA,X,INCX,TOL,ZETA)
      ELSE
C
         IF (FIRST) THEN
            FIRST = .FALSE.
            EPS = X02AJF()
         END IF
C
         SSQ = ONE
         SCALE = ZERO
         CALL F06FJF(NX,X,INCX,SCALE,SSQ)
         CALL F06FJF(NY,Y,INCY,SCALE,SSQ)
C
C        Treat  cases  where   SCALE = zero,   SCALE is negligible   and
C        ALPHA = zero  specially.  Note that
C
C           SCALE = max( abs( v( i ) ) ).
C
         IF ((SCALE.EQ.ZERO) .OR. (SCALE.LE.MAX(TOL,EPS*ABS(ALPHA))))
     *       THEN
            ZETA = ZERO
         ELSE IF (ALPHA.EQ.ZERO) THEN
            ZETA = ONE
            ALPHA = SCALE*SQRT(SSQ)
            CALL DSCAL(NX,-1/ALPHA,X,INCX)
            CALL DSCAL(NY,-1/ALPHA,Y,INCY)
         ELSE
            IF (SCALE.LT.ABS(ALPHA)) THEN
               BETA = ABS(ALPHA)*SQRT(1+SSQ*(SCALE/ALPHA)**2)
            ELSE
               BETA = SCALE*SQRT(SSQ+(ALPHA/SCALE)**2)
            END IF
            ZETA = SQRT((BETA+ABS(ALPHA))/BETA)
            IF (ALPHA.GT.ZERO) BETA = -BETA
            CALL DSCAL(NX,-1/(ZETA*BETA),X,INCX)
            CALL DSCAL(NY,-1/(ZETA*BETA),Y,INCY)
            ALPHA = BETA
         END IF
      END IF
C
      RETURN
C
C     End of F02WEX. ( SGRFG2 )
C
      END
