      SUBROUTINE F02XEX(NX,NY,ALPHA,X,INCX,Y,INCY,TOL,THETA)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C  F02XEX generates details of a generalized Householder reflection such
C  that
C
C     P*(   x   ) = (   0  ),   conjg( P' )*P = I,  aimag( beta ) = 0.0.
C       ( alpha )   ( beta ),
C       (   y   )   (   0  )
C
C  P is given in the form
C
C     P = I - gamma*(   w  )*( conjg( w' )  zeta  conjg( z' ) ),
C                   ( zeta )
C                   (   z  )
C
C  where  w is an nx element vector, z is an ny element vector, gamma is
C  a scalar such that
C
C     real ( gamma ) = 1.0,
C     aimag( gamma ) = aimag( alpha )/( beta - real( alpha ) )
C
C  and zeta is a real scalar that satisfies
C
C     1.0 .le. zeta .le. sqrt( 2.0 ).
C
C  Note that when alpha is real then gamma = 1.0.
C
C  gamma and zeta are returned in THETA as
C
C     THETA = ( zeta, aimag( gamma ) )
C
C  unless the vector v given by
C
C     v = ( x )
C         ( y )
C
C  is such that
C
C     max( abs( real( v( i ) ) ), abs( aimag( v( i ) ) ) ) .le.
C     max( tol, eps*max( abs( real ( alpha  ) ),
C                        abs( aimag( alpha  ) ) ) ),
C
C  where eps is the relative machine precision and tol is the user
C  supplied tolerance  TOL,  in which case  THETA is returned as 0.0, or
C  THETA  is such that  real( THETA ) .le. 0.0, in which case  P  can be
C  taken to be
C
C     P = I              when   THETA = 0.0,
C
C     P = ( THETA  0 )   when   real( THETA ) .le. 0.0, THETA .ne. 0.0.
C         (   0    I )
C
C  beta is overwritten on alpha with the imaginary part of alpha set to
C  zero, w is overwritten on x and z is overwritten on y.
C
C  The routine may be called with either or both  nx = 0, ny = 0.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 17-November-1987.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=(0.0D+0,0.0D+0))
C     .. Scalar Arguments ..
      COMPLEX*16        ALPHA, THETA
      DOUBLE PRECISION  TOL
      INTEGER           INCX, INCY, NX, NY
C     .. Array Arguments ..
      COMPLEX*16        X(*), Y(*)
C     .. Local Scalars ..
      COMPLEX*16        GAMMA
      DOUBLE PRECISION  BETA, EPS, SCALE, SSQ, ZETA
      LOGICAL           FIRST
C     .. Local Arrays ..
      COMPLEX*16        WORK(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          ZSCAL, ZDSCAL, F06HRF, F06KJF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DCMPLX, DCONJG, MAX, DBLE, SIGN,
     *                  SQRT
C     .. Save statement ..
      SAVE              EPS, FIRST
C     .. Data statements ..
      DATA              FIRST/.TRUE./
C     .. Executable Statements ..
      IF (NX.LT.1) THEN
         CALL F06HRF(NY,ALPHA,Y,INCY,TOL,THETA)
      ELSE IF (NY.LT.1) THEN
         CALL F06HRF(NX,ALPHA,X,INCX,TOL,THETA)
      ELSE
C
         IF (FIRST) THEN
            FIRST = .FALSE.
            EPS = X02AJF()
         END IF
C
         SSQ = ONE
         SCALE = DBLE(ZERO)
         CALL F06KJF(NX,X,INCX,SCALE,SSQ)
         CALL F06KJF(NY,Y,INCY,SCALE,SSQ)
C
C        Treat cases where  SCALE = zero,  SCALE is negligible
C        and  ALPHA = zero  specially.
C        Note that
C        SCALE = max( abs( real( v( i ) ) ), abs( aimag( v( i ) ) ) ).
C
         IF ((SCALE.EQ.DBLE(ZERO)) .OR. (SCALE.LE.MAX(TOL,
     *       EPS*MAX(ABS(DBLE(ALPHA)),ABS(DIMAG(ALPHA)))))) THEN
            IF (DIMAG(ALPHA).EQ.DBLE(ZERO)) THEN
               THETA = ZERO
            ELSE
               BETA = -SIGN(ABS(ALPHA),DBLE(ALPHA))
               THETA = DCONJG(ALPHA)/BETA
               ALPHA = BETA
            END IF
         ELSE IF (ALPHA.EQ.ZERO) THEN
            THETA = ONE
            BETA = SCALE*SQRT(SSQ)
            CALL ZDSCAL(NX,-1/BETA,X,INCX)
            CALL ZDSCAL(NY,-1/BETA,Y,INCY)
            ALPHA = BETA
         ELSE
            WORK(1) = ALPHA
            CALL F06KJF(1,WORK,1,SCALE,SSQ)
            BETA = SCALE*SQRT(SSQ)
            ZETA = SQRT((BETA+ABS(DBLE(ALPHA)))/BETA)
            IF (DBLE(ALPHA).GT.DBLE(ZERO)) BETA = -BETA
            IF (DIMAG(ALPHA).EQ.DBLE(ZERO)) THEN
               CALL ZDSCAL(NX,-1/(ZETA*BETA),X,INCX)
               CALL ZDSCAL(NY,-1/(ZETA*BETA),Y,INCY)
               THETA = ZETA
               ALPHA = BETA
            ELSE
               GAMMA = DCMPLX(ONE,DIMAG(ALPHA)/(BETA-DBLE(ALPHA)))
               CALL ZSCAL(NX,-1/(DCONJG(GAMMA)*ZETA*BETA),X,INCX)
               CALL ZSCAL(NY,-1/(DCONJG(GAMMA)*ZETA*BETA),Y,INCY)
               THETA = DCMPLX(ZETA,DIMAG(GAMMA))
               ALPHA = BETA
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F02XEX. ( CGRFG2 )
C
      END
