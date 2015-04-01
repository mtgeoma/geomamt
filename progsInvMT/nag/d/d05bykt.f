      SUBROUTINE D05BYK(LENP,IS,SW,LSW,GAMFUN,TEMPST,VANMAT)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ----------------------------------------------------------------
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     D05BYK calculates the values of the starting fractional weights
C     by solving a Vandermonde system. The right-hand side of the
C     system is assumed to have been already calculated. On entry SW
C     contains the values of the right-hand side and on exit these
C     values are overwritten by the values of the starting weights.
C     D05BYM.
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     ----------------------------------------------------------------
C
C     .. Scalar Arguments ..
      INTEGER           IS, LENP, LSW
C     .. Array Arguments ..
      DOUBLE PRECISION  GAMFUN(IS), SW(LSW,0:IS-1), TEMPST(IS),
     *                  VANMAT(IS)
C     .. Local Scalars ..
      INTEGER           M, N1, NSTAGE
C     .. External Subroutines ..
      EXTERNAL          D05BYM
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Executable Statements ..
C
      N1 = 2**(LENP+1)
C
      DO 20 M = 1, IS
         VANMAT(M) = SQRT(DBLE(M-1))
   20 CONTINUE
C
      DO 80 NSTAGE = 1, IS
C
C        ... Setting the right hand side, NSTAGE=1, IS ...
C
         DO 40 M = 1, IS
            TEMPST(M) = GAMFUN(M)*((NSTAGE-1)**(0.5D0*M))
   40    CONTINUE
C
         CALL D05BYM(VANMAT,TEMPST,IS)
C
         DO 60 M = 0, IS - 1
            SW(NSTAGE,M) = TEMPST(M+1)
   60    CONTINUE
C
   80 CONTINUE
C
      DO 140 NSTAGE = IS + 1, N1 + IS
C
C        ... Setting the right hand side, NSTAGE=IS+1, N+IS ...
C
         DO 100 M = 1, IS
            TEMPST(M) = GAMFUN(M)*(DBLE((NSTAGE-1))**(0.5D0*M)) -
     *                  SW(NSTAGE,M-1)
  100    CONTINUE
C
         CALL D05BYM(VANMAT,TEMPST,IS)
C
         DO 120 M = 0, IS - 1
            SW(NSTAGE,M) = TEMPST(M+1)
  120    CONTINUE
C
  140 CONTINUE
C
      RETURN
      END
