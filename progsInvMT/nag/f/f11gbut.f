      SUBROUTINE F11GBU(NEXT,IREVCM,ITN,PRECON,NORM,ITERM,TOL,ANORM,
     *                  BNORM,XNORM,SIGMAX,N,B,X,WGT,U,V,STPLHS,STPRHS,
     *                  INFOCH)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-----------------------------------------------------------------------
C
C     F11GBU - Termination criterion (Conjugate Gradient Method)
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO, TEN
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,TEN=10.0D0)
      INTEGER           NLOOP
      PARAMETER         (NLOOP=2)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ANORM, BNORM, SIGMAX, STPLHS, STPRHS, TOL, XNORM
      INTEGER           INFOCH, IREVCM, ITERM, ITN, N, NORM
      LOGICAL           NEXT, PRECON
C     .. Array Arguments ..
      DOUBLE PRECISION  B(N), U(N), V(N), WGT(N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPSLON, STPRHX, TOLEPS, TOLF, XNORM0
      INTEGER           LOOP, NRESTC, NRESTL
      LOGICAL           FTERM
C     .. External Functions ..
      DOUBLE PRECISION  F11BBU, X02AJF
      EXTERNAL          F11BBU, X02AJF
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX, SQRT
C     .. Save statement ..
      SAVE              EPSLON, STPRHX, TOLEPS, TOLF, XNORM0, LOOP,
     *                  FTERM, NRESTL, NRESTC
C     .. Executable Statements ..
C
C     First calls to this routine
C
      IF (NEXT) THEN
         IF (ITN.LE.0) THEN
            NRESTL = 0
            NRESTC = 0
            FTERM = .TRUE.
            EPSLON = X02AJF()
            TOLF = MAX(TOL,TEN*EPSLON,SQRT(DBLE(N))*EPSLON)
         END IF
         IF (FTERM) THEN
            TOLEPS = TOL/EPSLON
            XNORM0 = XNORM
            IF (ITERM.LE.1) THEN
               STPRHX = TOL*(BNORM+ANORM*XNORM0)
            ELSE
               STPRHX = ZERO
            END IF
         END IF
         IF (INFOCH.EQ.1) LOOP = 0
         IF (ITERM.LE.1) THEN
            STPRHS = TOL*(BNORM+ANORM*XNORM)
            STPRHX = MAX(STPRHX,STPRHS)
         ELSE
            STPRHS = TOL*SIGMAX*XNORM
         END IF
         IF (STPLHS.LE.MAX(STPRHS,STPRHX)) THEN
            IF (ITN.LE.0) THEN
               INFOCH = 0
            ELSE IF ((ITERM.GE.2) .OR. (STPLHS.LE.STPRHS)) THEN
               NEXT = .FALSE.
               INFOCH = 3
               IREVCM = 1
               CALL DCOPY(N,X,1,U,1)
            ELSE IF ((XNORM0/XNORM)/TOLEPS.GT.TWO) THEN
               NRESTL = NRESTL + 1
               INFOCH = -1
            ELSE
               INFOCH = 1
            END IF
         ELSE
            INFOCH = 1
         END IF
C
C     Subsequent calls to this routine
C
      ELSE
         NEXT = .TRUE.
         CALL DSCAL(N,-ONE,V,1)
         CALL DAXPY(N,ONE,B,1,V,1)
         IF (ITERM.LE.1) THEN
            STPLHS = F11BBU(NORM,N,V,WGT)
            STPRHS = TOLF*(BNORM+ANORM*XNORM)
         END IF
         LOOP = LOOP + 1
         IF (STPLHS.LE.STPRHS) THEN
            INFOCH = 0
         ELSE IF (LOOP.LT.NLOOP) THEN
            IF ((XNORM0/XNORM)/TOLEPS.GT.1.0625D0) THEN
               NRESTC = NRESTC + 1
               INFOCH = -2
            ELSE
               INFOCH = 2
            END IF
         ELSE
            INFOCH = 4
         END IF
      END IF
C
      FTERM = (INFOCH.LT.0)
C
C     End of subroutine F11GBU
C
      RETURN
      END
