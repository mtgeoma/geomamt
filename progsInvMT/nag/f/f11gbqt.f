      SUBROUTINE F11GBQ(NEXT,IREVCM,ITN,RESID,PRECON,NORM,ITERM,N,TOL,
     *                  SIGMAX,ANORM,XNORM0,BNORM,BNORM2,RNORM2,BETAIN,
     *                  BETA1,BETA2,PI,GAMMAB,RHO1,RHO2,ZETA,ZETAB,B,X,
     *                  V1,V2,W1,WGT,U,V,STPLHS,STPRHS,USEXC,INFOCH)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-----------------------------------------------------------------------
C
C     F11GBQ - Termination criterion (Lanczos Method (SYMMLQ))
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO, TEN
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,TEN=1.0D1)
      INTEGER           NLOOP
      PARAMETER         (NLOOP=2)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ANORM, BETA1, BETA2, BETAIN, BNORM, BNORM2,
     *                  GAMMAB, PI, RHO1, RHO2, RNORM2, SIGMAX, STPLHS,
     *                  STPRHS, TOL, XNORM0, ZETA, ZETAB
      INTEGER           INFOCH, IREVCM, ITERM, ITN, N, NORM, RESID
      LOGICAL           NEXT, PRECON, USEXC
C     .. Array Arguments ..
      DOUBLE PRECISION  B(N), U(N), V(N), V1(N), V2(N), W1(N), WGT(N),
     *                  X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPSLON, FACTOL, RCBNRM, RLBNRM, STPRHX, TOLC,
     *                  TOLEPS, TOLF, TOLL, XCBNRM, XLBNRM, XNORM
      INTEGER           LOOP
      LOGICAL           FTERM
C     .. External Functions ..
      DOUBLE PRECISION  A02ABF, F11BBU, X02AJF
      EXTERNAL          A02ABF, F11BBU, X02AJF
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX, SQRT
C     .. Save statement ..
      SAVE              FACTOL, STPRHX, TOLC, TOLEPS, TOLF, TOLL,
     *                  XCBNRM, XLBNRM, XNORM, LOOP, FTERM
C     .. Executable Statements ..
C
C     ITERM = 1 :  Using ||r|| <= tol*(||b|| + ||A||*||x||)
C
      IF (ITERM.LE.1) THEN
C
C        First call to this routine
C
         IF (NEXT) THEN
            RESID = 0
            IF (ITN.LE.0) THEN
               FTERM = .TRUE.
               TOLF = MAX(TOL,MAX(TEN,SQRT(DBLE(N)))*X02AJF())
            END IF
            IF (FTERM) THEN
               TOLEPS = TOL/X02AJF()
               STPRHX = TOL*(BNORM+ANORM*XNORM0)
               XLBNRM = ZERO
               XCBNRM = ZERO
            END IF
            IF (INFOCH.EQ.1) LOOP = 0
C
C           Compute ||x_bar_c[k+1]||, ||x_bar_l[k]||, ||r_bar_c[k+1]||,
C           ||r_bar_l[k]|| and the tolerances
C
            XLBNRM = A02ABF(XLBNRM,ZETA)
            XCBNRM = A02ABF(XLBNRM,ZETAB)
            RCBNRM = BETAIN*BETA2*ABS(PI/GAMMAB)
            RLBNRM = A02ABF(RHO1,RHO2)
            TOLL = TOL*(BETAIN+SIGMAX*XLBNRM)
            TOLC = TOL*(BETAIN+SIGMAX*XCBNRM)
            USEXC = ((RCBNRM/TOLC).LE.(RLBNRM/TOLL))
C
C           Use the approximations:
C
C           ||r_c[k+1]|| = ||(r_c_bar[k+1]/beta[k+2])*v[k+2])  or
C
C           ||r_l[k]|| = ||(rho[1]/beta[k+1])*v[k+1]) +
C                        (rho[2]/beta[k+2])*v[k+2])||
C
            IF (USEXC) THEN
               CALL DCOPY(N,X,1,U,1)
               CALL DAXPY(N,ZETAB,W1,1,U,1)
               XNORM = F11BBU(NORM,N,U,WGT)
               STPLHS = (BETAIN*ABS(PI/GAMMAB))*F11BBU(NORM,N,V2,WGT)
               STPRHS = TOL*(BNORM+ANORM*XNORM)
            ELSE
               XNORM = F11BBU(NORM,N,X,WGT)
               IF (ABS(RHO2*BETA1).GT.ZERO) THEN
                  CALL DCOPY(N,V2,1,V,1)
                  CALL DAXPY(N,(RHO1*BETA2)/(RHO2*BETA1),V1,1,V,1)
                  STPLHS = ABS(RHO2/BETA2)*F11BBU(NORM,N,V,WGT)
               ELSE
                  STPLHS = ABS(RHO1/BETA1)*F11BBU(NORM,N,V1,WGT)
               END IF
               STPRHS = TOL*(BNORM+ANORM*XNORM)
            END IF
            STPRHX = MAX(STPRHX,STPRHS)
            IF (STPLHS.LE.MAX(STPRHS,STPRHX)) THEN
               IF (STPLHS.LE.STPRHS) THEN
                  NEXT = .FALSE.
                  INFOCH = 3
                  IREVCM = 1
                  IF ( .NOT. USEXC) CALL DCOPY(N,X,1,U,1)
               ELSE IF (XNORM0/(XNORM*TOLEPS).GT.TWO) THEN
                  INFOCH = -1
                  IF (USEXC) CALL DCOPY(N,U,1,X,1)
               ELSE
                  INFOCH = 1
               END IF
            ELSE
               INFOCH = 1
            END IF
            FTERM = INFOCH .LT. 0
C
C        Subsequent calls to this routine
C
         ELSE
            NEXT = .TRUE.
            RESID = 2
            CALL DSCAL(N,-ONE,V,1)
            CALL DAXPY(N,ONE,B,1,V,1)
            STPLHS = F11BBU(NORM,N,V,WGT)
            STPRHS = TOLF*(BNORM+ANORM*XNORM)
            LOOP = LOOP + 1
            IF (STPLHS.LE.STPRHS) THEN
               INFOCH = 0
            ELSE IF (LOOP.LT.NLOOP) THEN
               IF (XNORM0/(XNORM*TOLEPS).GT.1.0625D0) THEN
                  INFOCH = -2
                  IF (USEXC) CALL DCOPY(N,U,1,X,1)
               ELSE
                  INFOCH = 2
               END IF
            ELSE
               INFOCH = 4
            END IF
         END IF
C
C     ITERM = 2 :  using ||r_bar|| <= tol*(beta[1]+||A_bar||*||Dx||)
C
      ELSE
C
C        First call to this routine
C
         IF (NEXT) THEN
            RESID = 0
            IF (ITN.LE.0) THEN
               EPSLON = X02AJF()
               TOLF = MAX(TOL,TEN*EPSLON,SQRT(DBLE(N))*EPSLON)
               IF ((BNORM2.LT.RNORM2) .OR. (RNORM2.EQ.ZERO)) THEN
                  FACTOL = ONE
               ELSE
                  FACTOL = BNORM2/RNORM2
               END IF
               XLBNRM = ZERO
               XCBNRM = ZERO
            END IF
C
C           Compute ||x_bar_c[k+1]||, ||x_bar_l[k]||, ||r_bar_c[k+1]||,
C           ||r_bar_l[k]|| and the tolerances
C
            XLBNRM = A02ABF(XLBNRM,ZETA)
            XCBNRM = A02ABF(XLBNRM,ZETAB)
            RCBNRM = BETAIN*BETA2*ABS(PI/GAMMAB)
            RLBNRM = A02ABF(RHO1,RHO2)
            TOLL = FACTOL*TOLF*(BETAIN+SIGMAX*XLBNRM)
            TOLC = FACTOL*TOLF*(BETAIN+SIGMAX*XCBNRM)
            USEXC = ((RCBNRM/TOLC).LE.(RLBNRM/TOLL))
            IF (USEXC) THEN
               STPLHS = RCBNRM
               STPRHS = TOLC
            ELSE
               STPLHS = RLBNRM
               STPRHS = TOLL
            END IF
            IF (STPLHS.LE.STPRHS) THEN
               NEXT = .FALSE.
               RESID = 2
               IF (USEXC) CALL DAXPY(N,ZETAB,W1,1,X,1)
               CALL DCOPY(N,X,1,U,1)
               INFOCH = 3
               IREVCM = 1
            ELSE
               INFOCH = 1
            END IF
C
C        Subsequent calls to this routine
C
         ELSE
            NEXT = .TRUE.
            INFOCH = 0
            CALL DSCAL(N,-ONE,V,1)
            CALL DAXPY(N,ONE,B,1,V,1)
         END IF
      END IF
C
C     End of subroutine F11GBQ
C
      RETURN
      END
