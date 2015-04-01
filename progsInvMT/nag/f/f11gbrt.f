      SUBROUTINE F11GBR(NEXT,FLOOP,IREVCM,ITN,PRECON,N,ALPHA,BETA1,
     *                  BETA2,PI,GAMMAB,RHO1,RHO2,ZETA,ZETAB,X,V1,V2,W1,
     *                  W2,U,V,INFO)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-----------------------------------------------------------------------
C
C     F11GBR - Body of the loop (Lanczos Method (SYMMLQ))
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA, BETA1, BETA2, GAMMAB, PI, RHO1, RHO2,
     *                  ZETA, ZETAB
      INTEGER           INFO, IREVCM, ITN, N
      LOGICAL           FLOOP, NEXT, PRECON
C     .. Array Arguments ..
      DOUBLE PRECISION  U(N), V(N), V1(N), V2(N), W1(N), W2(N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  C, DELTA, DELTAB, EPSLON, GAMMA, S, T
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DSCAL, F06BAF, F06FPF
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Save statement ..
      SAVE              DELTAB, EPSLON
C     .. Executable Statements ..
C
C     First call to F11GBR
C
      IF (NEXT) THEN
         ITN = ITN + 1
         NEXT = .FALSE.
         IREVCM = 1
         IF (FLOOP) THEN
            DELTAB = BETA2
            FLOOP = .FALSE.
         END IF
         CALL DCOPY(N,W2,1,U,1)
C
C     Subsequent calls to F11GBR
C
      ELSE
         IF (IREVCM.LE.1) THEN
            CALL DAXPY(N,-BETA2/BETA1,V1,1,V,1)
            ALPHA = DDOT(N,W2,1,V,1)
            CALL DAXPY(N,-ALPHA/BETA2,V2,1,V,1)
            CALL DCOPY(N,V2,1,V1,1)
            CALL DCOPY(N,V,1,V2,1)
            IF (PRECON) THEN
               IREVCM = 2
               CALL DCOPY(N,V,1,U,1)
            ELSE
               NEXT = .TRUE.
            END IF
         ELSE
            NEXT = .TRUE.
         END IF
C
C        Completion
C
         IF (NEXT) THEN
C
C           Compute the rotation parameters c[k], s[k] and compute
C           gamma[k]
C
            GAMMA = GAMMAB
            T = BETA2
            CALL F06BAF(GAMMA,T,C,S)
C
C           Compute gamma_bar[k+1], delta[k+1], zeta[k], update rho_1
C           and compute zeta[k], zeta_bar[k+1]
C
            GAMMAB = S*DELTAB - C*ALPHA
            DELTA = C*DELTAB + S*ALPHA
            ZETA = RHO1/GAMMA
            RHO1 = RHO2 - DELTA*ZETA
            ZETAB = RHO1/GAMMAB
C
C           Compute beta[k+2]**2, set beta[k+1] => BETA1 and
C           beta[k+2] => BETA2
C
            BETA1 = BETA2
            BETA2 = DDOT(N,V2,1,V,1)
C
C           beta[k+2]**2 <= 0.  M is not positive-definite: terminate
C
            IF (BETA2.LE.ZERO) THEN
               INFO = 6
C
C           beta[k+2]**2 > 0.  Continue execution
C
            ELSE
C
C              Compute beta[k+2], delta_bar[k+2], epsilon[k+2] and
C              update pi[k] and rho_2
C
               BETA2 = SQRT(BETA2)
               DELTAB = -C*BETA2
               EPSLON = S*BETA2
               PI = PI*S
               RHO2 = -EPSLON*ZETA
C
C              Normalize z[k+2]
C
               CALL DSCAL(N,ONE/BETA2,V,1)
            END IF
C
C           (w[k] w_bar[k+1]) = (w_bar[k] z[k+1]) * (c[k]  s[k])
C                                                   (s[k] -c[k])
C           w[k] => W1, w_bar[k+1] => W2
C
            CALL F06FPF(N,W1,1,W2,1,C,S)
C
C           Compute x_l[k] = x_l[k-1] + zeta[k]*w[k]
C
            CALL DAXPY(N,ZETA,W1,1,X,1)
C
C           Set w_bar[k+1] => W1, z[k+2] => W2
C
            IF (INFO.LE.0) THEN
               CALL DCOPY(N,W2,1,W1,1)
               CALL DCOPY(N,V,1,W2,1)
            END IF
         END IF
C
      END IF
C
C     End of subroutine F11GBR
C
      RETURN
      END
