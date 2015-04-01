      SUBROUTINE F11GBS(NEXT,IREVCM,PRECON,N,ALPHA,BETA1,BNORM2,RNORM2,
     *                  X,V1,V2,W1,W2,U,V,BETA2,ZETAB,INFO)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-----------------------------------------------------------------------
C
C     F11GBS - STAGE 2: Initialization (Lanczos Method (SYMMLQ))
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA, BETA1, BETA2, BNORM2, RNORM2, ZETAB
      INTEGER           INFO, IREVCM, N
      LOGICAL           NEXT, PRECON
C     .. Array Arguments ..
      DOUBLE PRECISION  U(N), V(N), V1(N), V2(N), W1(N), W2(N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  T
      LOGICAL           DOV1
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2
      EXTERNAL          DDOT, DNRM2
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Save statement ..
      SAVE              DOV1
C     .. Executable Statements ..
C
C     First call to F11GBS
C
      IF (NEXT) THEN
         NEXT = .FALSE.
         IREVCM = 1
         DOV1 = .TRUE.
         BNORM2 = DNRM2(N,V1,1)
         CALL DCOPY(N,X,1,U,1)
C
C     Subsequent calls to F11GBS
C
      ELSE
C
C        Compute v[1], beta[1], w_bar[1] = z[1]
C
         IF (DOV1) THEN
            IF (IREVCM.LE.1) THEN
               CALL DAXPY(N,-ONE,V,1,V1,1)
               IF (PRECON) THEN
                  IREVCM = 2
                  CALL DCOPY(N,V1,1,U,1)
               ELSE
                  DOV1 = .FALSE.
                  CALL DCOPY(N,V1,1,W1,1)
               END IF
            ELSE
               DOV1 = .FALSE.
               CALL DCOPY(N,V,1,W1,1)
            END IF
            IF ( .NOT. DOV1) THEN
               BETA1 = DDOT(N,V1,1,W1,1)
               RNORM2 = DNRM2(N,V1,1)
               IF (BETA1.LE.ZERO) THEN
                  NEXT = .TRUE.
                  IF (BETA1.LT.ZERO) INFO = 6
               ELSE
                  IREVCM = 1
                  BETA1 = SQRT(BETA1)
                  CALL DSCAL(N,ONE/BETA1,W1,1)
                  CALL DCOPY(N,W1,1,U,1)
               END IF
            END IF
C
C        Compute alpha[1], v[2]; ensure that v[2]'*z[1] = 0,
C        compute beta[2], z[2]
C
         ELSE
            IF (IREVCM.LE.1) THEN
               CALL DCOPY(N,V,1,V2,1)
               ALPHA = DDOT(N,V2,1,W1,1)
               CALL DAXPY(N,-ALPHA/BETA1,V1,1,V2,1)
               T = -DDOT(N,V2,1,W1,1)/DDOT(N,W1,1,W1,1)
               CALL DAXPY(N,T,W1,1,V2,1)
               IF (PRECON) THEN
                  IREVCM = 2
                  CALL DCOPY(N,V2,1,U,1)
               ELSE
                  NEXT = .TRUE.
                  CALL DCOPY(N,V2,1,W2,1)
               END IF
            ELSE
               NEXT = .TRUE.
               CALL DCOPY(N,V,1,W2,1)
               IREVCM = 1
            END IF
            IF (NEXT) THEN
               BETA2 = DDOT(N,V2,1,W2,1)
               ZETAB = BETA1/ALPHA
               IF (BETA2.LT.ZERO) THEN
                  INFO = 6
               ELSE
                  BETA2 = SQRT(BETA2)
                  IF (BETA2.GT.ZERO) CALL DSCAL(N,ONE/BETA2,W2,1)
               END IF
            END IF
         END IF
C
      END IF
C
C     End of F11GBS
C
      RETURN
      END
