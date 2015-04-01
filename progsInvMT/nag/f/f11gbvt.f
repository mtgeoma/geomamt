      SUBROUTINE F11GBV(NEXT,FLOOP,IREVCM,ITN,MAXITN,ITERM,PRECON,NORM,
     *                  N,X,R,P,W,WGT,U,V,TALPHA,TBETA,XNORM,STPLHS,
     *                  INFO)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-----------------------------------------------------------------------
C
C     F11GBV - Body of the loop (Conjugate Gradient Method)
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  STPLHS, TALPHA, TBETA, XNORM
      INTEGER           INFO, IREVCM, ITERM, ITN, MAXITN, N, NORM
      LOGICAL           FLOOP, NEXT, PRECON
C     .. Array Arguments ..
      DOUBLE PRECISION  P(N), R(N), U(N), V(N), W(N), WGT(N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA, BETA, DELTA, RHO, RHOX, XI, XIX
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, F11BBU
      EXTERNAL          DDOT, F11BBU
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Save statement ..
      SAVE              RHO, XI
C     .. Executable Statements ..
C
C     First call to this routine
C
      IF (NEXT) THEN
         IREVCM = 1
         ITN = ITN + 1
         IF (FLOOP) THEN
            FLOOP = .FALSE.
            RHO = TALPHA
            XI = TBETA
         END IF
         ALPHA = RHO/XI
         CALL DAXPY(N,ALPHA,P,1,X,1)
         CALL DAXPY(N,-ALPHA,W,1,R,1)
         IF (ITN.LE.MAXITN) THEN
            NEXT = .FALSE.
            IF (PRECON) THEN
               IREVCM = 2
            ELSE
               IREVCM = 1
            END IF
            CALL DCOPY(N,R,1,U,1)
         ELSE IF (PRECON) THEN
            NEXT = .FALSE.
            IREVCM = 2
            CALL DCOPY(N,R,1,U,1)
         ELSE
            NEXT = .TRUE.
            CALL DCOPY(N,R,1,U,1)
         END IF
C
C     Subsequent calls to this routine
C
      ELSE
         IF (IREVCM.LE.1) THEN
            NEXT = .TRUE.
         ELSE
            IF (ITN.LE.MAXITN) THEN
               IREVCM = 1
               CALL DCOPY(N,V,1,U,1)
            ELSE
               NEXT = .TRUE.
            END IF
         END IF
      END IF
C
C     Completion
C
      IF (NEXT) THEN
         RHOX = RHO
         XIX = XI
         RHO = DDOT(N,R,1,U,1)
         DELTA = DDOT(N,U,1,V,1)
         XNORM = F11BBU(NORM,N,X,WGT)
         IF (ITERM.LE.1) THEN
            STPLHS = F11BBU(NORM,N,R,WGT)
         ELSE
            STPLHS = F11BBU(NORM,N,U,WGT)
         END IF
         IF (ITN.LE.MAXITN) THEN
            IF (RHO.LE.ZERO) THEN
               IF (STPLHS.GT.ZERO) INFO = 6
            ELSE
               BETA = RHO/RHOX
               XI = DELTA - BETA*BETA*XIX
               IF (XI.EQ.ZERO) THEN
                  IF (STPLHS.GT.ZERO) INFO = 7
               ELSE
                  TALPHA = (BETA**2*XIX+XI)/RHO
                  TBETA = -XIX*BETA/(RHO*SQRT(RHOX/RHO))
                  CALL DSCAL(N,BETA,P,1)
                  CALL DAXPY(N,ONE,U,1,P,1)
                  CALL DSCAL(N,BETA,W,1)
                  CALL DAXPY(N,ONE,V,1,W,1)
               END IF
            END IF
         END IF
      END IF
C
C     End of subroutine F11GBV
C
      RETURN
      END
