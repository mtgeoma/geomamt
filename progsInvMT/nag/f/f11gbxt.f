      SUBROUTINE F11GBX(NEXT,IREVCM,PRECON,NORM,ITERM,BNORM,XNORM,N,B,X,
     *                  R,P,W,WGT,U,V,TALPHA,TBETA,STPLHS,INFOCH)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-----------------------------------------------------------------------
C
C     F11GBX - STAGE 2: Initialization (Conjugate Gradient Method)
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BNORM, STPLHS, TALPHA, TBETA, XNORM
      INTEGER           INFOCH, IREVCM, ITERM, N, NORM
      LOGICAL           NEXT, PRECON
C     .. Array Arguments ..
      DOUBLE PRECISION  B(N), P(N), R(N), U(N), V(N), W(N), WGT(N), X(N)
C     .. Local Scalars ..
      LOGICAL           FPASS
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, F11BBU
      EXTERNAL          DDOT, F11BBU
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY
C     .. Save statement ..
      SAVE              FPASS
C     .. Executable Statements ..
C
C     First calls to this routine
C
      IF (NEXT) THEN
         NEXT = .FALSE.
         IF ((INFOCH.EQ.1) .OR. (INFOCH.EQ.-1)) THEN
            FPASS = .TRUE.
            IREVCM = 1
            CALL DCOPY(N,B,1,R,1)
            CALL DCOPY(N,X,1,U,1)
         ELSE
            FPASS = .FALSE.
            CALL DCOPY(N,V,1,U,1)
            CALL DCOPY(N,U,1,R,1)
            IF (PRECON) THEN
               IREVCM = 2
            ELSE
               IREVCM = 1
            END IF
         END IF
C
C     Subsequent calls to this routine
C
      ELSE
         IF (FPASS) THEN
            FPASS = .FALSE.
            CALL DAXPY(N,-ONE,V,1,R,1)
            CALL DCOPY(N,R,1,U,1)
            IF (PRECON) THEN
               IREVCM = 2
            ELSE
               IREVCM = 1
            END IF
         ELSE
            IF (IREVCM.NE.1) THEN
               IREVCM = 1
               CALL DCOPY(N,V,1,U,1)
            ELSE
               NEXT = .TRUE.
               CALL DCOPY(N,U,1,P,1)
               CALL DCOPY(N,V,1,W,1)
               TALPHA = DDOT(N,R,1,P,1)
               TBETA = DDOT(N,P,1,W,1)
               XNORM = F11BBU(NORM,N,X,WGT)
               IF (INFOCH.GT.0) BNORM = F11BBU(NORM,N,B,WGT)
               IF (ITERM.LE.1) THEN
                  STPLHS = F11BBU(NORM,N,R,WGT)
               ELSE
                  STPLHS = F11BBU(NORM,N,P,WGT)
               END IF
               IF (TALPHA.LT.ZERO) THEN
                  INFOCH = 6
               ELSE IF (TALPHA.EQ.ZERO) THEN
                  INFOCH = 0
               ELSE IF (TBETA.EQ.ZERO) THEN
                  INFOCH = 7
               END IF
            END IF
         END IF
      END IF
C
C     End of subroutine F11GBX
C
      RETURN
      END
