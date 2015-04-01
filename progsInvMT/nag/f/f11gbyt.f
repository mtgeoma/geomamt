      SUBROUTINE F11GBY(IREVCM,IDATA,RDATA,U,V,WORK,LWORK,INFO)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-----------------------------------------------------------------------
C
C     F11GBY - Lanczos Method (SYMMLQ)
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0)
      INTEGER           STBEG, STNORM, STINIT, STLOOP, STTERM, STMONI,
     *                  STEND
      PARAMETER         (STBEG=0,STNORM=1,STINIT=2,STLOOP=3,STTERM=4,
     *                  STMONI=5,STEND=6)
C     .. Scalar Arguments ..
      INTEGER           INFO, IREVCM, LWORK
C     .. Array Arguments ..
      DOUBLE PRECISION  RDATA(20), U(*), V(*), WORK(LWORK)
      INTEGER           IDATA(20)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA, ANORM, BETA1, BETA2, BETAIN, BNORM,
     *                  BNORM2, GAMMAB, PI, RHO1, RHO2, RNORM2, SIGERC,
     *                  SIGERR, SIGMAC, SIGMAX, STPLHS, STPRHS, TOL,
     *                  XNORM, XNORM0, ZETA, ZETAB
      INTEGER           IB, ID, IE2, INFOCH, ITERM, ITN, ITS, IV1, IV2,
     *                  IW1, IW2, IWGT, IX, KILL, MAXITN, MAXITS, MONIT,
     *                  N, NORM, RESID, SIGCMP, SIGCMX, STAGE
      LOGICAL           FLOOP, NEXT, PRECON, USEXC
C     .. External Functions ..
      DOUBLE PRECISION  F11BBU
      EXTERNAL          F11BBU
C     .. External Subroutines ..
      EXTERNAL          DCOPY, F06FBF, F11BBY, F11GBP, F11GBQ, F11GBR,
     *                  F11GBS, F11GBW
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MOD
C     .. Save statement ..
      SAVE              ALPHA, ANORM, BETA1, BETA2, BETAIN, BNORM,
     *                  BNORM2, GAMMAB, PI, RHO1, RHO2, RNORM2, SIGERR,
     *                  SIGERC, SIGMAX, SIGMAC, STPLHS, STPRHS, TOL,
     *                  XNORM, XNORM0, ZETA, ZETAB, IB, ID, IE2, INFOCH,
     *                  ITERM, ITN, ITS, IV1, IV2, IW1, IW2, IWGT, IX,
     *                  MAXITN, MAXITS, MONIT, N, NORM, RESID, SIGCMP,
     *                  SIGCMX, STAGE, FLOOP, NEXT, PRECON, USEXC
C     .. Executable Statements ..
C
C     Initialize
C
      KILL = 0
      IF (IREVCM.GE.100) THEN
         KILL = 2
         IREVCM = IREVCM/100
         STAGE = STEND
      ELSE IF (IREVCM.GE.10) THEN
         KILL = 1
         IREVCM = IREVCM/10
      ELSE IF (IREVCM.EQ.0) THEN
         STAGE = STBEG
      END IF
C
   20 CONTINUE
C-----------------------------------------------------------------------
C
C     STAGE 0: initialization
C
C-----------------------------------------------------------------------
      IF (STAGE.EQ.STBEG) THEN
         NEXT = .TRUE.
         RESID = 0
C
C        Get all parameters
C
         PRECON = (IDATA(3).GE.1)
         SIGCMP = IDATA(4)
         NORM = IDATA(5)
         ITERM = IDATA(6)
         N = IDATA(7)
         MAXITN = IDATA(8)
         MAXITS = IDATA(9)
         MONIT = IDATA(10)
         TOL = RDATA(1)
         ANORM = RDATA(3)
         SIGMAX = RDATA(4)
C
C        Assign storage
C
         IWGT = 1
         IF (NORM.GE.3) THEN
            IX = IWGT + N
         ELSE
            IX = IWGT
         END IF
         IB = IX + N
         IV1 = IB + N
         IV2 = IV1 + N
         IW1 = IV2 + N
         IW2 = IW1 + N
         IF (SIGCMP.GE.2) THEN
            ID = IW2 + N
            IE2 = ID + MAXITS + 1
         ELSE
            ID = IW1
            IE2 = IW2
         END IF
C
C        Copy vectors
C
         CALL DCOPY(N,U,1,WORK(IX),1)
         CALL DCOPY(N,V,1,WORK(IB),1)
C
C        Initialize other parameters
C
         INFOCH = 1
         STPLHS = ZERO
         STPRHS = ZERO
         ITN = 0
         ITS = 0
         SIGCMX = -SIGCMP
         SIGERR = ZERO
C
C        Compute the norm of the right-hand side b
C
         BNORM = F11BBU(NORM,N,WORK(IB),WORK(IWGT))
C
C        Assign next stage
C
         IF ((ITERM.LE.1) .AND. (ANORM.LE.ZERO)) THEN
            STAGE = STNORM
         ELSE
            STAGE = STINIT
         END IF
C-----------------------------------------------------------------------
C
C     STAGE 1: Compute the 1-norm of the matrix A
C
C-----------------------------------------------------------------------
      ELSE IF (STAGE.EQ.STNORM) THEN
         CALL F11BBY(NEXT,IREVCM,NORM,ANORM,N,WORK(IV2),U,V)
         IF (NEXT) THEN
            STAGE = STINIT
         ELSE
            IREVCM = ABS(IREVCM)
         END IF
C-----------------------------------------------------------------------
C
C     STAGE 2: Initialize the computation:
C              Generate v[1], w_bar[1], alpha[1], beta[1]
C                       v[2], z[2], beta[2]
C
C-----------------------------------------------------------------------
      ELSE IF (STAGE.EQ.STINIT) THEN
         IF (NEXT) THEN
            RESID = 0
            INFOCH = 1
            PI = ONE
            RHO2 = ZERO
            FLOOP = .TRUE.
            CALL DCOPY(N,WORK(IB),1,WORK(IV1),1)
            XNORM = F11BBU(NORM,N,WORK(IX),WORK(IWGT))
            XNORM0 = XNORM
         END IF
         CALL F11GBS(NEXT,IREVCM,PRECON,N,ALPHA,BETA1,BNORM2,RNORM2,
     *               WORK(IX),WORK(IV1),WORK(IV2),WORK(IW1),WORK(IW2),U,
     *               V,BETA2,ZETAB,INFO)
         IF (NEXT) THEN
            IF (INFO.LE.0) THEN
               IF (BETA1.LE.ZERO) THEN
                  STAGE = STEND
                  RESID = 2
                  CALL DCOPY(N,WORK(IX),1,U,1)
                  CALL F06FBF(N,ZERO,V,1)
                  IF (ITERM.LE.1) THEN
                     STPRHS = TOL*(BNORM+ANORM*XNORM)
                  ELSE
                     STPRHS = TOL*(BETAIN+SIGMAX*XNORM)
                  END IF
               ELSE
                  BETAIN = BETA1
                  RHO1 = BETA1
                  RHO2 = ZERO
                  GAMMAB = ALPHA
                  ZETA = ZERO
                  STAGE = STTERM
               END IF
            ELSE
               STAGE = STEND
            END IF
            IF ((SIGCMX.NE.0) .AND. (BETA1.GT.ZERO)) THEN
               ITS = MAXITS
               CALL F11GBW(SIGCMX,ITS,ALPHA,RDATA(2),WORK(ID),WORK(IE2),
     *                     SIGMAC,SIGERC)
               IF (ITN.LE.0) THEN
                  SIGMAX = SIGMAC
                  SIGERR = SIGERC
               ELSE IF (SIGCMP.LE.1) THEN
                  SIGMAX = MAX(SIGMAX,SIGMAC)
               ELSE IF (SIGERC.LT.SIGERR) THEN
                  SIGMAX = SIGMAC
                  SIGERR = SIGERC
               END IF
            END IF
         END IF
C-----------------------------------------------------------------------
C
C     STAGE 3: Body of the loop
C
C-----------------------------------------------------------------------
      ELSE IF (STAGE.EQ.STLOOP) THEN
         CALL F11GBR(NEXT,FLOOP,IREVCM,ITN,PRECON,N,ALPHA,BETA1,BETA2,
     *               PI,GAMMAB,RHO1,RHO2,ZETA,ZETAB,WORK(IX),WORK(IV1),
     *               WORK(IV2),WORK(IW1),WORK(IW2),U,V,INFO)
         IF (NEXT) THEN
            IF (INFO.EQ.0) THEN
               IF (SIGCMX.GT.0) THEN
                  CALL F11GBW(SIGCMX,ITS,ALPHA,BETA1,WORK(ID),WORK(IE2),
     *                        SIGMAC,SIGERC)
                  IF (SIGCMP.LE.1) THEN
                     SIGMAX = MAX(SIGMAX,SIGMAC)
                  ELSE IF (SIGERC.LE.SIGERR) THEN
                     SIGMAX = SIGMAC
                     SIGERR = SIGERC
                  END IF
               END IF
               STAGE = STTERM
            ELSE
               STAGE = STEND
            END IF
         END IF
C-----------------------------------------------------------------------
C
C     STAGE 4: Check the termination condition
C
C-----------------------------------------------------------------------
      ELSE IF (STAGE.EQ.STTERM) THEN
         CALL F11GBQ(NEXT,IREVCM,ITN,RESID,PRECON,NORM,ITERM,N,TOL,
     *               SIGMAX,ANORM,XNORM0,BNORM,BNORM2,RNORM2,BETAIN,
     *               BETA1,BETA2,PI,GAMMAB,RHO1,RHO2,ZETA,ZETAB,WORK(IB)
     *               ,WORK(IX),WORK(IV1),WORK(IV2),WORK(IW1),WORK(IWGT),
     *               U,V,STPLHS,STPRHS,USEXC,INFOCH)
         IF (NEXT) THEN
            IF (INFOCH.EQ.0) THEN
               INFO = 0
               STAGE = STEND
            ELSE IF (INFOCH.LE.2) THEN
               IF (INFO.EQ.0) THEN
                  STAGE = STLOOP
                  IF (KILL.GE.1) THEN
                     IF (INFOCH.LE.1) THEN
                        INFO = 4
                     ELSE
                        INFO = 2
                     END IF
                     STAGE = STEND
                  ELSE IF (ITN.GE.MAXITN) THEN
                     INFO = 5
                     STAGE = STEND
                  ELSE IF ((ITN.GT.0) .AND. (MONIT.GT.0)) THEN
                     IF (MOD(ITN,MONIT).EQ.0) STAGE = STMONI
                  ELSE IF (INFOCH.LT.0) THEN
                     STAGE = STINIT
                  END IF
               ELSE
                  STAGE = STEND
               END IF
            ELSE IF (INFOCH.EQ.3) THEN
               STAGE = STTERM
            ELSE IF (INFOCH.EQ.4) THEN
               INFO = 2
               STAGE = STEND
            END IF
         END IF
C-----------------------------------------------------------------------
C
C     STAGE 5: Monitoring
C
C-----------------------------------------------------------------------
      ELSE IF (STAGE.EQ.STMONI) THEN
         IREVCM = 3
         CALL F11GBP(RESID,USEXC,ITERM,N,ZETAB,-(BETAIN*PI/GAMMAB),
     *               BETA1,BETA2,RHO1,RHO2,WORK(IX),WORK(IV1),WORK(IV2),
     *               WORK(IW1),U,V)
         IDATA(13) = ITN
         IDATA(14) = ITS
         RDATA(3) = ANORM
         RDATA(4) = SIGMAX
         RDATA(5) = STPLHS
         RDATA(6) = STPRHS
         RDATA(7) = SIGERR
         IF (INFOCH.GE.0) THEN
            STAGE = STLOOP
         ELSE
            STAGE = STINIT
         END IF
C-----------------------------------------------------------------------
C
C     STAGE 6: Termination
C
C-----------------------------------------------------------------------
      ELSE IF (STAGE.EQ.STEND) THEN
         IREVCM = 4
         IF (KILL.LE.1) THEN
            IF (INFO.LE.5) THEN
               IF ((ITN.LE.0) .AND. (BETAIN.EQ.ZERO)) THEN
                  CALL DCOPY(N,WORK(IX),1,U,1)
               ELSE
                  CALL F11GBP(RESID,USEXC,ITERM,N,ZETAB,
     *                        -(BETAIN*PI/GAMMAB),BETA1,BETA2,RHO1,RHO2,
     *                        WORK(IX),WORK(IV1),WORK(IV2),WORK(IW1),U,
     *                        V)
                  IF ((RESID.LE.0) .AND. (ITERM.LE.1)) THEN
                     XNORM = F11BBU(NORM,N,U,WORK(IWGT))
                     STPLHS = F11BBU(NORM,N,V,WORK(IWGT))
                     STPRHS = TOL*(BNORM+ANORM*MAX(XNORM,XNORM0))
                  END IF
               END IF
            ELSE
               CALL DCOPY(N,WORK(IX),1,U,1)
               CALL F06FBF(N,ZERO,V,1)
            END IF
         ELSE
            CALL DCOPY(N,WORK(IX),1,U,1)
            CALL F06FBF(N,ZERO,V,1)
         END IF
         IDATA(13) = ITN
         IDATA(14) = ITS
         RDATA(3) = ANORM
         RDATA(4) = SIGMAX
         RDATA(5) = STPLHS
         RDATA(6) = STPRHS
         RDATA(7) = SIGERR
C
      END IF
C
C     Loop over the stages
C
      IF (NEXT .AND. (IREVCM.LE.2)) GO TO 20
C
      IDATA(12) = IREVCM
C
C     End of subroutine F11GBY
C
      RETURN
      END
