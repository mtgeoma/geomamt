      SUBROUTINE F11GBZ(IREVCM,IDATA,RDATA,U,V,WORK,LWORK,INFO)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-----------------------------------------------------------------------
C
C     F11GBZ - Preconditioned Conjugate Gradient
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
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
      DOUBLE PRECISION  ANORM, BNORM, SIGERC, SIGERR, SIGMAC, SIGMAX,
     *                  STPLHS, STPRHS, TALPHA, TBETA, TOL, XNORM
      INTEGER           IB, ID, IE2, INFOCH, IP, IR, ITERM, ITN, ITS,
     *                  IW, IWGT, IX, KILL, MAXITN, MAXITS, MONIT, N,
     *                  NORM, SIGCMP, SIGCMX, STAGE
      LOGICAL           FLOOP, NEXT, PRECON
C     .. External Subroutines ..
      EXTERNAL          DCOPY, F11BBY, F11GBU, F11GBV, F11GBW, F11GBX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MOD
C     .. Save statement ..
      SAVE              ANORM, BNORM, SIGERR, SIGERC, SIGMAX, SIGMAC,
     *                  STPLHS, STPRHS, TALPHA, TBETA, TOL, XNORM, IB,
     *                  ID, IE2, INFOCH, IP, IR, ITERM, ITN, ITS, IW,
     *                  IWGT, IX, MAXITN, MAXITS, MONIT, N, NORM,
     *                  SIGCMP, SIGCMX, STAGE, FLOOP, NEXT, PRECON
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
         IR = IB + N
         IP = IR + N
         IW = IP + N
         IF (SIGCMP.GE.2) THEN
            ID = IW + N
            IE2 = ID + MAXITS + 1
         ELSE
            ID = IP
            IE2 = IW
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
C        Assign next stage
C
         IF ((ABS(ITERM).LE.1) .AND. (ANORM.LE.0)) THEN
            STAGE = STNORM
         ELSE
            STAGE = STINIT
         END IF
C-----------------------------------------------------------------------
C
C     STAGE 1: Compute the norm of the matrix A
C
C-----------------------------------------------------------------------
      ELSE IF (STAGE.EQ.STNORM) THEN
         CALL F11BBY(NEXT,IREVCM,NORM,ANORM,N,WORK(IW),U,V)
         IF (NEXT) THEN
            STAGE = STINIT
         ELSE
            IREVCM = ABS(IREVCM)
         END IF
C-----------------------------------------------------------------------
C
C     STAGE 2: Start the computation:
C              Compute r[0] = (b - A*x[0]),  z[0] = inv(M)*r[0],
C              rho[0] = r[0]'*z[0] (stored in TALPHA),  p[1] = z[0]
C
C-----------------------------------------------------------------------
      ELSE IF (STAGE.EQ.STINIT) THEN
         IF (NEXT .AND. (INFOCH.LT.0) .AND. (SIGCMP.EQ.2)
     *       .AND. (SIGCMX.NE.0)) SIGCMX = -2
         CALL F11GBX(NEXT,IREVCM,PRECON,NORM,ITERM,BNORM,XNORM,N,
     *               WORK(IB),WORK(IX),WORK(IR),WORK(IP),WORK(IW),
     *               WORK(IWGT),U,V,TALPHA,TBETA,STPLHS,INFOCH)
         IF (NEXT) THEN
            IF (INFOCH.EQ.0) THEN
               STAGE = STTERM
            ELSE IF (INFOCH.GT.1) THEN
               INFO = INFOCH
               STAGE = STEND
            ELSE
               INFOCH = 1
               FLOOP = .TRUE.
               IF (SIGCMX.NE.0) THEN
                  ITS = MAXITS
                  CALL F11GBW(SIGCMX,ITS,TBETA/TALPHA,RDATA(2),WORK(ID),
     *                        WORK(IE2),SIGMAC,SIGERC)
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
               STAGE = STTERM
            END IF
         END IF
C-----------------------------------------------------------------------
C
C     STAGE 3: Body of the loop
C
C-----------------------------------------------------------------------
      ELSE IF (STAGE.EQ.STLOOP) THEN
         CALL F11GBV(NEXT,FLOOP,IREVCM,ITN,MAXITN,ITERM,PRECON,NORM,N,
     *               WORK(IX),WORK(IR),WORK(IP),WORK(IW),WORK(IWGT),U,V,
     *               TALPHA,TBETA,XNORM,STPLHS,INFO)
         IF (NEXT) THEN
            IF ((INFO.EQ.0) .AND. (SIGCMX.GT.0)) THEN
               CALL F11GBW(SIGCMX,ITS,TALPHA,TBETA,WORK(ID),WORK(IE2),
     *                     SIGMAC,SIGERC)
               IF (SIGCMP.LE.1) THEN
                  SIGMAX = MAX(SIGMAX,SIGMAC)
               ELSE IF (SIGERC.LE.SIGERR) THEN
                  SIGMAX = SIGMAC
                  SIGERR = SIGERC
               END IF
            END IF
            STAGE = STTERM
         END IF
C-----------------------------------------------------------------------
C
C     STAGE 4: Check the termination condition
C
C-----------------------------------------------------------------------
      ELSE IF (STAGE.EQ.STTERM) THEN
         CALL F11GBU(NEXT,IREVCM,ITN,PRECON,NORM,ITERM,TOL,ANORM,BNORM,
     *               XNORM,SIGMAX,N,WORK(IB),WORK(IX),WORK(IWGT),U,V,
     *               STPLHS,STPRHS,INFOCH)
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
                     STAGE = STEND
                     INFO = 5
                  ELSE IF (( .NOT. FLOOP) .AND. (MONIT.GT.0)) THEN
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
         CALL DCOPY(N,WORK(IX),1,U,1)
         IF ((ABS(INFOCH).NE.2) .AND. (INFOCH.NE.4)) CALL DCOPY(N,
     *       WORK(IR),1,V,1)
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
         CALL DCOPY(N,WORK(IX),1,U,1)
         IF ((INFOCH.EQ.1) .OR. (KILL.GT.1)) CALL DCOPY(N,WORK(IR),1,V,
     *       1)
         IDATA(13) = ITN
         IDATA(14) = ITS
         RDATA(3) = ANORM
         RDATA(4) = SIGMAX
         RDATA(5) = STPLHS
         RDATA(6) = STPRHS
         RDATA(7) = SIGERR
      END IF
C
C     Loop over the stages
C
      IF (NEXT .AND. (IREVCM.LE.2)) GO TO 20
C
      IDATA(12) = IREVCM
C
C     End of subroutine F11GBZ
C
      RETURN
      END
