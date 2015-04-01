      SUBROUTINE D02PDV(F,HAVG,JFLSTP,TOOMCH,MAXFCN,WORK,IER,NREC)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C$$$$ SUBROUTINE STIFF $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C************************************************
C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
C************************************************
C
C  Purpose:      Diagnose stiffness. This depends on two things: whether
C                the step size is being restricted on grounds of
C                stability and whether the integration to TND can be
C                completed in no more than MAXFCN function evaluations.
C
C  Input:        HAVG, TOOMCH, MAXFCN, WORK(*)
C  Input/output: JFLSTP
C  Output:       IER, NREC
C  Workspace:    WORK(*)
C  External:     F
C
C  Common:       Initializes:   /JD02PD/ REC
C                Reads:         /AD02PD/ TND, NEQN
C                               /BD02PD/ T, H, NFCN, SVNFCN, OKSTP
C                               /CD02PD/ PRY, PRYP, PRTHRS, PRWT, PRSCR,
C                                        PRSTGS, PRYOLD
C                               /ED02PD/ COST
C                Alters:        /BD02PD/ NFCN
C                               /JD02PD/ REC
C
C     .. Parameters ..
      DOUBLE PRECISION  HALF
      PARAMETER         (HALF=0.5D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HAVG
      INTEGER           IER, JFLSTP, MAXFCN, NREC
      LOGICAL           TOOMCH
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(*)
C     .. Subroutine Arguments ..
      EXTERNAL          F
C     .. Scalars in Common ..
      DOUBLE PRECISION  COST, DIR, EXPON, H, HOLD, HSTRT, RS, RS1, RS2,
     *                  RS3, RS4, SAFETY, STBRAD, T, TANANG, TND, TOLD,
     *                  TOLR, TOOSML, TSTRT
      INTEGER           FLSTP, LNINTP, LSTSTG, MAXTRY, NEQN, NFCN, NSEC,
     *                  OKSTP, ORDER, PRERST, PRINTP, PRSCR, PRSTGS,
     *                  PRTHRS, PRWT, PRY, PRYOLD, PRYP, SVNFCN
      LOGICAL           FIRST, FSAL, LAST
C     .. Arrays in Common ..
      CHARACTER*80      REC(10)
C     .. Local Scalars ..
      DOUBLE PRECISION  AVGY, XTRAWK
      INTEGER           L
      LOGICAL           LOTSFL, STIF, UNSURE
C     .. External Subroutines ..
      EXTERNAL          D02PDU
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX, MOD
C     .. Common blocks ..
      COMMON            /AD02PD/TSTRT, TND, DIR, HSTRT, TOLR, NEQN
      COMMON            /BD02PD/T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP,
     *                  FLSTP, FIRST, LAST
      COMMON            /CD02PD/PRTHRS, PRERST, PRWT, PRYOLD, PRSCR,
     *                  PRY, PRYP, PRSTGS, PRINTP, LNINTP
      COMMON            /ED02PD/TOOSML, COST, SAFETY, EXPON, STBRAD,
     *                  TANANG, RS, RS1, RS2, RS3, RS4, ORDER, LSTSTG,
     *                  MAXTRY, NSEC, FSAL
      COMMON            /JD02PD/REC
C     .. Save statement ..
      SAVE              /AD02PD/, /BD02PD/, /CD02PD/, /ED02PD/, /JD02PD/
C     .. Executable Statements ..
C
      IF (MOD(OKSTP-10,40).EQ.0) THEN
         LOTSFL = JFLSTP .GE. 10
         JFLSTP = 0
      ELSE
         LOTSFL = .FALSE.
      END IF
C
C  If either too much work has been done or there are lots of failed
C  steps, test for stiffness.
C
      IF (TOOMCH .OR. LOTSFL) THEN
C
C  Regenerate weight vector
C
         DO 20 L = 1, NEQN
            AVGY = HALF*(ABS(WORK(PRY-1+L))+ABS(WORK(PRYOLD-1+L)))
            WORK(PRWT-1+L) = MAX(AVGY,WORK(PRTHRS-1+L))
   20    CONTINUE
C
C  D02PDU determines whether the problem is STIFF. In some circumstances
C  it is UNSURE. The decision depends on two things: whether the step
C  size is being restricted on grounds of stability and whether the
C  integration to TND can be completed in no more than MAXFCN function
C  evaluations. The last four arguments of D02PDU are vectors of length
C  NEQN used for working storage. Some storage in WORK(*) reserved for
C  the stages (there are a minimum of three such vectors reserved for
C  the METHDs implemented) and the scratch vector starting at PRSCR are
C  used for this purpose.
C
         CALL D02PDU(F,T,WORK(PRY),H,HAVG,TND,MAXFCN,WORK(PRWT),
     *               WORK(PRYP),WORK(PRERST),UNSURE,STIF,WORK(PRSTGS),
     *               WORK(PRSTGS+NEQN),WORK(PRSTGS+2*NEQN),WORK(PRSCR))
         IF ( .NOT. UNSURE) THEN
            IF (STIF) THEN
C
C  Predict how much eXTRA WorK will be needed to reach TND.
C
               XTRAWK = (COST*ABS((TND-T)/HAVG))/DBLE(SVNFCN+NFCN)
               IER = 4
               WRITE (REC(NREC+1),FMT='(A)')
     *          ' ** Your problem has been diagnosed as stiff.  If the '
               WRITE (REC(NREC+2),FMT='(A,D13.5)')
     *           ' ** situation persists, it will cost roughly ', XTRAWK
               WRITE (REC(NREC+3),FMT='(A)')
     *   ' ** times as much to reach TEND as it has cost to reach TNOW.'
               WRITE (REC(NREC+4),FMT='(A)')
     *          ' ** You should probably change to a code intended for '
               WRITE (REC(NREC+5),FMT='(A)') ' ** stiff problems. '
               NREC = NREC + 5
            END IF
         END IF
      END IF
C
      RETURN
      END
