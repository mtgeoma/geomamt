      SUBROUTINE D02PCF(F,TWANT,TGOT,YGOT,YPGOT,YMAX,WORK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C$$$$ SUBROUTINE UT $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C  If you are not familiar with the code D02PCF and how it is used in
C  conjunction with D02PVF to solve initial value problems, you should
C  study the document file RKSUITE.DOC carefully before proceeding
C  further. The following "Brief Reminder" is intended only to remind
C  you of the meaning, type, and size requirements of the arguments.
C
C  NAME DECLARED IN AN EXTERNAL STATEMENT IN THE CALLING PROGRAM:
C
C     F         - name of the subroutine for evaluating the differential
C                 equations.
C
C  The subroutine F must have the form
C
C  SUBROUTINE F(T,Y,YP)
C  DOUBLE PRECISION T,Y(*),YP(*)
C    Given input values of the independent variable T and the solution
C    components Y(*), for each L = 1,2,...,NEQ evaluate the differential
C    equation for the derivative of the Ith solution component and place
C    the value in YP(L).  Do not alter the input values of T and Y(*).
C  RETURN
C  END
C
C  INPUT VARIABLE
C
C     TWANT     - DOUBLE PRECISION
C                 The next value of the independent variable where a
C                 solution is desired.
C
C                 Constraints: TWANT must lie between the previous value
C                 of TGOT (TSTART on the first call) and TEND. TWANT can
C                 be equal to TEND, but it must be clearly
C                 distinguishable from the previous value of TGOT
C                 (TSTART on the first call) in the precision available.
C
C  OUTPUT VARIABLES
C
C     TGOT      - DOUBLE PRECISION
C                 A solution has been computed at this value of the
C                 independent variable.
C     YGOT(*)   - DOUBLE PRECISION array of length NEQ
C                 Approximation to the true solution at TGOT. Do not
C                 alter the contents of this array
C     YPGOT(*)  - DOUBLE PRECISION array of length NEQ
C                 Approximation to the first derivative of the true
C                 solution at TGOT.
C     YMAX(*)   - DOUBLE PRECISION array of length NEQ
C                 YMAX(L) is the largest magnitude of YGOT(L) computed
C                 at any time in the integration from TSTART to TGOT.
C                 Do not alter the contents of this array.
C
C  WORKSPACE
C
C     WORK(*)   - DOUBLE PRECISION array as used in D02PVF
C                 Do not alter the contents of this array.
C
C  OUTPUT VARIABLE
C
C     IFAIL     - INTEGER
C
C                       SUCCESS.  TGOT = TWANT.
C                 = 0 - Complete success.
C
C                       "SOFT" FAILURES
C                 = 2 - Warning:  You are using METHOD = 3 inefficiently
C                       by computing answers at many values of TWANT. If
C                       you really need answers at so many specific
C                       points, it would be more efficient to compute
C                       them with METHOD = 2.  To do this you would need
C                       to restart from TGOT, YGOT(*) by a call to
C                       D02PVF.  If you wish to continue as you are, you
C                       may.
C                 = 3 - Warning:  A considerable amount of work has been
C                       expended.  If you wish to continue on to TWANT,
C                       just call D02PCF again.
C                 = 4 - Warning:  It appears that this problem is
C                       "stiff". You really should change to another
C                       code that is intended for such problems, but if
C                       you insist, you can continue with D02PCF by
C                       calling it again.
C
C                       "HARD" FAILURES
C                 = 5 - You are asking for too much accuracy. You cannot
C                       continue integrating this problem.
C                 = 6 - The global error assessment may not be reliable
C                       beyond the current point in the integration. You
C                       cannot continue integrating this problem.
C
C                       "CATASTROPHIC" FAILURES
C                 = 1 - The nature of the catastrophe is reported on
C                       the standard output channel. Unless special
C                       provision was made in advance (see RKSUITE.DOC),
C                       the computation then comes to a STOP.
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02PCF')
      LOGICAL           ASK, TELL
      PARAMETER         (ASK=.TRUE.,TELL=.FALSE.)
      INTEGER           MINUS1, MINUS2
      PARAMETER         (MINUS1=-1,MINUS2=-2)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TGOT, TWANT
      INTEGER           IFAIL
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(*), YGOT(*), YMAX(*), YPGOT(*)
C     .. Subroutine Arguments ..
      EXTERNAL          F
C     .. Scalars in Common ..
      DOUBLE PRECISION  COST, CUBRMC, DIR, DWARF, EXPON, H, HOLD, HSTRT,
     *                  MCHEPS, RNDOFF, RS, RS1, RS2, RS3, RS4, SAFETY,
     *                  SQRRMC, STBRAD, T, TANANG, TINY, TND, TOLD,
     *                  TOLR, TOOSML, TSTRT
      INTEGER           FLSTP, LNINTP, LSTSTG, MAXTRY, METHD, MINTP,
     *                  NEQN, NFCN, NSEC, NSTAGE, OKSTP, ORDER, OUTCH,
     *                  PRERST, PRINTP, PRSCR, PRSTGS, PRTHRS, PRWT,
     *                  PRY, PRYOLD, PRYP, SVNFCN
      LOGICAL           FIRST, FSAL, INTP, LAST, UTASK
C     .. Arrays in Common ..
      DOUBLE PRECISION  A(13,13), B(13), BHAT(13), C(13), E(7), R(11,6)
      INTEGER           PTR(13)
      CHARACTER*80      REC(10)
C     .. Local Scalars ..
      DOUBLE PRECISION  HMIN, TLAST, TNOW, UTEND
      INTEGER           IER, JFAIL, L, NREC, STATE
      LOGICAL           BADERR, GOBACK
C     .. External Subroutines ..
      EXTERNAL          D02PDF, D02PDM, D02PDP, D02PWF, D02PXF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Common blocks ..
      COMMON            /AD02PD/TSTRT, TND, DIR, HSTRT, TOLR, NEQN
      COMMON            /BD02PD/T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP,
     *                  FLSTP, FIRST, LAST
      COMMON            /CD02PD/PRTHRS, PRERST, PRWT, PRYOLD, PRSCR,
     *                  PRY, PRYP, PRSTGS, PRINTP, LNINTP
      COMMON            /DD02PD/A, B, C, BHAT, R, E, PTR, NSTAGE, METHD,
     *                  MINTP, INTP
      COMMON            /ED02PD/TOOSML, COST, SAFETY, EXPON, STBRAD,
     *                  TANANG, RS, RS1, RS2, RS3, RS4, ORDER, LSTSTG,
     *                  MAXTRY, NSEC, FSAL
      COMMON            /GD02PD/MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC,
     *                  TINY, OUTCH
      COMMON            /HD02PD/UTASK
      COMMON            /JD02PD/REC
C     .. Save statement ..
      SAVE              /AD02PD/, /BD02PD/, /CD02PD/, /DD02PD/,
     *                  /ED02PD/, /GD02PD/, /HD02PD/, /JD02PD/, UTEND,
     *                  TLAST
C     .. Executable Statements ..
      IER = 0
      NREC = 0
      GOBACK = .FALSE.
      BADERR = .FALSE.
C
C  Is it permissible to call D02PCF?
C
      CALL D02PDM(ASK,'D02PVF',STATE)
      IF (STATE.EQ.1) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT='(A)')
     *   ' ** A catastrophic error has already been detected elsewhere.'
         GO TO 100
      END IF
      IF (STATE.EQ.MINUS1) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT='(A)')
     *     ' ** You have not called D02PVF, so you cannot use D02PCF.'
         GO TO 100
      END IF
      IF ( .NOT. UTASK) THEN
         IER = 1
         NREC = 2
         WRITE (REC,FMT='(A/A)')
     *  ' ** You have called D02PCF after you specified in D02PVF that '
     *     , ' ** you were going to use D02PDF. This is not permitted.'
         GO TO 100
      END IF
      CALL D02PDM(ASK,SRNAME,STATE)
      IF (STATE.EQ.5 .OR. STATE.EQ.6) THEN
         IER = 1
         NREC = 3
         WRITE (REC,FMT='(A/A/A)')
     *     ' ** D02PCF has already returned with IFAIL set to 5 or 6.',
     *     ' ** You cannot continue integrating this problem. You must '
     *     , ' ** call D02PVF to start another problem.'
         GO TO 100
      END IF
      STATE = MINUS2
      CALL D02PDM(TELL,SRNAME,STATE)
C
      IF (FIRST) THEN
C
C  First call.
C
C  A value of TND is specified in D02PVF. When INTP = .FALSE., as with
C  METHD = 3, output is obtained at the specified TWANT by resetting TND
C  to TWANT.  At this point, before the integration gets started, this
C  can be done with a simple assignment.  Later it is done with a call
C  to D02PWF. The original TND is SAVEd as a local variable UTEND.
C
         UTEND = TND
         IF ( .NOT. INTP) TND = TWANT
C
C  The last TGOT returned is SAVEd in the variable TLAST.  T (a variable
C  passed through the common block BD02PD) records how far the
C  integration has advanced towards the specified TND.  When output is
C  obtained by interpolation, the integration goes past the TGOT
C  returned (T is closer to the specified TND than TGOT).  Initialize
C  these variables and YMAX(*).
C
         TLAST = TSTRT
         TGOT = TSTRT
         DO 20 L = 1, NEQN
            YMAX(L) = ABS(WORK(PRY-1+L))
   20    CONTINUE
C
C  If the code is to find an on-scale initial step size H, a bound was
C  placed on H in D02PVF.  Here the first output point is used to
C  refine this bound.
C
         IF (HSTRT.EQ.ZERO) THEN
            H = MIN(ABS(H),ABS(TWANT-TSTRT))
            HMIN = MAX(TINY,TOOSML*MAX(ABS(TSTRT),ABS(TND)))
            H = MAX(H,HMIN)
         END IF
C
      ELSE
C
C  Subsequent call.
C
         IF (TLAST.EQ.UTEND) THEN
            IER = 1
            NREC = 3
            WRITE (REC,FMT='(A/A/A)')
     *  ' ** You have called D02PCF after reaching TEND. (Your last    '
     *        ,
     *   ' ** call to D02PCF resulted in TGOT = TEND.)  To start a new '
     *        , ' ** problem, you will need to call D02PVF.'
            GO TO 100
         END IF
C
      END IF
C
C  Check for valid TWANT.
C
      IF (DIR*(TWANT-TLAST).LE.ZERO) THEN
         IER = 1
         NREC = 4
         WRITE (REC,FMT='(A/A/A/A)')
     *    ' ** You have made a call to D02PCF with a TWANT that does   '
     *     , ' ** not lie between the previous value of TGOT (TSTART  ',
     *     ' ** on the first call) and TEND. This is not permitted. ',
     *     ' ** Check your program carefully.'
         GO TO 100
      END IF
      IF (DIR*(TWANT-UTEND).GT.ZERO) THEN
         HMIN = MAX(TINY,TOOSML*MAX(ABS(TWANT),ABS(UTEND)))
         IF (ABS(TWANT-UTEND).LT.HMIN) THEN
            IER = 1
            NREC = 5
            WRITE (REC,FMT='(A/A/A/A)')
     * ' ** You have made a call to D02PCF with a TWANT that does      '
     *        ,
     *     ' ** not lie between the previous value of TGOT (TSTART on  '
     *        ,
     *     ' ** the first call) and TEND. This is not permitted. TWANT '
     *        ,
     *     ' ** is very close to TEND, so you may have meant to set    '
     *        ,
     *     ' ** it to be TEND exactly.  Check your program carefully.  '
         ELSE
            IER = 1
            NREC = 4
            WRITE (REC,FMT='(A/A/A/A)')
     *    ' ** You have made a call to D02PCF with a TWANT that does   '
     *        ,
     *        ' ** not lie between the previous value of TGOT (TSTART  '
     *        ,
     *        ' ** on the first call) and TEND. This is not permitted. '
     *        , ' ** Check your program carefully.'
         END IF
         GO TO 100
      END IF
      IF ( .NOT. INTP) THEN
         HMIN = MAX(TINY,TOOSML*MAX(ABS(TLAST),ABS(TWANT)))
         IF (ABS(TWANT-TLAST).LT.HMIN) THEN
            IER = 1
            NREC = 4
            WRITE (REC,FMT='(A/A/A/A,D13.5,A)')
     *    ' ** You have made a call to D02PCF with a TWANT that is not '
     *        ,
     *        ' ** sufficiently different from the last value of TGOT  '
     *        ,
     *        ' ** (TSTART on the first call).  When using METHOD = 3, '
     *        , ' ** it must differ by at least ', HMIN, '.'
            GO TO 100
         END IF
C
C  We have a valid TWANT. There is no interpolation with this METHD and
C  therefore we step to TWANT exactly by resetting TND with a call to
C  D02PWF. On the first step this matter is handled differently as
C  explained above.
C
         IF ( .NOT. FIRST) THEN
            JFAIL = -1
            CALL D02PWF(TWANT,JFAIL)
            BADERR = JFAIL .GT. 0
            IF (BADERR) GO TO 100
         END IF
      END IF
C
C  Process output, decide whether to take another step.
C
   40 CONTINUE
C
      IF (INTP) THEN
C
C  Interpolation is possible with this METHD.  The integration has
C  already reached T. If this is past TWANT, GOBACK is set .TRUE. and
C  the answers are obtained by interpolation.
C
         GOBACK = DIR*(T-TWANT) .GE. ZERO
         IF (GOBACK) THEN
            JFAIL = -1
            CALL D02PXF(TWANT,'Both solution and derivative',NEQN,YGOT,
     *                  YPGOT,F,WORK,WORK(PRINTP),LNINTP,JFAIL)
            BADERR = JFAIL .GT. 0
            IF (BADERR) GO TO 100
            TGOT = TWANT
         END IF
      ELSE
C
C  Interpolation is not possible with this METHD, so output is obtained
C  by integrating to TWANT = TND.  Both YGOT(*) and YPGOT(*) are then
C  already loaded with the solution at TWANT by D02PDF.
C
         GOBACK = T .EQ. TWANT
         IF (GOBACK) TGOT = TWANT
      END IF
C
C  Updating of YMAX(*) is done here to account for the fact that when
C  interpolation is done, the integration goes past TGOT.  Note that
C  YGOT(*) is not defined until D02PDF is called.  YMAX(*) was
C  initialized at TSTRT from values stored in WORK(*), so only needs
C  to be updated for T different from TSTRT.
C
      IF (T.NE.TSTRT) THEN
         DO 60 L = 1, NEQN
            YMAX(L) = MAX(YMAX(L),ABS(YGOT(L)))
   60    CONTINUE
      END IF
C
C  If done, go to the exit point.
C
      IF (GOBACK) GO TO 100
C
C  Take a step with D02PDF in the direction of TND. On exit, the
C  solution is advanced to TNOW.  The way D02PDF is written, the
C  approximate solution at TNOW is available in both YGOT(*) and in
C  WORK(*).  If output is obtained by stepping to the end
C  (TNOW = TWANT = TND), YGOT(*) can be returned directly. If output
C  is obtained by interpolation, the subroutine D02PXF that does this
C  uses the values in WORK(*) for its computations and places the
C  approximate solution at TWANT in the array YGOT(*) for return to the
C  calling program. The approximate derivative is handled in the same
C  way. TNOW is output from D02PDF and is actually a copy of T declared
C  above in a common block.
C
      IF (IFAIL.LE.0) THEN
         JFAIL = -13
      ELSE
         JFAIL = IFAIL
      END IF
      CALL D02PDF(F,TNOW,YGOT,YPGOT,WORK,JFAIL)
      IER = JFAIL
C
C  A successful step by D02PDF is indicated by JFAIL = 0 or = 2.
      IF (JFAIL.EQ.0) THEN
         GO TO 40
      ELSE IF (JFAIL.EQ.2) THEN
C
C  Supplement the warning message written in D02PDF.
         NREC = 3
         WRITE (REC,FMT='(A,A/A,A/A)')
     *     ' ** The last message was produced on a call to D02PDF from '
     *     , 'D02PCF.',
     *   ' ** In D02PCF the appropriate action is to change to METHOD ='
     *     , ' 2,',
     *    ' ** or, if insufficient memory is available, to METHOD = 1. '
      ELSE IF (JFAIL.GT.1) THEN
         NREC = 1
         WRITE (REC,FMT='(A,A)')
     *     ' ** The last message was produced on a call to D02PDF from '
     *     , 'D02PCF.'
      ELSE
         BADERR = .TRUE.
      END IF
      TGOT = T
C
C  Update YMAX(*) before the return.
      DO 80 L = 1, NEQN
         YMAX(L) = MAX(YMAX(L),ABS(YGOT(L)))
   80 CONTINUE
C
C  Exit point for D02PCF.
C
  100 CONTINUE
C
      IF (BADERR) THEN
         IER = 1
         NREC = 4
         WRITE (REC,FMT='(A/A/A/A)')
     * ' ** An internal call by D02PCF to a subroutine resulted in an  '
     *     ,
     *     ' ** error that should not happen.  Check your program      '
     *     ,
     *     ' ** carefully for array sizes, correct number of arguments,'
     *     , ' ** type mismatches ... .'
      END IF
C
      IF (IER.NE.1) TLAST = TGOT
C
C  All exits are done here after a call to D02PDP to report
C  what happened and set IFAIL.
C
      CALL D02PDP(IER,SRNAME,NREC,IFAIL)
C
      RETURN
      END
