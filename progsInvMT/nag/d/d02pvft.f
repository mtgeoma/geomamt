      SUBROUTINE D02PVF(NEQ,TSTART,YSTART,TEND,TOL,THRES,METHOD,TASK,
     *                  ERRASS,HSTART,WORK,LENWRK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C$$$$ SUBROUTINE SETUP $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C  If you are not familiar with the code D02PVF and how it is used in
C  conjunction with D02PCF or D02PDF to solve initial value problems,
C  you should study the document file RKSUITE.DOC carefully before
C  attempting to use the code.
C  The following "Brief Reminder" is intended only to remind you of the
C  meaning, type, and size requirements of the arguments.
C
C  The environmental parameters OUTCH, MCHEPS, and DWARF are used in the
C  following description.  To find out their values
C
C       CALL D02PVX(OUTCH,MCHEPS,DWARF)
C
C  INPUT VARIABLES
C
C     NEQ       - INTEGER
C                 The number of differential equations in the system.
C                 Constraint: NEQ >= 1
C     TSTART    - DOUBLE PRECISION
C                 The initial value of the independent variable.
C     YSTART(*) - DOUBLE PRECISION array of length NEQ
C                 The vector of initial values of the solution
C                 components.
C     TEND      - DOUBLE PRECISION
C                 The integration proceeds from TSTART in the direction
C                 of TEND. You cannot go past TEND.
C                 Constraint: TEND must be clearly distinguishable from
C                 TSTART in the precision available.
C     TOL       - DOUBLE PRECISION
C                 The relative error tolerance.
C                 Constraint: 0.01D0 >= TOL >= 10*MCHEPS
C     THRES(*)  - DOUBLE PRECISION array of length NEQ
C                 THRES(L) is the threshold for the Ith solution
C                 component.
C                 Constraint: THRES(L) >= SQRT(DWARF)
C     METHOD    - INTEGER
C                 Specifies which Runge-Kutta pair is to be used.
C                  = 1 - use the (2,3) pair
C                  = 2 - use the (4,5) pair
C                  = 3 - use the (7,8) pair
C     TASK      - CHARACTER*(*)
C                 Only the first character of TASK is significant.
C                 TASK(1:1) = `U' or `u' - D02PCF is to be used
C                           = `C' or `c' - D02PDF is to be used
C                 Constraint: TASK(1:1) = `U'or `u' or`C' or `c'
C     ERRASS    - LOGICAL
C                 = .FALSE. - do not attempt to assess the true error.
C                 = .TRUE.  - assess the true error. Costs roughly twice
C                             as much as the integration with METHODs 2
C                             and 3, and three times with METHOD = 1.
C     HSTART    - DOUBLE PRECISION
C                 0.0D0     - select automatically the first step size.
C                 non-zero  - try HSTART for the first step.
C
C  WORKSPACE
C
C     WORK(*) - DOUBLE PRECISION array of length LENWRK
C               Do not alter the contents of this array after calling
C               D02PVF.
C
C  INPUT VARIABLES
C
C     LENWRK  - INTEGER
C               Length of WORK(*): How big LENWRK must be depends
C               on the task and how it is to be solved.
C
C               LENWRK = 32*NEQ is sufficient for all cases.
C
C               If storage is a problem, the least storage possible
C               in the various cases is:
C
C                 If TASK = `U' or `u', then
C                   if ERRASS = .FALSE. and
C                     METHOD = 1, LENWRK must be at least 10*NEQ
C                            = 2                          20*NEQ
C                            = 3                          16*NEQ
C                   if ERRASS = .TRUE. and
C                     METHOD = 1, LENWRK must be at least 15*NEQ
C                            = 2                          32*NEQ
C                            = 3                          21*NEQ
C
C                 If TASK = `C' or `c', then
C                   if ERRASS = .FALSE. and
C                     METHOD = 1, LENWRK must be at least 10*NEQ
C                            = 2                          14*NEQ
C                            = 3                          16*NEQ
C                   if ERRASS = .TRUE. and
C                     METHOD = 1, LENWRK must be at least 15*NEQ
C                            = 2                          26*NEQ
C                            = 3                          21*NEQ
C
C                 Warning:  To exploit the interpolation capability
C                 of METHODs 1 and 2, you have to call D02PXF.  This
C                 subroutine requires working storage in addition to
C                 that specified here.
C
C  In the event of a "catastrophic" failure to call D02PVF correctly,
C  the nature of the catastrophe is reported on the standard output
C  channel, regardless of the value of MESAGE.  Unless special provision
C  was made in advance (see RKSUITE.DOC), the computation then comes to
C  STOP.
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02PVF')
      INTEGER           MINUS1
      LOGICAL           TELL
      PARAMETER         (MINUS1=-1,TELL=.FALSE.)
      DOUBLE PRECISION  ONE, ZERO, PT01
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0,PT01=0.01D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HSTART, TEND, TOL, TSTART
      INTEGER           IFAIL, LENWRK, METHOD, NEQ
      LOGICAL           ERRASS
      CHARACTER*(*)     TASK
C     .. Array Arguments ..
      DOUBLE PRECISION  THRES(NEQ), WORK(LENWRK), YSTART(NEQ)
C     .. Scalars in Common ..
      DOUBLE PRECISION  COST, CUBRMC, DIR, DWARF, EXPON, H, HOLD, HSTRT,
     *                  LOCMAX, MAXERR, MCHEPS, RNDOFF, RS, RS1, RS2,
     *                  RS3, RS4, SAFETY, SQRRMC, STBRAD, T, TANANG,
     *                  TINY, TND, TOLD, TOLR, TOOSML, TSTRT
      INTEGER           FLSTP, GNFCN, LNINTP, LSTSTG, MAXTRY, METHD,
     *                  MINTP, NEQN, NFCN, NSEC, NSTAGE, OKSTP, ORDER,
     *                  OUTCH, PRERST, PRINTP, PRSCR, PRSTGS, PRTHRS,
     *                  PRWT, PRY, PRYOLD, PRYP, PRZERR, PRZERS, PRZSTG,
     *                  PRZY, PRZYNU, PRZYP, SVNFCN
      LOGICAL           ERASFL, ERASON, FIRST, FSAL, INTP, LAST, UTASK
C     .. Arrays in Common ..
      DOUBLE PRECISION  A(13,13), B(13), BHAT(13), C(13), E(7), R(11,6)
      INTEGER           PTR(13)
      CHARACTER*80      REC(10)
C     .. Local Scalars ..
      DOUBLE PRECISION  HMIN
      INTEGER           FREEPR, IER, L, LINTPL, LREQ, NREC, VECSTG
      LOGICAL           LEGALT, REQSTG
      CHARACTER         TASK1
C     .. External Subroutines ..
      EXTERNAL          D02PDM, D02PDP, D02PVY, D02PVZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SIGN
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
      COMMON            /FD02PD/MAXERR, LOCMAX, GNFCN, PRZSTG, PRZY,
     *                  PRZYP, PRZERS, PRZERR, PRZYNU, ERASON, ERASFL
      COMMON            /GD02PD/MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC,
     *                  TINY, OUTCH
      COMMON            /HD02PD/UTASK
      COMMON            /JD02PD/REC
C     .. Save statement ..
      SAVE              /AD02PD/, /BD02PD/, /CD02PD/, /DD02PD/,
     *                  /ED02PD/, /FD02PD/, /GD02PD/, /HD02PD/, /JD02PD/
C     .. Executable Statements ..
C
C  Clear previous flag values of subprograms in the suite.
C
      IER = MINUS1
      CALL D02PDM(TELL,SRNAME,IER)
C
      IER = 0
      NREC = 0
C
C  Fetch output channel and machine constants; initialise common
C  block /GD02PD/
C
      CALL D02PVY(METHOD)
C
C  Check for valid input of trivial arguments
      TASK1 = TASK(1:1)
      LEGALT = TASK1 .EQ. 'U' .OR. TASK1 .EQ. 'u' .OR. TASK1 .EQ.
     *         'C' .OR. TASK1 .EQ. 'c'
      IF ( .NOT. LEGALT) THEN
         IER = 1
         NREC = 3
         WRITE (REC,FMT='(A/A,A,A/A)')
     *     ' ** You have set the first character of ',
     *     ' ** TASK to be ''', TASK1, '''. It must be one of ',
     *     ' ** ''U'',''u'',''C'' or ''c''.'
      ELSE IF (NEQ.LT.1) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT='(A,I6,A)') ' ** You have set NEQ = ', NEQ,
     *     ' which is less than 1.'
      ELSE IF (METHOD.LT.1 .OR. METHOD.GT.3) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT='(A,I6,A)') ' ** You have set METHOD = ',
     *     METHOD, ' which is not 1, 2, or 3.'
      ELSE IF (TSTART.EQ.TEND) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT='(A,D13.5,A)')
     *     ' ** You have set TSTART = TEND = ', TSTART, '.'
      ELSE IF ((TOL.GT.PT01) .OR. (TOL.LT.RNDOFF)) THEN
         IER = 1
         NREC = 2
         WRITE (REC,FMT='(A,D13.5,A/A,D13.5,A)')
     *     ' ** You have set TOL = ', TOL,
     *     ' which is not permitted. The',
     *     ' ** range of permitted values is (', RNDOFF, ',0.01D0).'
      ELSE
         L = 1
   20    CONTINUE
         IF (THRES(L).LT.TINY) THEN
            IER = 1
            NREC = 2
            WRITE (REC,FMT='(A,I6,A,D13.5,A/A,D13.5,A)')
     *        ' ** You have set THRES(', L, ') to be ', THRES(L),
     *        ' which is ', ' ** less than the permitted minimum,',
     *        TINY, '.'
         END IF
         L = L + 1
         IF (IER.EQ.0 .AND. L.LE.NEQ) GO TO 20
      END IF
C
C  Return if error detected
C
      IF (IER.NE.0) GO TO 80
C
C  Set formula definitions and characteristics by means of arguments
C  in the call list and COMMON blocks /DD02PD/ and /ED02PD/
C
      CALL D02PVZ(METHOD,VECSTG,REQSTG,LINTPL)
C
C  Set option in /HD02PD/
      UTASK = TASK1 .EQ. 'U' .OR. TASK1 .EQ. 'u'
C
C  Initialise problem status in /AD02PD/ and /BD02PD/
      NEQN = NEQ
      TSTRT = TSTART
      TND = TEND
      T = TSTART
      TOLD = TSTART
      DIR = SIGN(ONE,TEND-TSTART)
C
C  In D02PDF the first step taken will have magnitude H. If
C  HSTRT = ABS(HSTART) is not equal to zero, H = HSTRT. If HSTRT is
C  equal to zero, the code is to find an on-scale initial step size
C  H. To start this process, H is set here to an upper bound on the
C  first step size that reflects the scale of the independent variable.
C  D02PCF has some additional information, namely the first output
C  point, that is used to refine this bound in D02PCF when UTASK is
C  .TRUE..  If HSTRT is not zero, but it is either too big or too small,
C  the input HSTART is ignored and HSTRT is set to zero to activate the
C  automatic determination of an on-scale initial step size.
C
      HSTRT = ABS(HSTART)
      HMIN = MAX(TINY,TOOSML*MAX(ABS(TSTART),ABS(TEND)))
      IF (HSTRT.GT.ABS(TEND-TSTART) .OR. HSTRT.LT.HMIN) HSTRT = ZERO
      IF (HSTRT.EQ.ZERO) THEN
         H = MAX(ABS(TEND-TSTART)/RS3,HMIN)
      ELSE
         H = HSTRT
      END IF
      HOLD = ZERO
      TOLR = TOL
      NFCN = 0
      SVNFCN = 0
      OKSTP = 0
      FLSTP = 0
      FIRST = .TRUE.
      LAST = .FALSE.
C
C  WORK(*) is partioned into a number of arrays using pointers. These
C  pointers are set in /CD02PD/.
      PRTHRS = 1
C                           the threshold values
      PRERST = PRTHRS + NEQ
C                           the error estimates
      PRWT = PRERST + NEQ
C                           the weights used in the local error test
      PRYOLD = PRWT + NEQ
C                           the previous value of the solution
      PRSCR = PRYOLD + NEQ
C                           scratch array used for the higher order
C                           approximate solution and for the previous
C                           value of the derivative of the solution
      PRY = PRSCR + NEQ
C                           the dependent variables
      PRYP = PRY + NEQ
C                           the derivatives
      PRSTGS = PRYP + NEQ
C                           intermediate stages held in an internal
C                           array STAGES(NEQ,VECSTG)
C
      FREEPR = PRSTGS + VECSTG*NEQ
C
C  Allocate storage for interpolation if the TASK = `U' or `u' was
C  specified. INTP and LINTPL returned by D02PVZ indicate whether there
C  is an interpolation scheme associated with the pair and how much
C  storage is required.
C
      PRINTP = 1
      LNINTP = 1
      IF (UTASK) THEN
         IF (INTP) THEN
            LNINTP = LINTPL*NEQ
            IF (REQSTG) THEN
               PRINTP = FREEPR
               FREEPR = PRINTP + LNINTP
            ELSE
               PRINTP = PRSTGS
               FREEPR = MAX(PRINTP+VECSTG*NEQ,PRINTP+LNINTP)
            END IF
         END IF
      END IF
C
C  Initialise state and allocate storage for global error assessment
C  using /FD02PD/
      GNFCN = 0
      MAXERR = ZERO
      LOCMAX = TSTART
      ERASON = ERRASS
      ERASFL = .FALSE.
      IF (ERRASS) THEN
C
C  Storage is required for the stages of a secondary integration. The
C  stages of the primary intergration can only be overwritten in the
C  cases where there is no interpolant or the interpolant does not
C  require information about the stages (e.g. METHOD 3 and METHOD 1,
C  respectively).
         IF ( .NOT. REQSTG) THEN
            PRZSTG = PRSTGS
         ELSE
            PRZSTG = FREEPR
            FREEPR = PRZSTG + VECSTG*NEQ
         END IF
         PRZY = FREEPR
         PRZYP = PRZY + NEQ
         PRZERS = PRZYP + NEQ
         PRZERR = PRZERS + NEQ
         PRZYNU = PRZERR + NEQ
         FREEPR = PRZYNU + NEQ
      ELSE
         PRZSTG = 1
         PRZY = 1
         PRZYP = 1
         PRZERS = 1
         PRZERR = 1
         PRZYNU = 1
      END IF
C
      LREQ = FREEPR - 1
C
C  Check for enough workspace and suitable range of integration
C
      IF (LENWRK.LT.LREQ) THEN
         IER = 1
         NREC = 2
         WRITE (REC,FMT='(A/A,I6,A,I6,A)')
     *    ' ** You have not supplied enough workspace. You gave LENWRK '
     *     , ' ** as', LENWRK, ', but it must be at least ', LREQ, '.'
      ELSE
         HMIN = MAX(TINY,TOOSML*MAX(ABS(TSTART),ABS(TEND)))
         IF (ABS(TEND-TSTART).LT.HMIN) THEN
            IER = 1
            NREC = 4
            WRITE (REC,FMT='(A/A/A,D13.5/A,D13.5,A)')
     *       ' ** You have set values for TEND and TSTART that are not '
     *        ,
     *   ' ** clearly distinguishable for the method and the precision '
     *        , ' ** of the computer being used. ABS(TEND-TSTART) is ',
     *        ABS(TEND-TSTART), ' ** but should be at least ', HMIN, '.'
         END IF
      END IF
C
C  Return if error detected
C
      IF (IER.NE.0) GO TO 80
C
C  Initialize elements of the workspace
      DO 40 L = 1, NEQ
         WORK(PRTHRS-1+L) = THRES(L)
         WORK(PRY-1+L) = YSTART(L)
   40 CONTINUE
C
C  Initialize the global error to zero when ERRASS = .TRUE.
      IF (ERRASS) THEN
         DO 60 L = 1, NEQ
            WORK(PRZERR-1+L) = ZERO
   60    CONTINUE
      END IF
C
   80 CONTINUE
C
      CALL D02PDP(IER,SRNAME,NREC,IFAIL)
C
      RETURN
      END
