      SUBROUTINE D02PDF(F,TNOW,YNOW,YPNOW,WORK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C$$$$ SUBROUTINE CT $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C  If you are not familiar with the code D02PDF and how it is used in
C  conjunction with D02PVF to solve initial value problems, you should
C  study the document file RKSUITE.DOC carefully before attempting to
C  use the code. The following "Brief Reminder" is intended only to
C  remind you of the meaning, type, and size requirements of the
C  arguments.
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
C     Using the input values of the independent variable T and the
C     solution components Y(*), for each L = 1,2,...,NEQ evaluate the
C     differential equation for the derivative of the Lth solution
C     component and place the value in YP(L).  Do not alter the input
C     values of T and Y(*).
C  RETURN
C  END
C
C  OUTPUT VARIABLES
C
C     TNOW      - DOUBLE PRECISION
C                 Current value of the independent variable.
C     YNOW(*)   - DOUBLE PRECISION array of length NEQ
C                 Approximation to the true solution at TNOW.
C     YPNOW(*)  - DOUBLE PRECISION array of length NEQ
C                 Approximation to the first derivative of the
C                 true solution at TNOW.
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
C                       SUCCESS.  A STEP WAS TAKEN TO TNOW.
C                 = 0 - Complete success.
C
C                       "SOFT" FAILURES
C                 = 2 - Warning:  You have obtained an answer by
C                       integrating to TEND (TNOW = TEND).  You have
C                       done this at least 100 times, and monitoring of
C                       the computation reveals that this way of getting
C                       output has degraded the efficiency of the code.
C                       If you really need answers at so many specific
C                       points, it would be more efficient to get them
C                       with D02PXF.  (If METHOD = 3, you would need to
C                       change METHOD and restart from TNOW, YNOW(*) by
C                       a call to D02PVF.)  If you wish to continue as
C                       you are, you may.
C                 = 3 - Warning:  A considerable amount of work has been
C                       expended. To continue the integration, just call
C                       D02PDF again.
C                 = 4 - Warning:  It appears that this problem is
C                       "stiff". You really should change to another
C                       code that is intended for such problems, but if
C                       you insist, you can continue with D02PDF by
C                       calling it again.
C
C                       "HARD" FAILURES
C                 = 5 - You are asking for too much accuracy. You cannot
C                       continue integrating this problem.
C                 = 6 - The global error assessment may not be reliable
C                       beyond the current point in the integration.
C                       You cannot continue integrating this problem.
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
      PARAMETER         (SRNAME='D02PDF')
      LOGICAL           ASK, TELL
      PARAMETER         (ASK=.TRUE.,TELL=.FALSE.)
      INTEGER           MINUS1, MINUS2
      PARAMETER         (MINUS1=-1,MINUS2=-2)
      INTEGER           MAXFCN
      PARAMETER         (MAXFCN=5000)
      DOUBLE PRECISION  ZERO, PT1, PT9, ONE, TWO, HUNDRD
      PARAMETER         (ZERO=0.0D+0,PT1=0.1D+0,PT9=0.9D+0,ONE=1.0D+0,
     *                  TWO=2.0D+0,HUNDRD=100.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TNOW
      INTEGER           IFAIL
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(*), YNOW(*), YPNOW(*)
C     .. Subroutine Arguments ..
      EXTERNAL          F
C     .. Scalars in Common ..
      DOUBLE PRECISION  COST, CUBRMC, DIR, DWARF, EXPON, H, HOLD, HSTRT,
     *                  LOCMAX, MAXERR, MCHEPS, RNDOFF, RS, RS1, RS2,
     *                  RS3, RS4, SAFETY, SQRRMC, STBRAD, T, TANANG,
     *                  TINY, TND, TOLD, TOLR, TOOSML, TSTRT
      INTEGER           FLSTP, GNFCN, LNINTP, LSTSTG, MAXTRY, NEQN,
     *                  NFCN, NSEC, OKSTP, ORDER, OUTCH, PRERST, PRINTP,
     *                  PRSCR, PRSTGS, PRTHRS, PRWT, PRY, PRYOLD, PRYP,
     *                  PRZERR, PRZERS, PRZSTG, PRZY, PRZYNU, PRZYP,
     *                  SVNFCN
      LOGICAL           ERASFL, ERASON, FIRST, FSAL, LAST, UTASK
C     .. Arrays in Common ..
      CHARACTER*80      REC(10)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA, BETA, ERR, ERROLD, HAVG, HMIN, HTRY, TAU,
     *                  TEMP1, TEMP2, YPNORM
      INTEGER           IER, JFLSTP, L, NREC, NTEND, POINT, STATE, YNEW,
     *                  YPOLD
      LOGICAL           CHKEFF, FAILED, MAIN, PHASE1, PHASE2, PHASE3,
     *                  TOOMCH
C     .. External Subroutines ..
      EXTERNAL          D02PDM, D02PDP, D02PDV, D02PDW, D02PDZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SIGN
C     .. Common blocks ..
      COMMON            /AD02PD/TSTRT, TND, DIR, HSTRT, TOLR, NEQN
      COMMON            /BD02PD/T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP,
     *                  FLSTP, FIRST, LAST
      COMMON            /CD02PD/PRTHRS, PRERST, PRWT, PRYOLD, PRSCR,
     *                  PRY, PRYP, PRSTGS, PRINTP, LNINTP
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
      SAVE              /AD02PD/, /BD02PD/, /CD02PD/, /ED02PD/,
     *                  /FD02PD/, /GD02PD/, /HD02PD/, /JD02PD/, JFLSTP,
     *                  NTEND, ERROLD, HAVG, PHASE2, YNEW, YPOLD, CHKEFF
C     .. Executable Statements ..
C
      IER = 0
      NREC = 0
C
C  Is it permissible to call D02PDF?
C
      CALL D02PDM(ASK,'D02PVF',STATE)
      IF (STATE.EQ.1) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT='(A)')
     *   ' ** A catastrophic error has already been detected elsewhere.'
         GO TO 180
      END IF
      IF (STATE.EQ.MINUS1) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT='(A)')
     *     ' ** You have not called D02PVF, so you cannot use D02PDF.'
         GO TO 180
      END IF
      IF (UTASK) THEN
         CALL D02PDM(ASK,'D02PCF',STATE)
         IF (STATE.NE.MINUS2) THEN
            IER = 1
            NREC = 2
            WRITE (REC,FMT='(A/A)')
     *  ' ** You have called D02PDF after you specified in D02PVF that '
     *        ,
     *        ' ** you were going to use D02PCF. This is not permitted.'
            UTASK = .FALSE.
            GO TO 180
         END IF
      END IF
      CALL D02PDM(ASK,SRNAME,STATE)
      IF (STATE.EQ.5 .OR. STATE.EQ.6) THEN
         IER = 1
         NREC = 3
         WRITE (REC,FMT='(A/A/A)')
     *     ' ** D02PDF has already returned with IFAIL set to 5 or 6.',
     *     ' ** You cannot continue integrating this problem. You must '
     *     , ' ** call D02PVF to start another problem.'
         GO TO 180
      END IF
C
      IF (FIRST) THEN
C
C  First call in an integration -- initialize everything.
C
         CHKEFF = .FALSE.
         NTEND = 0
         JFLSTP = 0
C
C  A scratch area of WORK(*) starting at PRSCR is used to hold two
C  arrays in this subroutine: the higher order approximate solution at
C  the end of a step and the approximate derivative of the solution at
C  the end of the last step. To make this clear, local pointers YNEW and
C  YPOLD are used.
C
         YNEW = PRSCR
         YPOLD = PRSCR
C
C  For this first step T was initialized to TSTRT in D02PVF and the
C  starting values YSTART(*) were loaded into the area of WORK(*)
C  reserved for the current solution approximation starting at location
C  PRY. The derivative is now computed and stored in WORK(*) starting
C  at PRYP. Subsequently these arrays are copied to the output vectors
C  YNOW(*) and YPNOW(*).
C
         CALL F(T,WORK(PRY),WORK(PRYP))
         NFCN = NFCN + 1
         DO 20 L = 1, NEQN
            YNOW(L) = WORK(PRY-1+L)
            YPNOW(L) = WORK(PRYP-1+L)
   20    CONTINUE
C
C  Set dependent variables for error assessment.
C
         IF (ERASON) THEN
            DO 40 L = 1, NEQN
               WORK(PRZY-1+L) = YNOW(L)
               WORK(PRZYP-1+L) = YPNOW(L)
   40       CONTINUE
         END IF
C
C  The weights for the control of the error depend on the size of the
C  solution at the beginning and at the end of the step. On the first
C  step we do not have all this information. Whilst determining the
C  initial step size we initialize the weight vector to the larger of
C  abs(Y(L)) and the threshold for this component.
C
         DO 60 L = 1, NEQN
            WORK(PRWT-1+L) = MAX(ABS(YNOW(L)),WORK(PRTHRS-1+L))
   60    CONTINUE
C
C  If HSTRT is equal to zero, the code is to find an on-scale initial
C  step size H.  D02PDF has an elaborate scheme of three phases for
C  finding such an H, and some preparations are made earlier. In D02PVF
C  an upper bound is placed on H that reflects the scale of the
C  independent variable. When UTASK is .TRUE., D02PCF refines this bound
C  using the first output point.  Here in D02PDF PHASE1 applies a rule
C  of thumb based on the error control, the order of the the formula,
C  and the size of the initial slope to get a crude approximation to an
C  on-scale H.  PHASE2 may reduce H in the course of taking the first
C  step.  PHASE3 repeatedly adjusts H and retakes the first step until H
C  is on scale.
C
C  A guess for the magnitude of the first step size H can be provided to
C  D02PVF as HSTART. If it is too big or too small, it is ignored and
C  the automatic determination of an on-scale initial step size is
C  activated.  If it is acceptable, H is set to HSTART in D02PVF.  Even
C  when H is supplied to D02PDF, PHASE3 of the scheme for finding an
C  on-scale initial step size is made active so that the code can deal
C  with a bad guess.
C
         PHASE1 = HSTRT .EQ. ZERO
         PHASE2 = PHASE1
         PHASE3 = .TRUE.
         IF (PHASE1) THEN
            H = ABS(H)
            YPNORM = ZERO
            DO 80 L = 1, NEQN
               IF (ABS(YNOW(L)).NE.ZERO) THEN
                  YPNORM = MAX(YPNORM,ABS(YPNOW(L))/WORK(PRWT-1+L))
               END IF
   80       CONTINUE
            TAU = TOLR**EXPON
            IF (H*YPNORM.GT.TAU) H = TAU/YPNORM
            HMIN = MAX(TINY,TOOSML*MAX(ABS(TSTRT),ABS(TND)))
            H = DIR*MAX(H,HMIN)
            PHASE1 = .FALSE.
         END IF
C
      ELSE
C
C Continuation call
C
         IF (LAST) THEN
            IER = 1
            NREC = 3
            WRITE (REC,FMT='(A,D13.5,A/A/A)')
     *        ' ** You have already reached TEND ( = ', TND, ').',
     *        ' ** To integrate further with the same problem you must '
     *        , ' ** call the routine D02PWF with a new value of TEND.'
            GO TO 180
         END IF
      END IF
C
C  Begin computation of a step here.
C
      FAILED = .FALSE.
C
  100 CONTINUE
      H = SIGN(ABS(H),DIR)
C
C  Reduce the step size if necessary so that the code will not step
C  past TND.  "Look ahead" to prevent unnecessarily small step sizes.
C
      LAST = DIR*((T+H)-TND) .GE. ZERO
      IF (LAST) THEN
         H = TND - T
      ELSE IF (DIR*((T+TWO*H)-TND).GE.ZERO) THEN
         H = (TND-T)/TWO
      END IF
C
C  When the integrator is at T and attempts a step of H, the function
C  defining the differential equations will be evaluated at a number of
C  arguments between T and T+H.  If H is too small, these arguments
C  cannot be clearly distinguished in the precision available.
C
      HMIN = MAX(TINY,TOOSML*MAX(ABS(T),ABS(T+H)))
      IF (ABS(H).LT.HMIN) THEN
         IER = 5
         NREC = 3
         WRITE (REC,FMT='(A/A,D13.5,A,D13.5/A)')
     *   ' ** In order to satisfy your error requirements D02PDF would '
     *     , ' ** have to use a step size of ', H, ' at TNOW = ', T,
     *     ' ** This is too small for the machine precision.'
         GO TO 180
      END IF
C
C  Monitor the impact of output on the efficiency of the integration.
C
      IF (CHKEFF) THEN
         NTEND = NTEND + 1
         IF (NTEND.GE.100 .AND. NTEND.GE.OKSTP/3) THEN
            IER = 2
            NREC = 6
            WRITE (REC,FMT='(A/A/A/A/A/A)')
     *        ' ** More than 100 output points have been obtained by ',
     *    ' ** integrating to TEND.  They have been sufficiently close '
     *        ,
     *  ' ** to one another that the efficiency of the integration has '
     *        ,
     *  ' ** been degraded. It would probably be (much) more efficient '
     *        ,
     *       ' ** to obtain output by interpolating with D02PXF (after '
     *        , ' ** changing to METHOD=2 if you are using METHOD = 3).'
            NTEND = 0
            GO TO 180
         END IF
      END IF
C
C  Check for stiffness and for too much work.  Stiffness can be
C  checked only after a successful step.
C
      IF ( .NOT. FAILED) THEN
C
C  Check for too much work.
         TOOMCH = NFCN .GT. MAXFCN
         IF (TOOMCH) THEN
            IER = 3
            NREC = 3
            WRITE (REC,FMT='(A,I6,A/A/A)') ' ** Approximately ', MAXFCN,
     *        ' function evaluations have been ',
     *        ' ** used to compute the solution since the integration ',
     *        ' ** started or since this message was last printed.'
C
C  After this warning message, NFCN is reset to permit the integration
C  to continue.  The total number of function evaluations in the primary
C  integration is SVNFCN + NFCN.
C
            SVNFCN = SVNFCN + NFCN
            NFCN = 0
         END IF
C
C  Check for stiffness.  NREC is passed on to D02PDV because when
C  TOOMCH = .TRUE. and stiffness is diagnosed, the message about too
C  much work is augmented inside D02PDV to explain that it is due to
C  stiffness.
C
         CALL D02PDV(F,HAVG,JFLSTP,TOOMCH,MAXFCN,WORK,IER,NREC)
C
         IF (IER.NE.0) GO TO 180
      END IF
C
C  Take a step.  Whilst finding an on-scale H (PHASE2 = .TRUE.), the
C  input value of H might be reduced (repeatedly), but it will not be
C  reduced
C  below HMIN.  The local error is estimated, a weight vector is formed,
C  and a weighted maximum norm, ERR, of the local error is returned.
C  The variable MAIN is input as .TRUE. to tell D02PDZ that this is the
C  primary, or "main", integration.
C
C  H resides in the common block /BD02PD/ which is used by both D02PDF
C  and D02PDZ; since it may be changed inside D02PDZ, a local copy is
C  made to ensure portability of the code.
C
      MAIN = .TRUE.
      HTRY = H
      CALL D02PDZ(F,NEQN,T,WORK(PRY),WORK(PRYP),WORK(PRSTGS),TOLR,HTRY,
     *            WORK(PRWT),WORK(YNEW),WORK(PRERST),ERR,MAIN,HMIN,
     *            WORK(PRTHRS),PHASE2)
      H = HTRY
C
C  Compare the norm of the local error to the tolerance.
C
      IF (ERR.GT.TOLR) THEN
C
C  Failed step.  Reduce the step size and try again.
C
C  First step: Terminate PHASE3 of the search for an on-scale step size.
C              The step size is not on scale, so ERR may not be
C              accurate; reduce H by a fixed factor. Failed attempts to
C              take the first step are not counted.
C  Later step: Use ERR to compute an "optimal" reduction of H. More than
C              one failure indicates a difficulty with the problem and
C              an ERR that may not be accurate, so reduce H by a fixed
C              factor.
C
         IF (FIRST) THEN
            PHASE3 = .FALSE.
            ALPHA = RS1
         ELSE
            FLSTP = FLSTP + 1
            JFLSTP = JFLSTP + 1
            IF (FAILED) THEN
               ALPHA = RS1
            ELSE
               ALPHA = SAFETY*(TOLR/ERR)**EXPON
               ALPHA = MAX(ALPHA,RS1)
            END IF
         END IF
         H = ALPHA*H
         FAILED = .TRUE.
         GO TO 100
      END IF
C
C  Successful step.
C
C  Predict a step size appropriate for the next step.  After the first
C  step the prediction can be refined using an idea of H.A. Watts that
C  takes account of how well the prediction worked on the previous step.
C
      BETA = (ERR/TOLR)**EXPON
      IF ( .NOT. FIRST) THEN
         TEMP1 = (ERR**EXPON)/H
         TEMP2 = (ERROLD**EXPON)/HOLD
         IF (TEMP1.LT.TEMP2*HUNDRD .AND. TEMP2.LT.TEMP1*HUNDRD) THEN
            BETA = BETA*(TEMP1/TEMP2)
         END IF
      END IF
      ALPHA = RS3
      IF (SAFETY.LT.BETA*ALPHA) ALPHA = SAFETY/BETA
C
C  On the first step a search is made for an on-scale step size. PHASE2
C  of the scheme comes to an end here because a step size has been found
C  that is both successful and has a credible local error estimate.
C  Except in the special case that the first step is also the last, the
C  step is repeated in PHASE3 as long as an increase greater than RS2
C  appears possible.  An increase as big as RS3 is permitted. A step
C  failure terminates PHASE3.
C
      IF (FIRST) THEN
         PHASE2 = .FALSE.
         PHASE3 = PHASE3 .AND. .NOT. LAST .AND. (ALPHA.GT.RS2)
         IF (PHASE3) THEN
            H = ALPHA*H
            GO TO 100
         END IF
      END IF
C
C  After getting on scale, step size changes are more restricted.
C
      ALPHA = MIN(ALPHA,RS)
      IF (FAILED) ALPHA = MIN(ALPHA,ONE)
      ALPHA = MAX(ALPHA,RS1)
      HOLD = H
      H = ALPHA*H
C
C  For the diagnosis of stiffness, an average accepted step size, HAVG,
C  must be computed and SAVEd.
C
      IF (FIRST) THEN
         HAVG = HOLD
      ELSE
         HAVG = PT9*HAVG + PT1*HOLD
      END IF
C
      FIRST = .FALSE.
      ERROLD = ERR
      TOLD = T
C
C  Take care that T is set to precisely TND when the end of the
C  integration is reached.
C
      IF (LAST) THEN
         T = TND
      ELSE
         T = T + HOLD
      END IF
C
C  Increment counter on accepted steps.  Note that successful steps
C  that are repeated whilst getting on scale are not counted.
C
      OKSTP = OKSTP + 1
C
C  Advance the current solution and its derivative.  (Stored in WORK(*)
C  with the first location being PRY and PRYP, respectively.)  Update
C  the previous solution and its derivative.  (Stored in WORK(*) with
C  the first location being PRYOLD and YPOLD, respectively.)  Note that
C  the previous derivative will overwrite YNEW(*).
C
      DO 120 L = 1, NEQN
         WORK(PRYOLD-1+L) = WORK(PRY-1+L)
         WORK(PRY-1+L) = WORK(YNEW-1+L)
         WORK(YPOLD-1+L) = WORK(PRYP-1+L)
  120 CONTINUE
C
      IF (FSAL) THEN
C
C  When FSAL = .TRUE., YP(*) is the last stage of the step.
C
         POINT = PRSTGS + (LSTSTG-1)*NEQN
         DO 140 L = 1, NEQN
            WORK(PRYP-1+L) = WORK(POINT-1+L)
  140    CONTINUE
      ELSE
C
C  Call F to evaluate YP(*).
C
         CALL F(T,WORK(PRY),WORK(PRYP))
         NFCN = NFCN + 1
      END IF
C
C  If global error assessment is desired, advance the secondary
C  integration from TOLD to T.
C
      IF (ERASON) THEN
         CALL D02PDW(F,NEQN,WORK(PRY),TOLR,WORK(PRWT),WORK(PRZY),
     *               WORK(PRZYP),WORK(PRZERR),WORK(PRZYNU),WORK(PRZERS),
     *               WORK(PRZSTG),IER)
         IF (IER.EQ.6) THEN
C
C  The global error estimating procedure has broken down. Treat it as a
C  failed step. The solution and derivative are reset to their values at
C  the beginning of the step since the last valid error assessment
C  refers to them.
C
            OKSTP = OKSTP - 1
            ERASFL = .TRUE.
            LAST = .FALSE.
            T = TOLD
            H = HOLD
            DO 160 L = 1, NEQN
               WORK(PRY-1+L) = WORK(PRYOLD-1+L)
               WORK(PRYP-1+L) = WORK(YPOLD-1+L)
  160       CONTINUE
            IF (OKSTP.GT.0) THEN
               NREC = 2
               WRITE (REC,FMT='(A/A,D13.5,A)')
     * ' ** The global error assessment may not be reliable for T past '
     *           , ' ** TNOW = ', T,
     *           '.  The integration is being terminated.'
            ELSE
               NREC = 2
               WRITE (REC,FMT='(A/A)')
     *   ' ** The global error assessment algorithm failed at the start'
     *           ,
     *      ' ** the integration.  The integration is being terminated.'
            END IF
            GO TO 180
         END IF
      END IF
C
C
C  Exit point for D02PDF
C
  180 CONTINUE
C
C  Set the output variables and flag that interpolation is permitted
C
      IF (IER.NE.1) THEN
         TNOW = T
         LAST = TNOW .EQ. TND
         CHKEFF = LAST
         DO 200 L = 1, NEQN
            YNOW(L) = WORK(PRY-1+L)
            YPNOW(L) = WORK(PRYP-1+L)
  200    CONTINUE
         IF (IER.EQ.0) THEN
            STATE = MINUS2
            CALL D02PDM(TELL,'D02PXF',STATE)
         END IF
      END IF
C
C  Call D02PDP to report what happened and set IFAIL
C
      CALL D02PDP(IER,SRNAME,NREC,IFAIL)
C
      RETURN
      END
