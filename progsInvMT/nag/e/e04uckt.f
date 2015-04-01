      SUBROUTINE E04UCK(FIRST,DEBUG,DONE,IMPRVD,INFORM,MAXF,NUMF,NOUT,
     *                  ALFMAX,EPSAF,G0,TARGTG,FTRY,GTRY,TOLABS,TOLREL,
     *                  TOLTNY,ALFA,ALFBST,FBEST,GBEST)
C     MARK 13 RE-ISSUE.  NAG COPYRIGHT 1988.
C     MARK 16 REVISED. IER-1080 (JUL 1993).
C
C     ==================================================================
C     E04UCK  finds a sequence of improving estimates of a minimizer of
C     the univariate function f(alpha) in the interval (0,ALFMAX].
C     f(alpha) is a smooth function such that  f(0) = 0  and  f'(0) < 0.
C     E04UCK  requires both  f(alpha)  and  f'(alpha) to be evaluated at
C     points in the interval.  Estimates of the minimizer are computed
C     using safeguarded cubic interpolation.
C
C     Reverse communication is used to allow the calling program to
C     evaluate f and f'.  Some of the parameters must be set or tested
C     by the calling program.  The remainder would ordinarily be local
C     variables.
C
C     Input parameters (relevant to the calling program)
C     --------------------------------------------------
C
C     FIRST         must be .TRUE. on the first entry.
C                   It is subsequently altered by E04UCK.
C
C     DEBUG         specifies whether detailed output is wanted.
C
C     MAXF          is an upper limit on the number of times E04UCK is
C                   to be entered consecutively with DONE = .FALSE.
C                   (following an initial entry with FIRST = .TRUE.).
C
C     ALFA          is the first estimate of a minimizer.  ALFA is
C                   subsequently altered by E04UCK (see below).
C
C     ALFMAX        is the upper limit of the interval to be searched.
C
C     EPSAF         is an estimate of the absolute precision in the
C                   computed value of f(0).
C
C     FTRY, GTRY    are the values of f, f'  at the new point
C                   ALFA = ALFBST + XTRY.
C
C     G0            is the value of f'(0).  G0 must be negative.
C
C     TOLABS,TOLREL define a function TOL(ALFA) = TOLREL*ALFA + TOLABS
C                   such that if f has already been evaluated at ALFA,
C                   it will not be evaluated closer than TOL(ALFA).
C                   These values may be reduced by E04UCK.
C
C     TARGTG        is the target value of abs(f'(ALFA)). The search
C                   is terminated when
C                    abs(f'(ALFA)) le TARGTG and f(ALFA) lt 0.
C
C     TOLTNY        is the smallest value that TOLABS is allowed to be
C                   reduced to.
C
C     Output parameters (relevant to the calling program)
C     ---------------------------------------------------
C
C     IMPRVD        is .TRUE. if the previous ALFA was the best point so
C                   far.  Any related quantities should be saved by the
C                   calling program (e.g., gradient arrays) before
C                   paying attention to the variable DONE.
C
C     DONE = .FALSE.  means the calling program should evaluate
C                      FTRY = f(ALFA),  GTRY = f'(ALFA)
C                   for the new trial ALFA, and re-enter E04UCK.
C
C     DONE = .TRUE.   means that no new ALFA was calculated.  The value
C                   of INFORM gives the result of the search as follows
C
C                   INFORM = 1 means the search has terminated
C                              successfully with ALFBST < ALFMAX.
C
C                   INFORM = 2 means the search has terminated
C                              successfully with ALFBST = ALFMAX.
C
C                   INFORM = 3 means that the search failed to find a
C                              point of sufficient decrease in MAXF
C                              functions, but a lower point was found.
C
C                   INFORM = 4 means ALFMAX is so small that a search
C                              should not have been attempted.
C
C                   INFORM = 5 is never set by E04UCK.
C
C                   INFORM = 6 means the search has failed to find a
C                              useful step.  The interval of uncertainty
C                              is [0,B] with B < 2*TOLABS. A minimizer
C                              lies very close to ALFA = 0, or f'(0) is
C                              not sufficiently accurate.
C
C                   INFORM = 7 if no better point could be found after
C                              MAXF  function calls.
C
C                   INFORM = 8 means the input parameters were bad.
C                              ALFMAX le TOLTNY  or G0 ge zero.
C                              No function evaluations were made.
C
C     NUMF          counts the number of times E04UCK has been entered
C                   consecutively with DONE = .FALSE. (i.e., with a new
C                   function value FTRY).
C
C     ALFA          is the point at which the next function FTRY and
C                   derivative GTRY must be computed.
C
C     ALFBST        should be accepted by the calling program as the
C                   approximate minimizer, whenever E04UCK returns
C                   INFORM = 1 or 2 (and possibly 3).
C
C     FBEST, GBEST  will be the corresponding values of f, f'.
C
C
C     The following parameters retain information between entries
C     -----------------------------------------------------------
C
C     BRAKTD        is .FALSE. if f and f' have not been evaluated at
C                   the far end of the interval of uncertainty.  In this
C                   case, the point B will be at ALFMAX + TOL(ALFMAX).
C
C     CRAMPD        is .TRUE. if ALFMAX is very small (le TOLABS).  If
C                   the search fails, this indicates that a zero step
C                   should be taken.
C
C     EXTRAP        is .TRUE. if XW lies outside the interval of
C                   uncertainty.  In this case, extra safeguards are
C                   applied to allow for instability in the polynomial
C                   fit.
C
C     MOVED         is .TRUE. if a better point has been found, i.e.,
C                   ALFBST gt 0.
C
C     WSET          records whether a second-best point has been
C                   determined it will always be .TRUE. when convergence
C                   is tested.
C
C     NSAMEA        is the number of consecutive times that the
C                   left-hand end point of the interval of uncertainty
C                   has remained the same.
C
C     NSAMEB        similarly for the right-hand end.
C
C     A, B, ALFBST  define the current interval of uncertainty.
C                   A minimizer lies somewhere in the interval
C                   [ALFBST + A, ALFBST + B].
C
C     ALFBST        is the best point so far.  It is always at one end
C                   of the interval of uncertainty.  hence we have
C                   either  A lt 0,  B = 0  or  A = 0,  B gt 0.
C
C     FBEST, GBEST  are the values of f, f' at the point ALFBST.
C
C     FACTOR        controls the rate at which extrapolated estimates
C                   of ALFA may expand into the interval of uncertainty.
C                   FACTOR is not used if a minimizer has been bracketed
C                   (i.e., when the variable BRAKTD is .TRUE.).
C
C     FW, GW        are the values of f, f' at the point ALFBST + XW.
C                   they are not defined until WSET is .TRUE..
C
C     XTRY          is the trial point in the shifted interval (A, B).
C
C     XW            is such that  ALFBST + XW  is the second-best point.
C                   it is not defined until  WSET  is .TRUE..
C                   in some cases,  XW  will replace a previous  XW
C                   that has a lower function but has just been excluded
C                   from the interval of uncertainty.
C
C
C     Systems Optimization Laboratory, Stanford University, California.
C     Original version February 1982.  Rev. May 1983.
C     Original f77 version 22-August-1985.
C     This version of E04UCK dated  14-Sep-92.
C     ==================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, POINT1, HALF
      PARAMETER         (ZERO=0.0D+0,POINT1=0.1D+0,HALF=0.5D+0)
      DOUBLE PRECISION  ONE, THREE, FIVE
      PARAMETER         (ONE=1.0D+0,THREE=3.0D+0,FIVE=5.0D+0)
      DOUBLE PRECISION  TEN, ELEVEN
      PARAMETER         (TEN=1.0D+1,ELEVEN=1.1D+1)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, ALFBST, ALFMAX, EPSAF, FBEST, FTRY, G0,
     *                  GBEST, GTRY, TARGTG, TOLABS, TOLREL, TOLTNY
      INTEGER           INFORM, MAXF, NOUT, NUMF
      LOGICAL           DEBUG, DONE, FIRST, IMPRVD
C     .. Local Scalars ..
      DOUBLE PRECISION  A, ABSR, ARTIFA, ARTIFB, B, DAUX, DTRY, FACTOR,
     *                  FW, GW, Q, R, S, SCALE, TOL, TOLMAX, TRUEA,
     *                  TRUEB, XMIDPT, XTRY, XW
      INTEGER           NSAMEA, NSAMEB
      LOGICAL           BADFUN, BRAKTD, CLOSEF, CRAMPD, EXTRAP, FITOK,
     *                  FOUND, MOVED, QUITF, QUITI, SETXW, WSET
C     .. Local Arrays ..
      CHARACTER*120     REC(7)
C     .. External Subroutines ..
      EXTERNAL          X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Save statement ..
      SAVE              BRAKTD, CRAMPD, EXTRAP, MOVED, WSET, NSAMEA,
     *                  NSAMEB, A, B, FACTOR, XTRY, XW, FW, GW, TOLMAX
C     .. Executable Statements ..
C
C     ------------------------------------------------------------------
C     Local variables
C     ===============
C
C     CLOSEF     is .TRUE. if the new function FTRY is within EPSAF of
C                FBEST (up or down).
C
C     FOUND      is .TRUE. if the sufficient decrease conditions hold at
C                ALFBST.
C
C     QUITF      is .TRUE. when  MAXF  function calls have been made.
C
C     QUITI      is .TRUE. when the interval of uncertainty is less than
C                2*TOL.
C     ------------------------------------------------------------------
C
      BADFUN = .FALSE.
      QUITF = .FALSE.
      QUITI = .FALSE.
      IMPRVD = .FALSE.
C
      IF (FIRST) THEN
C        ---------------------------------------------------------------
C        First entry.  Initialize various quantities, check input data
C        and prepare to evaluate the function at the initial ALFA.
C        ---------------------------------------------------------------
         FIRST = .FALSE.
         NUMF = 0
         ALFBST = ZERO
         BADFUN = ALFMAX .LE. TOLTNY .OR. G0 .GE. ZERO
         DONE = BADFUN
         MOVED = .FALSE.
C
         IF ( .NOT. DONE) THEN
            BRAKTD = .FALSE.
            CRAMPD = ALFMAX .LE. TOLABS
            EXTRAP = .FALSE.
            WSET = .FALSE.
            NSAMEA = 0
            NSAMEB = 0
C
            TOLMAX = TOLABS + TOLREL*ALFMAX
            A = ZERO
            B = ALFMAX + TOLMAX
            FACTOR = FIVE
            TOL = TOLABS
            XTRY = ALFA
            IF (DEBUG) THEN
               WRITE (REC,FMT=99999) G0, TOLABS, ALFMAX, TARGTG, TOLREL,
     *           EPSAF, CRAMPD
               CALL X04BAY(NOUT,4,REC)
            END IF
         END IF
      ELSE
C        ---------------------------------------------------------------
C        Subsequent entries. The function has just been evaluated at
C        ALFA = ALFBST + XTRY,  giving FTRY and GTRY.
C        ---------------------------------------------------------------
         IF (DEBUG) THEN
            WRITE (REC,FMT=99998) ALFA, FTRY, GTRY
            CALL X04BAY(NOUT,2,REC)
         END IF
C
         NUMF = NUMF + 1
         NSAMEA = NSAMEA + 1
         NSAMEB = NSAMEB + 1
C
         IF ( .NOT. BRAKTD) THEN
            TOLMAX = TOLABS + TOLREL*ALFMAX
            B = ALFMAX - ALFBST + TOLMAX
         END IF
C
C        See if the new step is better.  If ALFA is large enough that
C        FTRY can be distinguished numerically from zero,  the function
C        is required to be sufficiently negative.
C
         CLOSEF = ABS(FTRY-FBEST) .LE. EPSAF
         IF (CLOSEF) THEN
            IMPRVD = ABS(GTRY) .LE. ABS(GBEST)
         ELSE
            IMPRVD = FTRY .LT. FBEST
         END IF
C
         IF (IMPRVD) THEN
C
C           We seem to have an improvement.  The new point becomes the
C           origin and other points are shifted accordingly.
C
            FW = FBEST
            FBEST = FTRY
            GW = GBEST
            GBEST = GTRY
            ALFBST = ALFA
            MOVED = .TRUE.
C
            A = A - XTRY
            B = B - XTRY
            XW = ZERO - XTRY
            WSET = .TRUE.
            EXTRAP = XW .LT. ZERO .AND. GBEST .LT. ZERO .OR. XW .GT.
     *               ZERO .AND. GBEST .GT. ZERO
C
C           Decrease the length of the interval of uncertainty.
C
            IF (GTRY.LE.ZERO) THEN
               A = ZERO
               NSAMEA = 0
            ELSE
               B = ZERO
               NSAMEB = 0
               BRAKTD = .TRUE.
            END IF
         ELSE
C
C           The new function value is not better than the best point so
C           far.  The origin remains unchanged but the new point may
C           qualify as XW.  XTRY must be a new bound on the best point.
C
            IF (XTRY.LE.ZERO) THEN
               A = XTRY
               NSAMEA = 0
            ELSE
               B = XTRY
               NSAMEB = 0
               BRAKTD = .TRUE.
            END IF
C
C           If XW has not been set or FTRY is better than FW, update the
C           points accordingly.
C
            IF (WSET) THEN
               SETXW = FTRY .LT. FW .OR. .NOT. EXTRAP
            ELSE
               SETXW = .TRUE.
            END IF
C
            IF (SETXW) THEN
               XW = XTRY
               FW = FTRY
               GW = GTRY
               WSET = .TRUE.
               EXTRAP = .FALSE.
            END IF
         END IF
C
C        ---------------------------------------------------------------
C        Check the termination criteria.  WSET will always be .TRUE..
C        ---------------------------------------------------------------
         TOL = TOLABS + TOLREL*ALFBST
         TRUEA = ALFBST + A
         TRUEB = ALFBST + B
C
         FOUND = ABS(GBEST) .LE. TARGTG
         QUITF = NUMF .GE. MAXF
         QUITI = B - A .LE. TOL + TOL
C
         IF (QUITI .AND. .NOT. MOVED) THEN
C
C           The interval of uncertainty appears to be small enough,
C           but no better point has been found.  Check that changing
C           ALFA by B-A changes f by less than EPSAF.
C
            TOL = TOL/TEN
            TOLABS = TOL
            QUITI = ABS(FW) .LE. EPSAF .OR. TOL .LE. TOLTNY
         END IF
C
         DONE = QUITF .OR. QUITI .OR. FOUND
C
         IF (DEBUG) THEN
            WRITE (REC,FMT=99997) TRUEA, TRUEB, B - A, TOL, NSAMEA,
     *        NSAMEB, NUMF, BRAKTD, EXTRAP, CLOSEF, IMPRVD, FOUND,
     *        QUITI, ALFBST, FBEST, GBEST, ALFBST + XW, FW, GW
            CALL X04BAY(NOUT,7,REC)
         END IF
C
C        ---------------------------------------------------------------
C        Proceed with the computation of a trial steplength.
C        The choices are...
C        1. Parabolic fit using derivatives only, if the f values are
C           close.
C        2. Cubic fit for a minimizer, using both f and f'.
C        3. Damped cubic or parabolic fit if the regular fit appears to
C           be consistently overestimating the distance to a minimizer.
C        4. Bisection, geometric bisection, or a step of  TOL  if
C           choices 2 or 3 are unsatisfactory.
C        ---------------------------------------------------------------
         IF ( .NOT. DONE) THEN
            XMIDPT = HALF*(A+B)
            S = ZERO
            Q = ZERO
C
            IF (CLOSEF) THEN
C              ---------------------------------------------------------
C              Fit a parabola to the two best gradient values.
C              ---------------------------------------------------------
               S = GBEST
               Q = GBEST - GW
               IF (DEBUG) THEN
                  WRITE (REC,FMT=99995)
                  CALL X04BAF(NOUT,REC(1))
               END IF
            ELSE
C              ---------------------------------------------------------
C              Fit cubic through  FBEST  and  FW.
C              ---------------------------------------------------------
               IF (DEBUG) THEN
                  WRITE (REC,FMT=99996)
                  CALL X04BAF(NOUT,REC(1))
               END IF
               FITOK = .TRUE.
               R = THREE*(FBEST-FW)/XW + GBEST + GW
               ABSR = ABS(R)
               S = SQRT(ABS(GBEST))*SQRT(ABS(GW))
C
C              Compute  Q =  the square root of  R*R - GBEST*GW.
C              The method avoids unnecessary underflow and overflow.
C
               IF ((GW.LT.ZERO .AND. GBEST.GT.ZERO)
     *             .OR. (GW.GT.ZERO .AND. GBEST.LT.ZERO)) THEN
                  SCALE = ABSR + S
                  IF (SCALE.EQ.ZERO) THEN
                     Q = ZERO
                  ELSE
                     Q = SCALE*SQRT((ABSR/SCALE)**2+(S/SCALE)**2)
                  END IF
               ELSE IF (ABSR.GE.S) THEN
                  Q = SQRT(ABSR+S)*SQRT(ABSR-S)
               ELSE
                  FITOK = .FALSE.
               END IF
C
               IF (FITOK) THEN
C
C                 Compute a minimizer of the fitted cubic.
C
                  IF (XW.LT.ZERO) Q = -Q
                  S = GBEST - R - Q
                  Q = GBEST - GW - Q - Q
               END IF
            END IF
C           ------------------------------------------------------------
C           Construct an artificial interval  (ARTIFA, ARTIFB)  in which
C           the new estimate of a minimizer must lie.  Set a default
C           value of XTRY that will be used if the polynomial fit fails.
C           ------------------------------------------------------------
            ARTIFA = A
            ARTIFB = B
            IF ( .NOT. BRAKTD) THEN
C
C              A minimizer has not been bracketed.  Set an artificial
C              upper bound by expanding the interval  XW  by a suitable
C              FACTOR.
C
               XTRY = -FACTOR*XW
               ARTIFB = XTRY
               IF (ALFBST+XTRY.LT.ALFMAX) FACTOR = FIVE*FACTOR
C
            ELSE IF (EXTRAP) THEN
C
C              The points are configured for an extrapolation.
C              Set a default value of  XTRY  in the interval  (A, B)
C              that will be used if the polynomial fit is rejected.  In
C              the following,  DTRY  and  DAUX  denote the lengths of
C              the intervals  (A, B)  and  (0, XW)  (or  (XW, 0),  if
C              appropriate).  The value of  XTRY is the point at which
C              the exponents of  DTRY  and  DAUX  are approximately
C              bisected.
C
               DAUX = ABS(XW)
               DTRY = B - A
               IF (DAUX.GE.DTRY) THEN
                  XTRY = FIVE*DTRY*(POINT1+DTRY/DAUX)/ELEVEN
               ELSE
                  XTRY = HALF*SQRT(DAUX)*SQRT(DTRY)
               END IF
               IF (XW.GT.ZERO) XTRY = -XTRY
               IF (DEBUG) THEN
                  WRITE (REC,FMT=99993) XTRY, DAUX, DTRY
                  CALL X04BAF(NOUT,REC(1))
               END IF
C
C              Reset the artificial bounds.  If the point computed by
C              extrapolation is rejected,  XTRY will remain at the
C              relevant artificial bound.
C
               IF (XTRY.LE.ZERO) ARTIFA = XTRY
               IF (XTRY.GT.ZERO) ARTIFB = XTRY
            ELSE
C
C              The points are configured for an interpolation.  The
C              default value XTRY bisects the interval of uncertainty.
C              the artificial interval is just (A, B).
C
               XTRY = XMIDPT
               IF (DEBUG) THEN
                  WRITE (REC,FMT=99994) XTRY
                  CALL X04BAF(NOUT,REC(1))
               END IF
               IF (NSAMEA.GE.3 .OR. NSAMEB.GE.3) THEN
C
C                 If the interpolation appears to be overestimating the
C                 distance to a minimizer,  damp the interpolation.
C
                  FACTOR = FACTOR/FIVE
                  S = FACTOR*S
               ELSE
                  FACTOR = ONE
               END IF
            END IF
C           ------------------------------------------------------------
C           The polynomial fits give  (S/Q)*XW  as the new step.
C           Reject this step if it lies outside  (ARTIFA, ARTIFB).
C           ------------------------------------------------------------
            IF (Q.NE.ZERO) THEN
               IF (Q.LT.ZERO) S = -S
               IF (Q.LT.ZERO) Q = -Q
               IF (S*XW.GE.Q*ARTIFA .AND. S*XW.LE.Q*ARTIFB) THEN
C
C                 Accept the polynomial fit.
C
                  IF (ABS(S*XW).GE.Q*TOL) THEN
                     XTRY = (S/Q)*XW
                  ELSE
                     XTRY = ZERO
                  END IF
                  IF (DEBUG) THEN
                     WRITE (REC,FMT=99992) XTRY
                     CALL X04BAF(NOUT,REC(1))
                  END IF
               END IF
            END IF
         END IF
      END IF
C
C     ==================================================================
C
      IF ( .NOT. DONE) THEN
         ALFA = ALFBST + XTRY
         IF (BRAKTD .OR. ALFA.LT.ALFMAX-TOLMAX) THEN
C
C           The function must not be evaluated too close to A or B.
C           (It has already been evaluated at both those points.)
C
            IF (XTRY.LE.A+TOL .OR. XTRY.GE.B-TOL) THEN
               IF (HALF*(A+B).LE.ZERO) THEN
                  XTRY = -TOL
               ELSE
                  XTRY = TOL
               END IF
               ALFA = ALFBST + XTRY
            END IF
         ELSE
C
C           The step is close to, or larger than ALFMAX, replace it by
C           ALFMAX to force evaluation of  f  at the boundary.
C
            BRAKTD = .TRUE.
            XTRY = ALFMAX - ALFBST
            ALFA = ALFMAX
         END IF
      END IF
C
C     ------------------------------------------------------------------
C     Exit.
C     ------------------------------------------------------------------
      IF (DONE) THEN
         IF (BADFUN) THEN
            INFORM = 8
         ELSE IF (FOUND) THEN
            IF (ALFBST.LT.ALFMAX) THEN
               INFORM = 1
            ELSE
               INFORM = 2
            END IF
         ELSE IF (MOVED) THEN
            INFORM = 3
         ELSE IF (QUITF) THEN
            INFORM = 7
         ELSE IF (CRAMPD) THEN
            INFORM = 4
         ELSE
            INFORM = 6
         END IF
      END IF
C
      IF (DEBUG) THEN
         WRITE (REC,FMT=99991)
         CALL X04BAY(NOUT,2,REC)
      END IF
      RETURN
C
C
C     End of E04UCK. (SRCHC)
C
99999 FORMAT (/'     G0  TOLABS  ALFMAX        ',1P,2D22.14,D16.8,/' T',
     *       'ARGTG  TOLREL   EPSAF        ',1P,2D22.14,D16.8,/' CRAMP',
     *       'D                        ',L3)
99998 FORMAT (/' ALFA    FTRY    GTRY          ',1P,2D22.14,D16.8)
99997 FORMAT (/' A       B       B - A   TOL   ',1P,2D22.14,2D16.8,
     *       /' NSAMEA  NSAMEB  NUMF          ',3I3,/' BRAKTD  EXTRAP ',
     *       ' CLOSEF  IMPRVD',4L3,/' FOUND   QUITI                 ',
     *       2L3,/' ALFBST  FBEST   GBEST         ',1P,3D22.14,/' ALFA',
     *       'W   FW      GW            ',1P,3D22.14)
99996 FORMAT (' Cubic.   ')
99995 FORMAT (' Parabola.')
99994 FORMAT (' Bisection.              XMIDPT',1P,D22.14)
99993 FORMAT (' Geo. bisection. XTRY,DAUX,DTRY',1P,3D22.14)
99992 FORMAT (' Polynomial fit accepted.  XTRY',1P,D22.14)
99991 FORMAT (' ----------------------------------------------------',/)
      END
