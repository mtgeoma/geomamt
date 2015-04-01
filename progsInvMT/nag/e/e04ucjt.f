      SUBROUTINE E04UCJ(FIRST,DEBUG,DONE,IMPRVD,INFORM,MAXF,NUMF,NOUT,
     *                  ALFMAX,ALFSML,EPSAF,G0,TARGTG,FTRY,TOLABS,
     *                  TOLREL,TOLTNY,ALFA,ALFBST,FBEST)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     MARK 16A REVISED. IER-997 (JUN 1993).
C
C     ==================================================================
C     E04UCJ  finds a sequence of improving estimates of a minimizer of
C     the univariate function f(alpha) in the interval (0,ALFMAX].
C     f(alpha) is a smooth function such that  f(0) = 0  and  f'(0) < 0.
C     E04UCJ  requires  f(alpha) (but not f'(alpha)) to be evaluated
C     in the interval.  New estimates of a minimizer are computed using
C     safeguarded quadratic interpolation.
C
C     Reverse communication is used to allow the calling program to
C     evaluate f.  Some of the parameters must be set or tested by the
C     calling program.  The remainder would ordinarily be local
C     variables.
C
C     Input parameters (relevant to the calling program)
C     --------------------------------------------------
C
C     FIRST         must be .TRUE. on the first entry.
C                   It is subsequently altered by E04UCJ.
C
C     DEBUG         specifies whether detailed output is wanted.
C
C     MAXF          is an upper limit on the number of times E04UCJ is
C                   to be entered consecutively with DONE = .FALSE.
C                   (following an initial entry with FIRST = .TRUE.).
C
C     ALFA          is the first estimate of a minimizer.  ALFA is
C                   subsequently altered by E04UCJ (see below).
C
C     ALFMAX        is the upper limit of the interval to be searched.
C
C     ALFSML        is intended to prevent inefficiency when a minimizer
C                   is very small, for cases where the calling program
C                   would prefer to redefine f'(ALFA).  ALFSML is
C                   allowed to be zero.  Early termination will occur if
C                   E04UCJ determines that a minimizer lies somewhere in
C                   the interval [0, ALFSML) (but not if ALFMAX is
C                   smaller that ALFSML).
C
C     EPSAF         is an estimate of the absolute precision in the
C                   computed value of f(0).
C
C     FTRY          the value of f at the new point
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
C                   calling program (e.g., arrays) before paying
C                   attention to the variable DONE.
C
C     DONE = .FALSE.  means the calling program should evaluate FTRY
C                   for the new trial step ALFA, and reenter E04UCJ.
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
C                   INFORM = 5 means that the search was terminated
C                              because of ALFSML (see above).
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
C                              ALFMAX le TOLTNY  or  G0 ge zero.
C                              No function evaluations were made.
C
C     NUMF          counts the number of times E04UCJ has been entered
C                   consecutively with DONE = .FALSE. (i.e., with a new
C                   function value FTRY).
C
C     ALFA          is the point at which the next function FTRY must
C                   be computed.
C
C     ALFBST        should be accepted by the calling program as the
C                   approximate minimizer, whenever E04UCJ returns
C                   INFORM = 1, 2 or 3.
C
C     FBEST         will be the corresponding value of f.
C
C     The following parameters retain information between entries
C     -----------------------------------------------------------
C
C     BRAKTD        is .FALSE. if f has not been evaluated at the far
C                   end of the interval of uncertainty.  In this case,
C                   the point B will be at ALFMAX + TOL(ALFMAX).
C
C     CRAMPD        is .TRUE. if ALFMAX is very small (le TOLABS).  If
C                   the search fails, this indicates that a zero step
C                   should be taken.
C
C     EXTRAP        is .TRUE. if ALFBST has MOVED at least once and XV
C                   lies outside the interval of uncertainty.  In this
C                   case, extra safeguards are applied to allow for
C                   instability in the polynomial fit.
C
C     MOVED         is .TRUE. if a better point has been found, i.e.,
C                   ALFBST gt 0.
C
C     VSET          records whether a third-best point has been defined.
C
C     WSET          records whether a second-best point has been
C                   defined.  It will always be .TRUE. by the time the
C                   convergence test is applied.
C
C     NSAMEA        is the number of consecutive times that the
C                   left-hand end point of the interval of uncertainty
C                   has remained the same.
C
C     NSAMEB        similarly for the right-hand end.
C
C     A, B, ALFBST  define the current interval of uncertainty.
C                   A minimizer lies somewhere in the  interval
C                   [ALFBST + A, ALFBST + B].
C
C     ALFBST        is the best point so far.  It lies strictly within
C                   [ATRUE,BTRUE]  (except when ALFBST has not been
C                   MOVED, in which case it lies at the left-hand end
C                   point).  Hence we have A .le. 0 and B .gt. 0.
C
C     FBEST         is the value of f at the point ALFBST.
C
C     FA            is the value of f at the point ALFBST + A.
C
C     FACTOR        controls the rate at which extrapolated estimates of
C                   ALFA  may expand into the interval of uncertainty.
C                   FACTOR is not used if a minimizer has been bracketed
C                   (i.e., when the variable BRAKTD is .TRUE.).
C
C     FV, FW        are the values of f at the points ALFBST + XV  and
C                   ALFBST + XW.  They are not defined until  VSET  or
C                   WSET  are .TRUE..
C
C     XTRY          is the trial point within the shifted interval
C                   (A, B).  The new trial function value must be
C                   computed at the point ALFA = ALFBST + XTRY.
C
C     XV            is such that ALFBST + XV is the third-best point.
C                   It is not defined until VSET is .TRUE..
C
C     XW            is such that ALFBST + XW is the second-best point.
C                   It is not defined until WSET is .TRUE..  In some
C                   cases,  XW will replace a previous XW that has a
C                   lower function but has just been excluded from
C                   (A,B).
C
C     Systems Optimization Laboratory, Stanford University, California.
C     Original version February 1982.  Rev. May 1983.
C     Original F77 version 22-August-1985.
C     This version of E04UCJ dated  24-Oct-91.
C     ==================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, POINT1, HALF
      PARAMETER         (ZERO=0.0D+0,POINT1=0.1D+0,HALF=0.5D+0)
      DOUBLE PRECISION  ONE, TWO, FIVE
      PARAMETER         (ONE=1.0D+0,TWO=2.0D+0,FIVE=5.0D+0)
      DOUBLE PRECISION  TEN, ELEVEN
      PARAMETER         (TEN=1.0D+1,ELEVEN=1.1D+1)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, ALFBST, ALFMAX, ALFSML, EPSAF, FBEST,
     *                  FTRY, G0, TARGTG, TOLABS, TOLREL, TOLTNY
      INTEGER           INFORM, MAXF, NOUT, NUMF
      LOGICAL           DEBUG, DONE, FIRST, IMPRVD
C     .. Local Scalars ..
      DOUBLE PRECISION  A, ARTIFA, ARTIFB, B, DAUX, DTRY, ENDPNT, FA,
     *                  FACTOR, FV, FW, GV, GW, Q, S, TOL, TOLMAX,
     *                  TRUEA, TRUEB, XMIDPT, XTRY, XV, XW
      INTEGER           NSAMEA, NSAMEB
      LOGICAL           BADFUN, BRAKTD, CLOSEF, CRAMPD, EXTRAP, FOUND,
     *                  MOVED, QUITF, QUITFZ, QUITI, QUITS, SETXV, VSET,
     *                  WSET, XINXW
C     .. Local Arrays ..
      CHARACTER*120     REC(7)
C     .. External Subroutines ..
      EXTERNAL          X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Save statement ..
      SAVE              BRAKTD, CRAMPD, EXTRAP, MOVED, VSET, WSET,
     *                  NSAMEA, NSAMEB, A, B, FA, FACTOR, XTRY, XW, FW,
     *                  XV, FV, TOLMAX
C     .. Executable Statements ..
C
C     ------------------------------------------------------------------
C     Local variables
C     ===============
C
C     CLOSEF     is .TRUE. if the worst function FV is within EPSAF of
C                FBEST (up or down).
C
C     FOUND      is .TRUE. if the sufficient decrease conditions holds
C                at ALFBST.
C
C     QUITF      is .TRUE. when  MAXF  function calls have been made.
C
C     QUITFZ     is .TRUE. when the three best function values are
C                within EPSAF of each other, and the new point satisfies
C                FBEST le FTRY le FBEST+EPSAF.
C
C     QUITI      is .TRUE. when the interval of uncertainty is less than
C                2*TOL.
C
C     QUITS      is .TRUE. as soon as ALFA is too small to be useful;
C                i.e., BTRUE le ALFSML.
C
C     XINXW      is .TRUE. if XTRY is in (XW,0) or (0,XW).
C     ------------------------------------------------------------------
C
      IMPRVD = .FALSE.
      BADFUN = .FALSE.
      QUITF = .FALSE.
      QUITFZ = .FALSE.
      QUITS = .FALSE.
      QUITI = .FALSE.
C
      IF (FIRST) THEN
C        ---------------------------------------------------------------
C        First entry.  Initialize various quantities, check input data
C        and prepare to evaluate the function at the initial step ALFA.
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
            VSET = .FALSE.
            WSET = .FALSE.
            NSAMEA = 0
            NSAMEB = 0
C
            TOLMAX = TOLREL*ALFMAX + TOLABS
            A = ZERO
            B = ALFMAX + TOLMAX
            FA = ZERO
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
C        Subsequent entries.  The function has just been evaluated at
C        ALFA = ALFBST + XTRY,  giving FTRY.
C        ---------------------------------------------------------------
         IF (DEBUG) THEN
            WRITE (REC,FMT=99998) ALFA, FTRY
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
C        Check if XTRY is in the interval (XW,0) or (0,XW).
C
         IF (WSET) THEN
            XINXW = ZERO .LT. XTRY .AND. XTRY .LE. XW .OR. XW .LE.
     *              XTRY .AND. XTRY .LT. ZERO
         ELSE
            XINXW = .FALSE.
         END IF
C
         IMPRVD = FTRY .LT. FBEST
         IF (VSET) THEN
            CLOSEF = ABS(FBEST-FV) .LE. EPSAF
         ELSE
            CLOSEF = .FALSE.
         END IF
C
         IF (IMPRVD) THEN
C
C           We seem to have an improvement.  The new point becomes the
C           origin and other points are shifted accordingly.
C
            IF (WSET) THEN
               XV = XW - XTRY
               FV = FW
               VSET = .TRUE.
            END IF
C
            XW = ZERO - XTRY
            FW = FBEST
            WSET = .TRUE.
            FBEST = FTRY
            ALFBST = ALFA
            MOVED = .TRUE.
C
            A = A - XTRY
            B = B - XTRY
            EXTRAP = .NOT. XINXW
C
C           Decrease the length of (A,B).
C
            IF (XTRY.GE.ZERO) THEN
               A = XW
               FA = FW
               NSAMEA = 0
            ELSE
               B = XW
               NSAMEB = 0
               BRAKTD = .TRUE.
            END IF
         ELSE IF (CLOSEF .AND. FTRY-FBEST.LT.EPSAF) THEN
C
C           Quit if there has been no progress and FTRY, FBEST, FW
C           and FV are all within EPSAF of each other.
C
            QUITFZ = .TRUE.
         ELSE
C
C           The new function value is no better than the current best
C           point.  XTRY must an end point of the new (A,B).
C
            IF (XTRY.LT.ZERO) THEN
               A = XTRY
               FA = FTRY
               NSAMEA = 0
            ELSE
               B = XTRY
               NSAMEB = 0
               BRAKTD = .TRUE.
            END IF
C
C           The origin remains unchanged but XTRY may qualify as XW.
C
            IF (WSET) THEN
               IF (FTRY.LT.FW) THEN
                  XV = XW
                  FV = FW
                  VSET = .TRUE.
C
                  XW = XTRY
                  FW = FTRY
                  IF (MOVED) EXTRAP = XINXW
               ELSE IF (MOVED) THEN
                  IF (VSET) THEN
                     SETXV = FTRY .LT. FV .OR. .NOT. EXTRAP
                  ELSE
                     SETXV = .TRUE.
                  END IF
C
                  IF (SETXV) THEN
                     IF (VSET .AND. XINXW) THEN
                        XW = XV
                        FW = FV
                     END IF
                     XV = XTRY
                     FV = FTRY
                     VSET = .TRUE.
                  END IF
               ELSE
                  XW = XTRY
                  FW = FTRY
               END IF
            ELSE
               XW = XTRY
               FW = FTRY
               WSET = .TRUE.
            END IF
         END IF
C
C        ---------------------------------------------------------------
C        Check the termination criteria.
C        ---------------------------------------------------------------
         TOL = TOLABS + TOLREL*ALFBST
         TRUEA = ALFBST + A
         TRUEB = ALFBST + B
C
         FOUND = MOVED .AND. ABS(FA-FBEST) .LE. -A*TARGTG
         QUITF = NUMF .GE. MAXF
         QUITI = B - A .LE. TOL + TOL
         QUITS = TRUEB .LE. ALFSML
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
         DONE = QUITF .OR. QUITFZ .OR. QUITS .OR. QUITI .OR. FOUND
C
         IF (DEBUG) THEN
            WRITE (REC,FMT=99997) TRUEA, TRUEB, B - A, TOL, NSAMEA,
     *        NSAMEB, NUMF, BRAKTD, EXTRAP, CLOSEF, IMPRVD, FOUND,
     *        QUITI, QUITFZ, QUITS, ALFBST, FBEST, ALFBST + XW, FW
            CALL X04BAY(NOUT,7,REC)
            IF (VSET) THEN
               WRITE (REC,FMT=99996) ALFBST + XV, FV
               CALL X04BAY(NOUT,2,REC)
            END IF
         END IF
C
C        ---------------------------------------------------------------
C        Proceed with the computation of an estimate of a minimizer.
C        The choices are...
C        1. Parabolic fit using function values only.
C        2. Damped parabolic fit if the regular fit appears to be
C           consistently overestimating the distance to a minimizer.
C        3. Bisection, geometric bisection, or a step of TOL if the
C           parabolic fit is unsatisfactory.
C        ---------------------------------------------------------------
         IF ( .NOT. DONE) THEN
            XMIDPT = HALF*(A+B)
            S = ZERO
            Q = ZERO
C
C           ============================================================
C           Fit a parabola.
C           ============================================================
C           See if there are two or three points for the parabolic fit.
C
            GW = (FW-FBEST)/XW
            IF (VSET .AND. MOVED) THEN
C
C              Three points available.  Use FBEST, FW and FV.
C
               GV = (FV-FBEST)/XV
               S = GV - (XV/XW)*GW
               Q = TWO*(GV-GW)
               IF (DEBUG) THEN
                  WRITE (REC,FMT=99994)
                  CALL X04BAF(NOUT,REC(1))
               END IF
            ELSE
C
C              Only two points available.  Use FBEST, FW and G0.
C
               IF (MOVED) THEN
                  S = G0 - TWO*GW
               ELSE
                  S = G0
               END IF
               Q = TWO*(G0-GW)
               IF (DEBUG) THEN
                  WRITE (REC,FMT=99995)
                  CALL X04BAF(NOUT,REC(1))
               END IF
            END IF
C
C           ------------------------------------------------------------
C           Construct an artificial interval (ARTIFA, ARTIFB) in which
C           the new estimate of the steplength must lie.  Set a default
C           value of  XTRY  that will be used if the polynomial fit is
C           rejected. In the following, the interval (A,B) is considered
C           the sum of two intervals of lengths  DTRY  and  DAUX, with
C           common end point the best point (zero).  DTRY is the length
C           of the interval into which the default XTRY will be placed
C           and ENDPNT denotes its non-zero end point.  The magnitude of
C           XTRY is computed so that the exponents of DTRY and DAUX are
C           approximately bisected.
C           ------------------------------------------------------------
            ARTIFA = A
            ARTIFB = B
            IF ( .NOT. BRAKTD) THEN
C
C              A minimizer has not yet been bracketed.
C              Set an artificial upper bound by expanding the interval
C              XW  by a suitable FACTOR.
C
               XTRY = -FACTOR*XW
               ARTIFB = XTRY
               IF (ALFBST+XTRY.LT.ALFMAX) FACTOR = FIVE*FACTOR
            ELSE IF (VSET .AND. MOVED) THEN
C
C              Three points exist in the interval of uncertainty.
C              Check if the points are configured for an extrapolation
C              or an interpolation.
C
               IF (EXTRAP) THEN
C
C                 The points are configured for an extrapolation.
C
                  IF (XW.LT.ZERO) ENDPNT = B
                  IF (XW.GT.ZERO) ENDPNT = A
               ELSE
C
C                 If the interpolation appears to be overestimating the
C                 distance to a minimizer,  damp the interpolation step.
C
                  IF (NSAMEA.GE.3 .OR. NSAMEB.GE.3) THEN
                     FACTOR = FACTOR/FIVE
                     S = FACTOR*S
                  ELSE
                     FACTOR = ONE
                  END IF
C
C                 The points are configured for an interpolation.  The
C                 artificial interval will be just (A,B).  Set ENDPNT so
C                 that XTRY lies in the larger of the intervals (A,B)
C                 and  (0,B).
C
                  IF (XMIDPT.GT.ZERO) THEN
                     ENDPNT = B
                  ELSE
                     ENDPNT = A
                  END IF
C
C                 If a bound has remained the same for three iterations,
C                 set ENDPNT so that  XTRY  is likely to replace the
C                 offending bound.
C
                  IF (NSAMEA.GE.3) ENDPNT = A
                  IF (NSAMEB.GE.3) ENDPNT = B
               END IF
C
C              Compute the default value of  XTRY.
C
               DTRY = ABS(ENDPNT)
               DAUX = B - A - DTRY
               IF (DAUX.GE.DTRY) THEN
                  XTRY = FIVE*DTRY*(POINT1+DTRY/DAUX)/ELEVEN
               ELSE
                  XTRY = HALF*SQRT(DAUX)*SQRT(DTRY)
               END IF
               IF (ENDPNT.LT.ZERO) XTRY = -XTRY
               IF (DEBUG) THEN
                  WRITE (REC,FMT=99992) XTRY, DAUX, DTRY
                  CALL X04BAF(NOUT,REC(1))
               END IF
C
C              If the points are configured for an extrapolation set the
C              artificial bounds so that the artificial interval lies
C              within (A,B).  If the polynomial fit is rejected,  XTRY
C              will remain at the relevant artificial bound.
C
               IF (EXTRAP) THEN
                  IF (XTRY.LE.ZERO) THEN
                     ARTIFA = XTRY
                  ELSE
                     ARTIFB = XTRY
                  END IF
               END IF
            ELSE
C
C              The gradient at the origin is being used for the
C              polynomial fit.  Set the default XTRY to one tenth XW.
C
               IF (EXTRAP) THEN
                  XTRY = -XW
               ELSE
                  XTRY = XW/TEN
               END IF
               IF (DEBUG) THEN
                  WRITE (REC,FMT=99993) XTRY
                  CALL X04BAF(NOUT,REC(1))
               END IF
            END IF
C
C           ------------------------------------------------------------
C           The polynomial fits give (S/Q)*XW as the new step.  Reject
C           this step if it lies outside (ARTIFA, ARTIFB).
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
                     WRITE (REC,FMT=99991) XTRY
                     CALL X04BAF(NOUT,REC(1))
                  END IF
               END IF
            END IF
         END IF
      END IF
C     ==================================================================
C
      IF ( .NOT. DONE) THEN
         ALFA = ALFBST + XTRY
         IF (BRAKTD .OR. ALFA.LT.ALFMAX-TOLMAX) THEN
C
C           The function must not be evaluated too close to A or B.
C           (It has already been evaluated at both those points.)
C
            XMIDPT = HALF*(A+B)
            IF (XTRY.LE.A+TOL .OR. XTRY.GE.B-TOL) THEN
               IF (XMIDPT.LE.ZERO) THEN
                  XTRY = -TOL
               ELSE
                  XTRY = TOL
               END IF
            END IF
C
            IF (ABS(XTRY).LT.TOL) THEN
               IF (XMIDPT.LE.ZERO) THEN
                  XTRY = -TOL
               ELSE
                  XTRY = TOL
               END IF
            END IF
            ALFA = ALFBST + XTRY
         ELSE
C
C           The step is close to or larger than ALFMAX, replace it by
C           ALFMAX to force evaluation of the function at the boundary.
C
            BRAKTD = .TRUE.
            XTRY = ALFMAX - ALFBST
            ALFA = ALFMAX
         END IF
      END IF
C     ------------------------------------------------------------------
C     Exit.
C     ------------------------------------------------------------------
      IF (DONE) THEN
         IF (BADFUN) THEN
            INFORM = 8
         ELSE IF (QUITS) THEN
            INFORM = 5
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
         WRITE (REC,FMT=99990)
         CALL X04BAY(NOUT,2,REC)
      END IF
      RETURN
C
C
C     End of  E04UCJ. (SRCHQ)
C
99999 FORMAT (/'     G0  TOLABS  ALFMAX        ',1P,2D22.14,D16.8,/' T',
     *       'ARGTG  TOLREL   EPSAF        ',1P,2D22.14,D16.8,/' CRAMP',
     *       'D                        ',L3)
99998 FORMAT (/' ALFA    FTRY                  ',1P,2D22.14)
99997 FORMAT (/' A       B       B - A   TOL   ',1P,2D22.14,2D16.8,
     *       /' NSAMEA  NSAMEB  NUMF          ',3I3,/' BRAKTD  EXTRAP ',
     *       ' CLOSEF  IMPRVD',4L3,/' FOUND   QUITI   QUITFZ  QUITS ',
     *       4L3,/' ALFBST  FBEST                 ',1P,2D22.14,/' ALFA',
     *       'W   FW                    ',1P,2D22.14)
99996 FORMAT (' ALFAV   FV                    ',1P,2D22.14,/)
99995 FORMAT (' Parabolic fit,    two points. ')
99994 FORMAT (' Parabolic fit,  three points. ')
99993 FORMAT (' Exponent reduced.  Trial point',1P,D22.14)
99992 FORMAT (' Geo. bisection. XTRY,DAUX,DTRY',1P,3D22.14)
99991 FORMAT (' Polynomial fit accepted.  XTRY',1P,D22.14)
99990 FORMAT (' ----------------------------------------------------',/)
      END
