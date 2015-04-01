      SUBROUTINE E04MFS(FIRSTV,N,NCLIN,ISTATE,BIGALF,BIGBND,PNORM,
     *                  HITLOW,MOVE,ONBND,UNBNDD,ALFA,ALFAP,JHIT,ANORM,
     *                  AP,AX,BL,BU,FEATOL,FEATLU,P,X)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1563 (JUN 1995).
C
C     ******************************************************************
C     E04MFS  finds a step ALFA such that the point x + ALFA*p reaches
C     one of the linear constraints (including bounds).
C
C     In this version of E04MFS, when x is infeasible, the number of
C     infeasibilities will never increase.  If the number stays the
C     same, the sum of infeasibilities will decrease.  If the number
C     decreases by one or more,  the sum of infeasibilities will usually
C     decrease also, but occasionally it will increase after the step
C     ALFA  is taken.  (Convergence is still assured because the number
C     has decreased.)
C
C     Three possible steps are computed as follows:
C
C     alfaf = the maximum step that can be taken without violating
C              one of the constraints that are currently satisfied.
C
C     ALFAI = reaches a linear constraint that is currently violated.
C              Usually this will be the furthest such constraint along
C              p, subject to the angle between the constraint normal and
C              p being reasonably close to the maximum value among
C              infeasible constraints,  but if FIRSTV = .true. it will
C              be the first one along p.  The latter case applies only
C              when the problem has been determined to be infeasible,
C              and the sum of infeasibilities are being minimized.
C              (ALFAI is not defined when x is feasible.)
C
C     ALFAI is needed occasionally when infeasible, to prevent
C     going unnecessarily far when alfaf is quite large.  It will
C     always come into effect when x is about to become feasible.
C     (The sum of infeasibilities will decrease initially as ALFA
C     increases from zero, but may start increasing for larger steps.
C     Choosing a large ALFAI allows several elements of  x  to
C     become feasible at the same time.
C
C     In the end, we take  ALFA = alfaf  if x is feasible, or if
C     ALFAI > ALFAP (where  ALFAP  is the perturbed step from pass 1).
C     Otherwise,  we take  ALFA = ALFAI.
C
C     Input parameters
C     ----------------
C     BIGALF defines what should be treated as an unbounded step.
C     BIGBND provides insurance for detecting unboundedness.
C            If ALFA reaches a bound as large as BIGBND, it is
C            classed as an unbounded step.
C     FEATOL is the array of current feasibility tolerances used by
C            E04MFH.  Typically in the range 0.5*TOLX to 0.99*TOLX,
C            where TOLX is the FEATOL specified by the user.
C     TOLINC (in common) is used to determine STEPMN (see below),
C            the minimum positive step.
C     ISTATE is set as follows:
C            ISTATE(j) = -2  if a'x .lt. BL - FEATOL
C                      = -1  if a'x .gt. BU + FEATOL
C                      =  0  if a'x is not in the working set
C                      =  1  if a'x is in the working set at BL
C                      =  2  if a'x is in the working set at BU
C                      =  3  if a'x is in the working set (an equality)
C                      =  4  if x(j) is temporarily fixed.
C            values -2 and -1 do not occur once feasible.
C     BL     the lower bounds on the variables.
C     BU     the upper bounds on ditto.
C     x      the values of       ditto.
C     p      the search direction.
C
C
C     Output Parameters
C     -----------------
C     HITLOW  = true  if a lower bound restricted ALFA.
C             = false otherwise.
C     MOVE    = true  if  EXACT ge STEPMN  (defined at end of code).
C     ONBND   = true  if  ALFA = EXACT.  This means that the step  ALFA
C                     moves x  exactly onto one of its constraints,
C                     namely  bound.
C             = false if the exact step would be too small
C                     ( EXACT .lt. STEPMN ).
C               (with these definitions,  MOVE = ONBND).
C     UNBNDD  = true  if ALFA = BIGALF.  JHIT may possibly be zero.
C               The parameters HITLOW, MOVE, ONBND, BOUND and EXACT
C               should not be used.
C     JHIT    = the index (if any) such that constraint JHIT reaches
C               a bound.
C     BOUND   = the bound value BL(JHIT) or BU(JHIT) corresponding
C               to HITLOW.
C     EXACT   = the step that would take constraint JHIT exactly onto
C               BOUND.
C     ALFA    = an allowable, positive step.
C               if UNBNDD is true,  ALFA = STEPMX.
C               otherwise,          ALFA = max( STEPMN, EXACT ).
C
C
C     E04MFS is based on MINOS 5.2 routine M5CHZR, which implements the
C     expand procedure to deal with degeneracy. The step alfaf is
C     chosen as in the two-pass approach of Paula Harris (1973), except
C     that this version insists on returning a positive step, ALFA.
C     Two features make this possible:
C
C        1. FEATOL increases slightly each iteration.
C
C        2. The blocking constraint, when added to the working set,
C           retains the value Ax(JHIT) + ALFA * Ap(JHIT),
C           even if this is not exactly on the blocking bound.
C
C     For infeasible variables moving towards their bound, we require
C     the rate of change of the chosen constraint to be at least GAMMA
C     times as large as the biggest available.  This still gives us
C     freedom in pass 2.
C     GAMMA = 0.1 and 0.01 seemed to inhibit phase 1 somewhat.
C     GAMMA = 0.001 seems to be safe.
C
C
C     Original version written by PEG,  19-April 1988.
C     This version of  E04MFS  dated   6-Jul-1988.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      DOUBLE PRECISION  GAMMA
      PARAMETER         (GAMMA=1.0D-3)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, ALFAP, BIGALF, BIGBND, PNORM
      INTEGER           JHIT, N, NCLIN
      LOGICAL           FIRSTV, HITLOW, MOVE, ONBND, UNBNDD
C     .. Array Arguments ..
      DOUBLE PRECISION  ANORM(*), AP(*), AX(*), BL(N+NCLIN),
     *                  BU(N+NCLIN), FEATLU(N+NCLIN), FEATOL(N+NCLIN),
     *                  P(N), X(N)
      INTEGER           ISTATE(N+NCLIN)
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPSPT3, EPSPT5, EPSPT8, EPSPT9, TOLINC, TOLX0
      INTEGER           IDEGEN, ITNFIX, KDEGEN, NDEGEN
C     .. Arrays in Common ..
      INTEGER           NFIX(2)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALFAI, ATP, ATPABS, ATPMXF, ATPMXI, ATPSCD, ATX,
     *                  BIGLOW, BIGUPP, BOUND, DELTA, EXACT, RES,
     *                  STEPMN, TOLPIV
      INTEGER           I, J, JHITF, JHITI, JS
      LOGICAL           BLOCKF, BLOCKI
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Common blocks ..
      COMMON            /CE04MF/TOLX0, TOLINC, IDEGEN, KDEGEN, NDEGEN,
     *                  ITNFIX, NFIX
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
C     .. Executable Statements ..
C
C     TOLPIV is a tolerance to exclude negligible elements of a'p.
C
      BIGLOW = -BIGBND
      BIGUPP = BIGBND
C
      TOLPIV = EPSPT9*PNORM
C
C     ------------------------------------------------------------------
C     First pass -- find steps to perturbed constraints, so that
C     ALFAP will be slightly larger than the true step.
C     In degenerate cases, this strategy gives us some freedom in the
C     second pass.  The general idea follows that described by P.M.J.
C     Harris, p.21 of Mathematical Programming 5, 1 (1973), 1--28.
C     ------------------------------------------------------------------
      ATPMXI = ZERO
      ALFAP = BIGALF
C
      DO 20 J = 1, N + NCLIN
         JS = ISTATE(J)
C
         IF (JS.LE.0) THEN
            DELTA = FEATOL(J)
C
            IF (J.LE.N) THEN
               ATX = X(J)
               ATP = P(J)
               ATPABS = ABS(ATP)
               ATPSCD = ATPABS
            ELSE
               I = J - N
               ATX = AX(I)
               ATP = AP(I)
               ATPABS = ABS(ATP)
               ATPSCD = ATPABS/(ONE+ANORM(I))
            END IF
C
            IF (ATPSCD.LE.TOLPIV) THEN
C              ---------------------------------------------------------
C              This constraint appears to be constant along p.  It is
C              not used to compute the step.  Give the residual a value
C              that can be spotted in the debug output.
C              ---------------------------------------------------------
               RES = -ONE
C
            ELSE IF (ATP.LE.ZERO .AND. JS.NE.-2) THEN
C              ---------------------------------------------------------
C              a'x  is decreasing and the lower bound is not violated.
C              ---------------------------------------------------------
C              test for smaller ALFAP.
C              if the upper bound is violated. test for bigger ATP.
C
               IF (BL(J).GT.BIGLOW) THEN
                  RES = ATX - BL(J) + DELTA
C
                  IF (RES.LT.ALFAP*ATPABS) ALFAP = RES/ATPABS
               END IF
C
               IF (JS.EQ.-1) ATPMXI = MAX(ATPMXI,ATPSCD)
C
            ELSE IF (ATP.GT.ZERO .AND. JS.NE.-1) THEN
C              ---------------------------------------------------------
C              a'x  is increasing and the upper bound is not violated.
C              ---------------------------------------------------------
C              test for smaller ALFAP.
C              if the lower bound is violated. test for bigger ATP.
C
               IF (BU(J).LT.BIGUPP) THEN
                  RES = BU(J) - ATX + DELTA
C
                  IF (RES.LT.ALFAP*ATP) ALFAP = RES/ATP
               END IF
C
               IF (JS.EQ.-2) ATPMXI = MAX(ATPMXI,ATPSCD)
            END IF
C
         END IF
   20 CONTINUE
C
C     ------------------------------------------------------------------
C     Second pass.
C     For feasible variables, recompute steps without perturbation.
C     amongst constraints that are closer than ALFAP, choose the one
C     That makes the largest angle with the search direction.
C     For infeasible variables, find the largest step subject to a'p
C     being no smaller than GAMMA * max(a'p).
C     ------------------------------------------------------------------
      IF (FIRSTV) THEN
         ALFAI = BIGALF
      ELSE
         ALFAI = ZERO
      END IF
C
      ATPMXF = ZERO
      ATPMXI = GAMMA*ATPMXI
      JHITF = 0
      JHITI = 0
C
      DO 40 J = 1, N + NCLIN
         JS = ISTATE(J)
C
         IF (JS.LE.0) THEN
C
            IF (J.LE.N) THEN
               ATX = X(J)
               ATP = P(J)
               ATPABS = ABS(ATP)
               ATPSCD = ATPABS
            ELSE
               I = J - N
               ATX = AX(I)
               ATP = AP(I)
               ATPABS = ABS(ATP)
               ATPSCD = ATPABS/(ONE+ANORM(I))
            END IF
C
            IF (ATPSCD.LE.TOLPIV) THEN
C              ---------------------------------------------------------
C              this constraint appears to be constant along p.  it is
C              not used to compute the step.  give the residual a value
C              that can be spotted in the debug output.
C              ---------------------------------------------------------
               RES = -ONE
C
            ELSE IF (ATP.LE.ZERO .AND. JS.NE.-2) THEN
C              ---------------------------------------------------------
C              a'x  is decreasing.
C              ---------------------------------------------------------
C              test for bigger a'p if the lower bound is satisfied.
C              test for smaller alfaf.
C
               IF (ATPSCD.GT.ATPMXF) THEN
C
                  IF (BL(J).GT.BIGLOW) THEN
                     RES = ATX - BL(J)
C
                     IF (RES.LE.ALFAP*ATPABS) THEN
                        ATPMXF = ATPSCD
                        JHITF = J
                     END IF
                  END IF
               END IF
C
               IF (JS.EQ.-1) THEN
C
C                 the upper bound is violated.
C                 test for bigger or smaller ALFAI,  depending on the
C                 value of FIRSTV.
C
                  IF (FIRSTV) THEN
                     RES = ATX - BU(J)
C
                     IF (RES.LE.ALFAI*ATPABS) THEN
                        ALFAI = RES/ATPABS
                        JHITI = J
                     END IF
C
                  ELSE IF (ATPSCD.GE.ATPMXI) THEN
                     RES = ATX - BU(J)
C
                     IF (RES.GT.ALFAI*ATPABS) THEN
                        ALFAI = RES/ATPABS
                        JHITI = J
                     END IF
                  END IF
               END IF
C
            ELSE IF (ATP.GT.ZERO .AND. JS.NE.-1) THEN
C              ---------------------------------------------------------
C              a'x  is increasing and the upper bound is not violated.
C              ---------------------------------------------------------
C              test for smaller ALFAP.
C
               IF (ATPSCD.GT.ATPMXF) THEN
C
                  IF (BU(J).LT.BIGUPP) THEN
                     RES = BU(J) - ATX
C
                     IF (RES.LE.ALFAP*ATP) THEN
                        ATPMXF = ATPSCD
                        JHITF = J
                     END IF
                  END IF
               END IF
C
               IF (JS.EQ.-2) THEN
C
C                 the lower bound is violated.
C                 test for bigger or smaller ALFAI,  depending on the
C                 value of FIRSTV.
C
                  IF (FIRSTV) THEN
                     RES = BL(J) - ATX
C
                     IF (RES.LE.ALFAI*ATP) THEN
                        ALFAI = RES/ATP
                        JHITI = J
                     END IF
                  ELSE IF (ATPSCD.GE.ATPMXI) THEN
                     RES = BL(J) - ATX
C
                     IF (RES.GT.ALFAI*ATP) THEN
                        ALFAI = RES/ATP
                        JHITI = J
                     END IF
                  END IF
               END IF
            END IF
C
         END IF
   40 CONTINUE
C
C     ------------------------------------------------------------------
C     See if a feasible and/or infeasible constraint blocks.
C     ------------------------------------------------------------------
      BLOCKF = JHITF .GT. 0
      BLOCKI = JHITI .GT. 0
      UNBNDD = .NOT. (BLOCKF .OR. BLOCKI)
C
      IF (UNBNDD) GO TO 60
C
      IF (BLOCKF) THEN
C        ---------------------------------------------------------------
C        A constraint is hit which is currently feasible.
C        The corresponding step alfaf is not used, so no need to get it,
C        but we know that alfaf .le. ALFAP, the step from pass 1.
C        ---------------------------------------------------------------
         JHIT = JHITF
         IF (JHIT.LE.N) THEN
            ATP = P(JHIT)
         ELSE
            ATP = AP(JHIT-N)
         END IF
         HITLOW = ATP .LT. ZERO
      END IF
C
C     If there is a choice between alfaf and ALFAI, it is probably best
C     to take ALFAI.  However, we can't if ALFAI is bigger than ALFAP.
C
      IF (BLOCKI .AND. ALFAI.LE.ALFAP) THEN
C        ---------------------------------------------------------------
C        An infeasible variable reaches its violated bound.
C        ---------------------------------------------------------------
         JHIT = JHITI
         IF (JHIT.LE.N) THEN
            ATP = P(JHIT)
         ELSE
            ATP = AP(JHIT-N)
         END IF
         HITLOW = ATP .GT. ZERO
      END IF
C
      IF (JHIT.LE.N) THEN
         ATX = X(JHIT)
      ELSE
         ATX = AX(JHIT-N)
      END IF
C
C     ------------------------------------------------------------------
C     Try to step exactly onto bound, but make sure the exact step
C     is sufficiently positive.  (Exact will be alfaf or ALFAI.)
C     Since FEATOL increases by  TOLINC  each iteration, we know that
C     a step as large as  STEPMN  (below) will not cause any feasible
C     variables to become infeasible (where feasibility is measured
C     by the current FEATOL).
C     ------------------------------------------------------------------
      IF (HITLOW) THEN
         BOUND = BL(JHIT)
      ELSE
         BOUND = BU(JHIT)
      END IF
C
      UNBNDD = ABS(BOUND) .GE. BIGBND
      IF (UNBNDD) GO TO 60
C
      STEPMN = TOLINC*FEATLU(JHIT)/ABS(ATP)
      EXACT = (BOUND-ATX)/ATP
      ALFA = MAX(STEPMN,EXACT)
      ONBND = ALFA .EQ. EXACT
      MOVE = EXACT .GE. STEPMN
      IF ( .NOT. MOVE) NDEGEN = NDEGEN + 1
C
      RETURN
C     ------------------------------------------------------------------
C     Unbounded.
C     ------------------------------------------------------------------
   60 ALFA = BIGALF
      MOVE = .TRUE.
      ONBND = .FALSE.
C
      RETURN
C
C     End of  E04MFS.  (CMCHZR)
C
      END
