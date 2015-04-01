      SUBROUTINE E04UCG(FIRSTV,HITLOW,ISTATE,INFORM,JADD,N,NCTOTL,
     *                  NUMINF,ALFA,PALFA,ATPHIT,BIGALF,BIGBND,PNORM,
     *                  ANORM,AP,AX,BL,BU,FEATOL,P,X)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1078 (JUL 1993).
C     MARK 17 REVISED. IER-1600 (JUN 1995).
C
C     ******************************************************************
C     E04UCG finds a step ALFA such that the point x + ALFA*P reaches
C     one of the linear constraints (including bounds).  Two possible
C     steps are defined as follows...
C
C     ALFA1   is the maximum step that can be taken without violating
C             one of the linear constraints that is currently satisfied.
C     ALFA2   reaches a linear constraint that is currently violated.
C             Usually this will be the furthest such constraint along P,
C             but if FIRSTV = .TRUE. it will be the first one along P.
C             This is used only when the problem has been determined to
C             be infeasible, and the sum of infeasibilities are being
C             minimized.  (ALFA2  is not defined if NUMINF = 0.)
C
C     ALFA will usually be the minimum of ALFA1 and ALFA2.
C     ALFA could be negative (since we allow inactive constraints
C     to be violated by as much as FEATOL).  In such cases, a
C     third possible step is computed, to find the nearest satisfied
C     constraint (perturbed by FEATOL) along the direction  - P.
C     ALFA  will be reset to this step if it is shorter.  This is the
C     only case for which the final step  ALFA  does not move X exactly
C     onto a constraint (the one denoted by JADD).
C
C     Constraints in the working set are ignored  (ISTATE(j) ge 1).
C
C     JADD    denotes which linear constraint is reached.
C
C     HITLOW  indicates whether it is the lower or upper bound that
C             has restricted ALFA.
C
C     Values of ISTATE(j)....
C
C     - 2         - 1         0           1          2         3
C     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
C
C     The values -2 and -1 do not occur once a feasible point has been
C     found.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original Fortran 66 version written  May 1980.
C     This version of  E04UCG  dated  10-June-1986.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, ATPHIT, BIGALF, BIGBND, PALFA, PNORM
      INTEGER           INFORM, JADD, N, NCTOTL, NUMINF
      LOGICAL           FIRSTV, HITLOW
C     .. Array Arguments ..
      DOUBLE PRECISION  ANORM(*), AP(*), AX(*), BL(NCTOTL), BU(NCTOTL),
     *                  FEATOL(NCTOTL), P(N), X(N)
      INTEGER           ISTATE(NCTOTL)
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPSPT3, EPSPT5, EPSPT8, EPSPT9
C     .. Local Scalars ..
      DOUBLE PRECISION  ABSATP, ALFA1, ALFA2, APMAX1, APMAX2, ATP, ATP1,
     *                  ATP2, ATX, PALFA1, PALFA2, RES, ROWNRM
      INTEGER           I, J, JADD1, JADD2, JS, JSAVE1, JSAVE2
      LOGICAL           HLOW1, HLOW2, LASTV, NEGSTP, STEP2
C     .. External Subroutines ..
      EXTERNAL          E04UCH
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Common blocks ..
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
C     .. Executable Statements ..
C
      INFORM = 0
C
C     ------------------------------------------------------------------
C     First pass -- find steps to perturbed constraints, so that
C     PALFA1 will be slightly larger than the true step, and
C     PALFA2 will be slightly smaller than it should be.
C     In degenerate cases, this strategy gives us some freedom in the
C     second pass.  The general idea follows that described by P.M.J.
C     Harris, p.21 of Mathematical Programming 5, 1 (1973), 1--28.
C     ------------------------------------------------------------------
C
      NEGSTP = .FALSE.
      CALL E04UCH(FIRSTV,NEGSTP,BIGALF,BIGBND,PNORM,JADD1,JADD2,PALFA1,
     *            PALFA2,ISTATE,N,NCTOTL,ANORM,AP,AX,BL,BU,FEATOL,P,X)
C
      JSAVE1 = JADD1
      JSAVE2 = JADD2
C
C     ------------------------------------------------------------------
C     Second pass -- recompute step-lengths without perturbation.
C     Amongst constraints that are less than the perturbed steps,
C     choose the one (of each type) that makes the largest angle
C     with the search direction.
C     ------------------------------------------------------------------
      ALFA1 = BIGALF
      ALFA2 = ZERO
      IF (FIRSTV) ALFA2 = BIGALF
C
      APMAX1 = ZERO
      APMAX2 = ZERO
      ATP1 = ZERO
      ATP2 = ZERO
      HLOW1 = .FALSE.
      HLOW2 = .FALSE.
      LASTV = .NOT. FIRSTV
C
      DO 20 J = 1, NCTOTL
         JS = ISTATE(J)
         IF (JS.LE.0) THEN
            IF (J.LE.N) THEN
               ATX = X(J)
               ATP = P(J)
               ROWNRM = ONE
            ELSE
               I = J - N
               ATX = AX(I)
               ATP = AP(I)
               ROWNRM = ANORM(I) + ONE
            END IF
C
            IF (ABS(ATP).LE.EPSPT9*ROWNRM*PNORM) THEN
C
C              This constraint appears to be constant along P.  It is
C              not used to compute the step.  Give the residual a value
C              that can be spotted in the debug output.
C
               RES = -ONE
            ELSE IF (ATP.LE.ZERO .AND. JS.NE.-2) THEN
C              ---------------------------------------------------------
C              a'x  is decreasing.
C              ---------------------------------------------------------
C              The lower bound is satisfied.  Test for smaller ALFA1.
C
               ABSATP = -ATP
               IF (BL(J).GT.(-BIGBND)) THEN
                  RES = ATX - BL(J)
                  IF (PALFA1*ABSATP.GE.RES .OR. J.EQ.JSAVE1) THEN
                     IF (APMAX1*ROWNRM*PNORM.LT.ABSATP) THEN
                        APMAX1 = ABSATP/(ROWNRM*PNORM)
                        ALFA1 = RES/ABSATP
                        JADD1 = J
                        ATP1 = ATP
                        HLOW1 = .TRUE.
                     END IF
                  END IF
               END IF
C
               IF (JS.EQ.-1) THEN
C
C                 The upper bound is violated.  Test for either a bigger
C                 or smaller ALFA2,  depending on the value of FIRSTV.
C
                  RES = ATX - BU(J)
                  IF ((FIRSTV .AND. PALFA2*ABSATP.GE.RES .OR.
     *                LASTV .AND. PALFA2*ABSATP.LE.RES)
     *                .OR. J.EQ.JSAVE2) THEN
                     IF (APMAX2*ROWNRM*PNORM.LT.ABSATP) THEN
                        APMAX2 = ABSATP/(ROWNRM*PNORM)
                        IF (ABSATP.GE.ONE) THEN
                           ALFA2 = RES/ABSATP
                        ELSE IF (RES.LT.BIGALF*ABSATP) THEN
                           ALFA2 = RES/ABSATP
                        ELSE
                           ALFA2 = BIGALF
                        END IF
                        JADD2 = J
                        ATP2 = ATP
                        HLOW2 = .FALSE.
                     END IF
                  END IF
               END IF
            ELSE IF (ATP.GT.ZERO .AND. JS.NE.-1) THEN
C              ---------------------------------------------------------
C              a'x  is increasing and the upper bound is not violated.
C              ---------------------------------------------------------
C              Test for smaller ALFA1.
C
               IF (BU(J).LT.BIGBND) THEN
                  RES = BU(J) - ATX
                  IF (PALFA1*ATP.GE.RES .OR. J.EQ.JSAVE1) THEN
                     IF (APMAX1*ROWNRM*PNORM.LT.ATP) THEN
                        APMAX1 = ATP/(ROWNRM*PNORM)
                        ALFA1 = RES/ATP
                        JADD1 = J
                        ATP1 = ATP
                        HLOW1 = .FALSE.
                     END IF
                  END IF
               END IF
C
               IF (JS.EQ.-2) THEN
C
C                 The lower bound is violated.  Test for a new ALFA2.
C
                  RES = BL(J) - ATX
                  IF ((FIRSTV .AND. PALFA2*ATP.GE.RES .OR. LASTV .AND.
     *                PALFA2*ATP.LE.RES) .OR. J.EQ.JSAVE2) THEN
                     IF (APMAX2*ROWNRM*PNORM.LT.ATP) THEN
                        APMAX2 = ATP/(ROWNRM*PNORM)
                        IF (ATP.GE.ONE) THEN
                           ALFA2 = RES/ATP
                        ELSE IF (RES.LT.BIGALF*ATP) THEN
                           ALFA2 = RES/ATP
                        ELSE
                           ALFA2 = BIGALF
                        END IF
                        JADD2 = J
                        ATP2 = ATP
                        HLOW2 = .TRUE.
                     END IF
                  END IF
               END IF
            END IF
C
         END IF
   20 CONTINUE
C
C     ==================================================================
C     Determine ALFA, the step to be taken.
C     ==================================================================
C     In the infeasible case, check whether to take the step ALFA2
C     rather than ALFA1...
C
      STEP2 = NUMINF .GT. 0 .AND. JADD2 .GT. 0
C
C     We do so if ALFA2 is less than ALFA1 or (if FIRSTV is false)
C     lies in the range  (ALFA1, PALFA1)  and has a smaller value of
C     ATP.
C
      STEP2 = STEP2 .AND. (ALFA2.LT.ALFA1 .OR. LASTV .AND. ALFA2.LE.
     *        PALFA1 .AND. APMAX2.GE.APMAX1)
C
      IF (STEP2) THEN
         ALFA = ALFA2
         PALFA = PALFA2
         JADD = JADD2
         ATPHIT = ATP2
         HITLOW = HLOW2
      ELSE
         ALFA = ALFA1
         PALFA = PALFA1
         JADD = JADD1
         ATPHIT = ATP1
         HITLOW = HLOW1
C
C        If ALFA1 is negative, the constraint to be added (JADD)
C        remains unchanged, but ALFA may be shortened to the step
C        to the nearest perturbed satisfied constraint along  - P.
C
         NEGSTP = ALFA .LT. ZERO
         IF (NEGSTP) THEN
            CALL E04UCH(FIRSTV,NEGSTP,BIGALF,BIGBND,PNORM,JADD1,JADD2,
     *                  PALFA1,PALFA2,ISTATE,N,NCTOTL,ANORM,AP,AX,BL,BU,
     *                  FEATOL,P,X)
C
            ALFA = -MIN(ABS(ALFA),PALFA1)
         END IF
      END IF
C
C     Test for undefined or infinite step.
C
      IF (JADD.EQ.0) THEN
         ALFA = BIGALF
         PALFA = BIGALF
      END IF
C
      IF (ALFA.GE.BIGALF) INFORM = 3
C
      RETURN
C
C
C     End of  E04UCG. (CMALF)
C
      END
