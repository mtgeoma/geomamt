      SUBROUTINE E04MFR(JOB,MSGLVL,N,NCLIN,NMOVED,ITER,NUMINF,ISTATE,BL,
     *                  BU,FEATOL,FEATLU,X)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1562 (JUN 1995).
C
C     ******************************************************************
C     E04MFR performs most of the manoeuvres associated with degeneracy.
C     the degeneracy-resolving strategy operates in the following way.
C
C     Over a cycle of iterations, the feasibility tolerance FEATOL
C     increases slightly (from TOLX0 to TOLX1 in steps of TOLINC).
C     This ensures that all steps taken will be positive.
C
C     After KDEGEN consecutive iterations, variables within
C     FEATOL of their bounds are set exactly on their bounds and x is
C     recomputed to satisfy the general constraints in the working set.
C     FEATOL is then reduced to TOLX0 for the next cycle of iterations.
C
C     FEATLU  is the array of user-supplied feasibility tolerances.
C     FEATOL  is the array of current feasibility tolerances.
C
C     If JOB = 'I', E04MFR initializes the parameters in
C     common block CE04MF:
C
C     TOLX0   is the minimum (scaled) feasibility tolerance.
C     TOLINC  is the scaled increment to the current FEATOL.
C     IDEGEN  is the expand frequency. It is the frequency of resetting
C             FEATOL to (scaled) TOLX0.
C     KDEGEN  is the expand frequency (specified by the user).
C             it is the frequency of resetting FEATOL to (scaled) TOLX0.
C     NDEGEN  counts the number of degenerate steps (incremented
C             by E04MFS).
C     ITNFIX  is the last iteration at which a JOB = 'E' or 'O' entry
C             caused an x to be put on a constraint.
C     NFIX(j) counts the number of times a JOB = 'O' entry has
C             caused the variables to be placed on the working set,
C             where j=1 if infeasible, j=2 if feasible.
C
C     TOLX0*FEATLU and TOLX1*FEATLU are both close to the feasibility
C     Tolerance FEATLU specified by the user.  (They must both be less
C     than FEATLU.)
C
C
C     If JOB = 'E',  E04MFR has been called after a cycle of KDEGEN
C     iterations.  Constraints in the working set are examined to see if
C     any are off their bounds by an amount approaching FEATOL.  NMOVED
C     returns how many.  If NMOVED is positive,  x  is moved onto the
C     constraints in the working set.  It is assumed that the calling
C     routine will then continue iterations.
C
C
C     If JOB = 'O',  E04MFR is being called after a subproblem has been
C     judged optimal, infeasible or unbounded.  Constraint violations
C     are examined as above.
C
C     19-Apr-1988. Original version based on MINOS routine M5DGEN.
C     15-Apr-1994. Expand frequency allowed to expand. This allows
C                  small initial values of KDEGEN.
C     05-Jul-1994. Current version.
C
C     This version of  E04MFR  dated 05-Jul-1994.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, POINT6
      PARAMETER         (ZERO=0.0D+0,POINT6=0.6D+0)
C     .. Scalar Arguments ..
      INTEGER           ITER, MSGLVL, N, NCLIN, NMOVED, NUMINF
      CHARACTER         JOB
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(N+NCLIN), BU(N+NCLIN), FEATLU(N+NCLIN),
     *                  FEATOL(N+NCLIN), X(N)
      INTEGER           ISTATE(N+NCLIN)
C     .. Scalars in Common ..
      DOUBLE PRECISION  TOLINC, TOLX0
      INTEGER           IDEGEN, IPRINT, ISUMM, ITNFIX, KDEGEN, LINES1,
     *                  LINES2, NDEGEN, NOUT
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
      INTEGER           NFIX(2)
C     .. Local Scalars ..
      DOUBLE PRECISION  D, EPSMCH, TOLX1, TOLZ
      INTEGER           IS, J, MAXFIX
      CHARACTER*80      REC
C     .. External Subroutines ..
      EXTERNAL          X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
      COMMON            /CE04MF/TOLX0, TOLINC, IDEGEN, KDEGEN, NDEGEN,
     *                  ITNFIX, NFIX
C     .. Save statement ..
      SAVE              /AX02ZA/, TOLX1, TOLZ
C     .. Executable Statements ..
C
      NMOVED = 0
      IF (JOB.EQ.'i' .OR. JOB.EQ.'I') THEN
C        ---------------------------------------------------------------
C        Job = 'Initialize'.
C        Initialize at the start of each linear problem.
C        KDEGEN  is the expand frequency      and
C        FEATLU  are the user-supplied feasibility tolerances.
C        They are not changed.
C        ---------------------------------------------------------------
         EPSMCH = WMACH(3)
C
         NDEGEN = 0
         ITNFIX = 0
         NFIX(1) = 0
         NFIX(2) = 0
         TOLX0 = 0.5D+0
         TOLX1 = 0.99D+0
         TOLZ = EPSMCH**POINT6
C
         IDEGEN = KDEGEN
         IF (KDEGEN.LT.9999999) THEN
            TOLINC = (TOLX1-TOLX0)/IDEGEN
         ELSE
            TOLINC = ZERO
         END IF
C
         DO 20 J = 1, N + NCLIN
            FEATOL(J) = TOLX0*FEATLU(J)
   20    CONTINUE
      ELSE
C        ---------------------------------------------------------------
C        JOB = 'End of cycle' or 'Optimal'.
C        initialize local variables MAXFIX and TOLZ.
C        ---------------------------------------------------------------
         MAXFIX = 2
C
         IF (JOB.EQ.'o' .OR. JOB.EQ.'O') THEN
C           ------------------------------------------------------------
C           JOB = 'Optimal'.
C           return with NMOVED = 0 if the last call was at the same
C           iteration,  or if there have already been MAXFIX calls with
C           the same state of feasibility.
C           ------------------------------------------------------------
            IF (ITNFIX.EQ.ITER) RETURN
            IF (NUMINF.GT.0) THEN
               J = 1
            ELSE
               J = 2
            END IF
C
            IF (NFIX(J).GE.MAXFIX) RETURN
            NFIX(J) = NFIX(J) + 1
         END IF
C
C        Increase the expand frequency.
C        Reset FEATOL to its minimum value.
C
         IDEGEN = IDEGEN + 10
         IF (KDEGEN.LT.9999999) THEN
            TOLINC = (TOLX1-TOLX0)/IDEGEN
            IDEGEN = IDEGEN + ITER
         ELSE
            TOLINC = ZERO
         END IF
C
         DO 40 J = 1, N + NCLIN
            FEATOL(J) = TOLX0*FEATLU(J)
   40    CONTINUE
C
C        Count the number of times a variable is moved a nontrivial
C        distance onto its bound.
C
         ITNFIX = ITER
C
         DO 60 J = 1, N
            IS = ISTATE(J)
            IF (IS.GT.0 .AND. IS.LT.4) THEN
               IF (IS.EQ.1) THEN
                  D = ABS(X(J)-BL(J))
               ELSE
                  D = ABS(X(J)-BU(J))
               END IF
C
               IF (D.GT.TOLZ) NMOVED = NMOVED + 1
            END IF
   60    CONTINUE
C
         IF (NMOVED.GT.0) THEN
C
C           Some variables were moved onto their bounds.
C
            IF (MSGLVL.GT.0) THEN
               WRITE (REC,FMT=99999) ITER, NMOVED
               CALL X04BAF(IPRINT,REC)
            END IF
         END IF
      END IF
C
C     End of E04MFR.  (CMDGEN)
C
99999 FORMAT (' Itn',I6,' --',I7,'  variables moved to their bounds.')
      END
