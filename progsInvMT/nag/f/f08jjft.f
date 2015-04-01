      SUBROUTINE F08JJF(RANGE,ORDER,N,VL,VU,IL,IU,ABSTOL,D,E,M,NSPLIT,W,
     *                  IBLOCK,ISPLIT,WORK,IWORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             DSTEBZ(RANGE,ORDER,N,VL,VU,IL,IU,ABSTOL,D,E,M,
     *                  NSPLIT,W,IBLOCK,ISPLIT,WORK,IWORK,INFO)
C
C  Purpose
C  =======
C
C  DSTEBZ computes the eigenvalues of a symmetric tridiagonal
C  matrix.  The user may ask for all eigenvalues, all eigenvalues
C  in the interval (VL, VU], or the IL-th through IU-th eigenvalues.
C
C  See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
C  Matrix", Report CS41, Computer Science Dept., Stanford
C  University, July 21, 1966
C
C  Arguments
C  =========
C
C  RANGE   (input) CHARACTER
C          Specifies which eigenvalues are to be found.
C          = 'A': ("All")   all eigenvalues will be found.
C          = 'V': ("Value") all eigenvalues in the half-open interval
C                           (VL, VU] will be found.
C          = 'I': ("Index") the IL-th through IU-th eigenvalues (of the
C                           entire matrix) will be found.
C
C  ORDER   (input) CHARACTER
C          Specifies the order in which the eigenvalues and their
C          block numbers will be stored in W and IBLOCK:
C          = 'B': ("By Block") the eigenvalues will be grouped by
C                              split-off block (see IBLOCK, ISPLIT) and
C                              ordered from smallest to largest within
C                              the block.
C          = 'E': ("Entire matrix")
C                              the eigenvalues for the entire matrix
C                              will be ordered from smallest to
C                              largest.
C
C  N       (input) INTEGER
C          The dimension of the tridiagonal matrix T.
C
C  VL      (input) DOUBLE PRECISION
C          If RANGE='V', the lower bound of the interval to be searched
C          for eigenvalues.  Eigenvalues less than or equal to VL will
C          not be returned.  Not referenced if RANGE='A' or 'I'.
C
C  VU      (input) DOUBLE PRECISION
C          If RANGE='V', the upper bound of the interval to be searched
C          for eigenvalues.  Eigenvalues greater than VU will not be
C          returned.  VU must be greater than VL.  Not referenced if
C          RANGE='A' or 'I'.
C
C  IL      (input) INTEGER
C          If RANGE='I', the index (from smallest to largest) of the
C          smallest eigenvalue to be returned.  IL must be at least 1.
C          Not referenced if RANGE='A' or 'V'.
C
C  IU      (input) INTEGER
C          If RANGE='I', the index (from smallest to largest) of the
C          largest eigenvalue to be returned.  IU must be at least IL
C          and no greater than N.  Not referenced if RANGE='A' or 'V'.
C
C  ABSTOL  (input) DOUBLE PRECISION
C          The absolute tolerance for the eigenvalues.  An eigenvalue
C          (or cluster) is considered to be located if it has been
C          determined to lie in an interval whose width is ABSTOL or
C          less.  If ABSTOL is less than or equal to zero, then ULP*|T|
C          will be used, where |T| means the 1-norm of T.
C
C  D       (input) DOUBLE PRECISION array, dimension (N)
C          The diagonal entries of the tridiagonal matrix T.  To avoid
C          overflow, the matrix must be scaled so that its largest
C          entry is no greater than overflow**(1/2) * underflow**(1/4)
C          in absolute value, and for greatest accuracy, it should not
C          be much smaller than that.
C
C  E       (input) DOUBLE PRECISION array, dimension (N-1)
C          The offdiagonal entries of the tridiagonal matrix T must be
C          in elements 1 through N-1, of the array E.
C          To avoid overflow, the matrix must be scaled so that its
C          largest entry is no greater than overflow**(1/2) *
C          underflow**(1/4) in absolute value, and for greatest
C          accuracy, it should not be much smaller than that.
C
C  M       (output) INTEGER
C          The actual number of eigenvalues found.  It will be in the
C          range 0 to N.  (See also the description of INFO=2,3.)
C
C  NSPLIT  (output) INTEGER
C          The number of diagonal blocks which T is considered to
C          consist of.  It will be between 1 and N.
C
C  W       (output) DOUBLE PRECISION array, dimension (N)
C          On exit, the first M elements of W will contain the
C          eigenvalues.  (DSTEBZ may use the remaining N-M elements as
C          workspace.)
C
C  IBLOCK  (output) INTEGER array, dimension (N)
C          At each row/column j where E(j) is zero or small, the
C          matrix T is considered to split into a block diagonal
C          matrix.  On exit, IBLOCK(i) specifies which block (from 1 to
C          the number of blocks) the eigenvalue W(i) belongs to.
C          (DSTEBZ may use the remaining N-M elements as workspace.)
C          NOTE:  in the (theoretically impossible) event that
C          bisection does not converge for all eigenvalues, INFO is set
C          to 1 or 3, and the ones for which it did not are identified
C          by a *negative* block number.
C
C  ISPLIT  (output) INTEGER array, dimension (N)
C          The splitting points, at which T breaks up into submatrices.
C          The first submatrix consists of rows/columns 1 to ISPLIT(1),
C          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
C          etc., and the NSPLIT-th consists of rows/columns
C          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.
C          (Only the first NSPLIT elements will actually be used, but
C          since the user cannot know a priori what value NSPLIT will
C          have, N words must be reserved for ISPLIT.)
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension (4*N)
C          Workspace.
C
C  IWORK   (workspace) INTEGER array, dimension (3*N)
C          Workspace.
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C          > 0: This case includes a number of "impossible" errors.
C          = 1-3: This includes a number of cases, any or all of
C                 which might conceivably occur.  The output arrays are
C                 defined, but the results may not be as accurate as
C                 desired, or eigenvalues may be missing.
C                 1,3: Bisection failed to converge for some
C                      eigenvalues; these eigenvalues are flagged by a
C                      negative block number.  The effect is that the
C                      eigenvalues may not be as accurate as the
C                      absolute and relative tolerances.  This is
C                      generally caused by arithmetic which is less
C                      accurate than DLAMCH says.
C                 2,3: RANGE='I' only: Not all of the eigenvalues IL:IU
C                      were found.
C                      Effect: M < IU+1-IL
C                      Cause:  non-monotonic arithmetic, causing the
C                              Sturm sequence to be non-monotonic.
C                      Cure:   recalculate, using RANGE='A', and pick
C                              out eigenvalues IL:IU.  In some cases,
C                              increasing the PARAMETER "FUDGE" may
C                              make things work.
C          = 4: RANGE='I', and the Gershgorin interval initially used
C               was too small.  No eigenvalues were computed.
C               Probable cause: your machine has sloppy floating-point
C                               arithmetic.
C               Cure: Increase the PARAMETER "FUDGE", recompile, and
C                     try again.
C
C  Internal Parameters
C  ===================
C
C  RELFAC  DOUBLE PRECISION, default = 2.0e0
C          The relative tolerance.  An interval (a,b] lies within
C          "relative tolerance" if  b-a < RELFAC*ulp*max(|a|,|b|),
C          where "ulp" is the machine precision (distance from 1 to
C          the next larger floating point number.)
C
C  FUDGE   DOUBLE PRECISION, default = 2
C          A "fudge factor" to widen the Gershgorin intervals.  Ideally,
C          a value of 1 should work, but on machines with sloppy
C          arithmetic, this needs to be larger.  The default for
C          publicly released versions should be large enough to handle
C          the worst machine around.  Note that this has no effect
C          on accuracy of the solution.
C
C-----------------------------------------------------------------------
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO, HALF
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,HALF=1.0D0/TWO)
      DOUBLE PRECISION  FUDGE, RELFAC
      PARAMETER         (FUDGE=2.0D0,RELFAC=2.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ABSTOL, VL, VU
      INTEGER           IL, INFO, IU, M, N, NSPLIT
      CHARACTER         ORDER, RANGE
C     .. Array Arguments ..
      DOUBLE PRECISION  D(*), E(*), W(*), WORK(*)
      INTEGER           IBLOCK(*), ISPLIT(*), IWORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ATOLI, BNORM, GL, GU, PIVMIN, RTOLI, SAFEMN,
     *                  TMP1, TMP2, TNORM, ULP, WKILL, WL, WLU, WU, WUL
      INTEGER           IB, IBEGIN, IDISCL, IDISCU, IE, IEND, IINFO, IM,
     *                  IN, IOFF, IORDER, IOUT, IRANGE, ITMAX, ITMP1,
     *                  IW, IWOFF, J, JB, JDISC, JE, NB, NWL, NWU
      LOGICAL           NCNVRG, TOOFEW
C     .. Local Arrays ..
      INTEGER           IDUMMA(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      INTEGER           X02BHF
      EXTERNAL          X02AJF, X02AMF, X02BHF
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07ZAZ, F08JJZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, INT, LOG, MAX, MIN, SQRT
C     .. Executable Statements ..
C
      INFO = 0
C
C     Decode RANGE
C
      IF ((RANGE.EQ.'A' .OR. RANGE.EQ.'a')) THEN
         IRANGE = 1
      ELSE IF ((RANGE.EQ.'V' .OR. RANGE.EQ.'v')) THEN
         IRANGE = 2
      ELSE IF ((RANGE.EQ.'I' .OR. RANGE.EQ.'i')) THEN
         IRANGE = 3
      ELSE
         IRANGE = 0
      END IF
C
C     Decode ORDER
C
      IF ((ORDER.EQ.'B' .OR. ORDER.EQ.'b')) THEN
         IORDER = 2
      ELSE IF ((ORDER.EQ.'E' .OR. ORDER.EQ.'e')) THEN
         IORDER = 1
      ELSE
         IORDER = 0
      END IF
C
C     Check for Errors
C
      IF (IRANGE.LE.0) THEN
         INFO = -1
      ELSE IF (IORDER.LE.0) THEN
         INFO = -2
      ELSE IF (N.LT.0) THEN
         INFO = -3
      ELSE IF (IRANGE.EQ.2 .AND. VL.GE.VU) THEN
         INFO = -5
      ELSE IF (IRANGE.EQ.3 .AND. (IL.LT.1 .OR. IL.GT.MAX(1,N))) THEN
         INFO = -6
      ELSE IF (IRANGE.EQ.3 .AND. (IU.LT.MIN(N,IL) .OR. IU.GT.N)) THEN
         INFO = -7
      END IF
C
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08JJF/DSTEBZ',-INFO)
         RETURN
      END IF
C
C     Initialize error flags
C
      INFO = 0
      NCNVRG = .FALSE.
      TOOFEW = .FALSE.
C
C     Quick return if possible
C
      M = 0
      IF (N.EQ.0) RETURN
C
C     Simplifications:
C
      IF (IRANGE.EQ.3 .AND. IL.EQ.1 .AND. IU.EQ.N) IRANGE = 1
C
C     Get machine constants
C     NB is the minimum vector length for vector bisection, or 0
C     if only scalar is to be done.
C
      SAFEMN = X02AMF()
      ULP = X02AJF()*X02BHF()
      RTOLI = ULP*RELFAC
      CALL F07ZAZ(1,'F08JJF',NB,0)
      IF (NB.LE.1) NB = 0
C
C     Special Case when N=1
C
      IF (N.EQ.1) THEN
         NSPLIT = 1
         ISPLIT(1) = 1
         IF (IRANGE.EQ.2 .AND. (VL.GE.D(1) .OR. VU.LT.D(1))) THEN
            M = 0
         ELSE
            W(1) = D(1)
            IBLOCK(1) = 1
            M = 1
         END IF
         RETURN
      END IF
C
C     Compute Splitting Points
C
      NSPLIT = 1
      WORK(N) = ZERO
      PIVMIN = ONE
C
      DO 20 J = 2, N
         TMP1 = E(J-1)**2
         IF (ABS(D(J)*D(J-1))*ULP**2+SAFEMN.GT.TMP1) THEN
            ISPLIT(NSPLIT) = J - 1
            NSPLIT = NSPLIT + 1
            WORK(J-1) = ZERO
         ELSE
            WORK(J-1) = TMP1
            PIVMIN = MAX(PIVMIN,TMP1)
         END IF
   20 CONTINUE
      ISPLIT(NSPLIT) = N
      PIVMIN = PIVMIN*SAFEMN
C
C     Compute Interval and ATOLI
C
      IF (IRANGE.EQ.3) THEN
C
C        RANGE='I': Compute the interval containing eigenvalues
C                   IL through IU.
C
C        Compute Gershgorin interval for entire (split) matrix
C        and use it as the initial interval
C
         GU = D(1)
         GL = D(1)
         TMP1 = ZERO
C
         DO 40 J = 1, N - 1
            TMP2 = SQRT(WORK(J))
            GU = MAX(GU,D(J)+TMP1+TMP2)
            GL = MIN(GL,D(J)-TMP1-TMP2)
            TMP1 = TMP2
   40    CONTINUE
C
         GU = MAX(GU,D(N)+TMP1)
         GL = MIN(GL,D(N)-TMP1)
         TNORM = MAX(ABS(GL),ABS(GU))
         GL = GL - FUDGE*TNORM*ULP*N - FUDGE*TWO*PIVMIN
         GU = GU + FUDGE*TNORM*ULP*N + FUDGE*PIVMIN
C
C        Compute Iteration parameters
C
         ITMAX = INT((LOG(TNORM+PIVMIN)-LOG(PIVMIN))/LOG(TWO)) + 2
         IF (ABSTOL.LE.ZERO) THEN
            ATOLI = ULP*TNORM
         ELSE
            ATOLI = ABSTOL
         END IF
C
         WORK(N+1) = GL
         WORK(N+2) = GL
         WORK(N+3) = GU
         WORK(N+4) = GU
         WORK(N+5) = GL
         WORK(N+6) = GU
         IWORK(1) = -1
         IWORK(2) = -1
         IWORK(3) = N + 1
         IWORK(4) = N + 1
         IWORK(5) = IL - 1
         IWORK(6) = IU
C
         CALL F08JJZ(3,ITMAX,N,2,2,NB,ATOLI,RTOLI,PIVMIN,D,E,WORK,
     *               IWORK(5),WORK(N+1),WORK(N+5),IOUT,IWORK,W,IBLOCK,
     *               IINFO)
C
         IF (IWORK(6).EQ.IU) THEN
            WL = WORK(N+1)
            WLU = WORK(N+3)
            NWL = IWORK(1)
            WU = WORK(N+4)
            WUL = WORK(N+2)
            NWU = IWORK(4)
         ELSE
            WL = WORK(N+2)
            WLU = WORK(N+4)
            NWL = IWORK(2)
            WU = WORK(N+3)
            WUL = WORK(N+1)
            NWU = IWORK(3)
         END IF
C
         IF (NWL.LT.0 .OR. NWL.GE.N .OR. NWU.LT.1 .OR. NWU.GT.N) THEN
            INFO = 4
            RETURN
         END IF
      ELSE
C
C        RANGE='A' or 'V' -- Set ATOLI
C
         TNORM = MAX(ABS(D(1))+ABS(E(1)),ABS(D(N))+ABS(E(N-1)))
C
         DO 60 J = 2, N - 1
            TNORM = MAX(TNORM,ABS(D(J))+ABS(E(J-1))+ABS(E(J)))
   60    CONTINUE
C
         IF (ABSTOL.LE.ZERO) THEN
            ATOLI = ULP*TNORM
         ELSE
            ATOLI = ABSTOL
         END IF
C
         IF (IRANGE.EQ.2) THEN
            WL = VL
            WU = VU
         END IF
      END IF
C
C     Find Eigenvalues -- Loop Over Blocks and recompute NWL and NWU.
C     NWL accumulates the number of eigenvalues .le. WL,
C     NWU accumulates the number of eigenvalues .le. WU
C
      M = 0
      IEND = 0
      INFO = 0
      NWL = 0
      NWU = 0
C
      DO 140 JB = 1, NSPLIT
         IOFF = IEND
         IBEGIN = IOFF + 1
         IEND = ISPLIT(JB)
         IN = IEND - IOFF
C
         IF (IN.EQ.1) THEN
C
C           Special Case -- IN=1
C
            IF (IRANGE.EQ.1 .OR. WL.GE.D(IBEGIN)-PIVMIN) NWL = NWL + 1
            IF (IRANGE.EQ.1 .OR. WU.GE.D(IBEGIN)-PIVMIN) NWU = NWU + 1
            IF (IRANGE.EQ.1 .OR. (WL.LT.D(IBEGIN)
     *          -PIVMIN .AND. WU.GE.D(IBEGIN)-PIVMIN)) THEN
               M = M + 1
               W(M) = D(IBEGIN)
               IBLOCK(M) = JB
            END IF
         ELSE
C
C           General Case -- IN > 1
C
C           Compute Gershgorin Interval
C           and use it as the initial interval
C
            GU = D(IBEGIN)
            GL = D(IBEGIN)
            TMP1 = ZERO
C
            DO 80 J = IBEGIN, IEND - 1
               TMP2 = ABS(E(J))
               GU = MAX(GU,D(J)+TMP1+TMP2)
               GL = MIN(GL,D(J)-TMP1-TMP2)
               TMP1 = TMP2
   80       CONTINUE
C
            GU = MAX(GU,D(IEND)+TMP1)
            GL = MIN(GL,D(IEND)-TMP1)
            BNORM = MAX(ABS(GL),ABS(GU))
            GL = GL - FUDGE*BNORM*ULP*IN - FUDGE*PIVMIN
            GU = GU + FUDGE*BNORM*ULP*IN + FUDGE*PIVMIN
C
            IF (IRANGE.GT.1) THEN
               IF (GU.LT.WL) THEN
                  NWL = NWL + IN
                  NWU = NWU + IN
                  GO TO 140
               END IF
               GL = MAX(GL,WL)
               GU = MIN(GU,WU)
               IF (GL.GE.GU) GO TO 140
            END IF
C
C           Set Up Initial Interval
C
            WORK(N+1) = GL
            WORK(N+IN+1) = GU
            CALL F08JJZ(1,0,IN,IN,1,NB,ATOLI,RTOLI,PIVMIN,D(IBEGIN),
     *                  E(IBEGIN),WORK(IBEGIN),IDUMMA,WORK(N+1),
     *                  WORK(N+2*IN+1),IM,IWORK,W(M+1),IBLOCK(M+1),
     *                  IINFO)
C
            NWL = NWL + IWORK(1)
            NWU = NWU + IWORK(IN+1)
            IWOFF = M - IWORK(1)
C
C           Compute Eigenvalues
C
            ITMAX = INT((LOG(GU-GL+PIVMIN)-LOG(PIVMIN))/LOG(TWO)) + 2
            CALL F08JJZ(2,ITMAX,IN,IN,1,NB,ATOLI,RTOLI,PIVMIN,D(IBEGIN),
     *                  E(IBEGIN),WORK(IBEGIN),IDUMMA,WORK(N+1),
     *                  WORK(N+2*IN+1),IOUT,IWORK,W(M+1),IBLOCK(M+1),
     *                  IINFO)
C
C           Copy Eigenvalues Into W and IBLOCK
C           Use -JB for block number for unconverged eigenvalues.
C
            DO 120 J = 1, IOUT
               TMP1 = HALF*(WORK(J+N)+WORK(J+IN+N))
C
C              Flag non-convergence.
C
               IF (J.GT.IOUT-IINFO) THEN
                  NCNVRG = .TRUE.
                  IB = -JB
               ELSE
                  IB = JB
               END IF
               DO 100 JE = IWORK(J) + 1 + IWOFF, IWORK(J+IN) + IWOFF
                  W(JE) = TMP1
                  IBLOCK(JE) = IB
  100          CONTINUE
  120       CONTINUE
C
            M = M + IM
         END IF
  140 CONTINUE
C
C     If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU
C     If NWL+1 < IL or NWU > IU, discard extra eigenvalues.
C
      IF (IRANGE.EQ.3) THEN
         IM = 0
         IDISCL = IL - 1 - NWL
         IDISCU = NWU - IU
C
         IF (IDISCL.GT.0 .OR. IDISCU.GT.0) THEN
            DO 160 JE = 1, M
               IF (W(JE).LE.WLU .AND. IDISCL.GT.0) THEN
                  IDISCL = IDISCL - 1
               ELSE IF (W(JE).GE.WUL .AND. IDISCU.GT.0) THEN
                  IDISCU = IDISCU - 1
               ELSE
                  IM = IM + 1
                  W(IM) = W(JE)
                  IBLOCK(IM) = IBLOCK(JE)
               END IF
  160       CONTINUE
            M = IM
         END IF
         IF (IDISCL.GT.0 .OR. IDISCU.GT.0) THEN
C
C           Code to deal with effects of bad arithmetic:
C           Some low eigenvalues to be discarded are not in (WL,WLU],
C           or high eigenvalues to be discarded are not in (WUL,WU]
C           so just kill off the smallest IDISCL/largest IDISCU
C           eigenvalues, by simply finding the smallest/largest
C           eigenvalue(s).
C
C           (If N(w) is monotone non-decreasing, this should never
C               happen.)
C
            IF (IDISCL.GT.0) THEN
               WKILL = WU
               DO 200 JDISC = 1, IDISCL
                  IW = 0
                  DO 180 JE = 1, M
                     IF (IBLOCK(JE).NE.0 .AND. (W(JE)
     *                   .LT.WKILL .OR. IW.EQ.0)) THEN
                        IW = JE
                        WKILL = W(JE)
                     END IF
  180             CONTINUE
                  IBLOCK(IW) = 0
  200          CONTINUE
            END IF
            IF (IDISCU.GT.0) THEN
C
               WKILL = WL
               DO 240 JDISC = 1, IDISCU
                  IW = 0
                  DO 220 JE = 1, M
                     IF (IBLOCK(JE).NE.0 .AND. (W(JE)
     *                   .GT.WKILL .OR. IW.EQ.0)) THEN
                        IW = JE
                        WKILL = W(JE)
                     END IF
  220             CONTINUE
                  IBLOCK(IW) = 0
  240          CONTINUE
            END IF
            IM = 0
            DO 260 JE = 1, M
               IF (IBLOCK(JE).NE.0) THEN
                  IM = IM + 1
                  W(IM) = W(JE)
                  IBLOCK(IM) = IBLOCK(JE)
               END IF
  260       CONTINUE
            M = IM
         END IF
         IF (IDISCL.LT.0 .OR. IDISCU.LT.0) THEN
            TOOFEW = .TRUE.
         END IF
      END IF
C
C     If ORDER='B', do nothing -- the eigenvalues are already sorted
C        by block.
C     If ORDER='E' or 'A', sort the eigenvalues from smallest to largest
C
      IF (IORDER.EQ.1 .AND. NSPLIT.GT.1) THEN
         DO 300 JE = 1, M - 1
            IE = 0
            TMP1 = W(JE)
            DO 280 J = JE + 1, M
               IF (W(J).LT.TMP1) THEN
                  IE = J
                  TMP1 = W(J)
               END IF
  280       CONTINUE
C
            IF (IE.NE.0) THEN
               ITMP1 = IBLOCK(IE)
               W(IE) = W(JE)
               IBLOCK(IE) = IBLOCK(JE)
               W(JE) = TMP1
               IBLOCK(JE) = ITMP1
            END IF
  300    CONTINUE
      END IF
C
      INFO = 0
      IF (NCNVRG) INFO = INFO + 1
      IF (TOOFEW) INFO = INFO + 2
      RETURN
C
C     End of F08JJF (DSTEBZ)
C
      END
