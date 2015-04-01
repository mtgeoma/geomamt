      SUBROUTINE F08JJZ(IJOB,NITMAX,N,MMAX,MINP,NBMIN,ABSTOL,RELTOL,
     *                  PIVMIN,D,E,E2,NVAL,AB,C,MOUT,NAB,WORK,IWORK,
     *                  INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             DLAEBZ(IJOB,NITMAX,N,MMAX,MINP,NBMIN,ABSTOL,
C    *                  RELTOL,PIVMIN,D,E,E2,NVAL,AB,C,MOUT,NAB,WORK,
C    *                  IWORK,INFO)
C
C  Purpose
C  =======
C
C  DLAEBZ contains the iteration loops which compute and use the
C  function N(w), which is the count of eigenvalues of a symmetric
C  tridiagonal matrix T less than or equal to its argument  w.  It
C  performs a choice of two types of loops:
C
C  IJOB=1, followed by
C  IJOB=2: It takes as input a list of intervals and returns a list of
C          sufficiently small intervals whose union contains the same
C          eigenvalues as the union of the original intervals.
C          The input intervals are (AB(j,1),AB(j,2)], j=1,...,MINP.
C          The output interval (AB(j,1),AB(j,2)] will contain
C          eigenvalues NAB(j,1)+1,...,NAB(j,2), where 1 <= j <= MOUT.
C
C  IJOB=3: It performs a binary search in each input interval
C          (AB(j,1),AB(j,2)] for a point  w(j)  such that
C          N(w(j))=NVAL(j), and uses  C(j)  as the starting point of
C          the search.  If such a w(j) is found, then on output
C          AB(j,1)=AB(j,2)=w.  If no such w(j) is found, then on output
C          (AB(j,1),AB(j,2)] will be a small interval containing the
C          point where N(w) jumps through NVAL(j), unless that point
C          lies outside the initial interval.
C
C  Note that the intervals are in all cases half-open intervals,
C  i.e., of the form  (a,b] , which includes  b  but not  a .
C
C  See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
C  Matrix", Report CS41, Computer Science Dept., Stanford
C  University, July 21, 1966
C
C  Note: the arguments are, in general, *not* checked for unreasonable
C  values.
C
C  Arguments
C  =========
C
C  IJOB    (input) INTEGER
C          Specifies what is to be done:
C          = 1:  Compute NAB for the initial intervals.
C          = 2:  Perform bisection iteration to find eigenvalues of T.
C          = 3:  Perform bisection iteration to invert N(w), i.e.,
C                to find a point which has a specified number of
C                eigenvalues of T to its left.
C          Other values will cause DLAEBZ to return with INFO=-1.
C
C  NITMAX  (input) INTEGER
C          The maximum number of "levels" of bisection to be
C          performed, i.e., an interval of width W will not be made
C          smaller than 2^(-NITMAX) * W.  If not all intervals
C          have converged after NITMAX iterations, then INFO is set
C          to the number of non-converged intervals.
C
C  N       (input) INTEGER
C          The dimension n of the tridiagonal matrix T.  It must be at
C          least 1.
C
C  MMAX    (input) INTEGER
C          The maximum number of intervals.  If more than MMAX intervals
C          are generated, then DLAEBZ will quit with INFO=MMAX+1.
C
C  MINP    (input) INTEGER
C          The initial number of intervals.  It may not be greater than
C          MMAX.
C
C  NBMIN   (input) INTEGER
C          The smallest number of intervals that should be processed
C          using a vector loop.  If zero, then only the scalar loop
C          will be used.
C
C  ABSTOL  (input) DOUBLE PRECISION
C          The minimum (absolute) width of an interval.  When an
C          interval is narrower than ABSTOL, or than RELTOL times the
C          larger (in magnitude) endpoint, then it is considered to be
C          sufficiently small, i.e., converged.  This must be at least
C          zero.
C
C  RELTOL  (input) DOUBLE PRECISION
C          The minimum relative width of an interval.  When an interval
C          is narrower than ABSTOL, or than RELTOL times the larger (in
C          magnitude) endpoint, then it is considered to be
C          sufficiently small, i.e., converged.  Note: this should
C          always be at least radix*machine epsilon.
C
C  PIVMIN  (input) DOUBLE PRECISION
C          The minimum absolute value of a "pivot" in the Sturm
C          sequence loop.  This *must* be at least  max |e(j)**2| *
C          safe_min  and at least safe_min, where safe_min is at least
C          the smallest number that can divide one without overflow.
C
C  D       (input) DOUBLE PRECISION array, dimension (N)
C          The diagonal entries of the tridiagonal matrix T.  To avoid
C          underflow, the matrix should be scaled so that its largest
C          entry is no greater than  overflow**(1/2) * underflow**(1/4)
C          in absolute value.  To assure the most accurate computation
C          of small eigenvalues, the matrix should be scaled to be
C          not much smaller than that, either.
C
C  E       (input) DOUBLE PRECISION array, dimension (N)
C          The offdiagonal entries of the tridiagonal matrix T in
C          positions 1 through N-1.  E(N) is arbitrary.
C          To avoid underflow, the
C          matrix should be scaled so that its largest entry is no
C          greater than  overflow**(1/2) * underflow**(1/4) in absolute
C          value.  To assure the most accurate computation of small
C          eigenvalues, the matrix should be scaled to be not much
C          smaller than that, either.
C
C  E2      (input) DOUBLE PRECISION array, dimension (N)
C          The squares of the offdiagonal entries of the tridiagonal
C          matrix T.  E2(N) is ignored.
C
C  NVAL    (input/output) INTEGER array, dimension (MINP)
C          If IJOB=1 or 2, not referenced.
C          If IJOB=3, the desired values of N(w).  The elements of NVAL
C          will be reordered to correspond with the intervals in AB.
C          Thus, NVAL(j) on output will not, in general be the same as
C          NVAL(j) on input, but it will correspond with the interval
C          (AB(j,1),AB(j,2)] on output.
C
C  AB      (input/output) DOUBLE PRECISION array, dimension (MMAX,2)
C          The endpoints of the intervals.  AB(j,1) is  a(j), the left
C          endpoint of the j-th interval, and AB(j,2) is b(j), the
C          right endpoint of the j-th interval.  The input intervals
C          will, in general, be modified, split, and reordered by the
C          calculation.
C
C  C       (input/workspace) DOUBLE PRECISION array, dimension (MMAX)
C          If IJOB=1, ignored.
C          If IJOB=2, workspace.
C          If IJOB=3, then on input C(j) should be initialized to the
C          first search point in the binary search.
C
C  MOUT    (output) INTEGER
C          If IJOB=1, the number of eigenvalues in the intervals.
C          If IJOB=2 or 3, the number of intervals output.
C          If IJOB=3, MOUT will equal MINP.
C
C  NAB     (input/output) INTEGER array, dimension (MMAX,2)
C          If IJOB=1, then on output NAB(i,j) will be set to N(AB(i,j)).
C          If IJOB=2, then on input, NAB(i,j) should be set.  It must
C             satisfy the condition:
C             N(AB(i,1)) <= NAB(i,1) <= NAB(i,2) <= N(AB(i,2)),
C             which means that in interval i only eigenvalues
C             NAB(i,1)+1,...,NAB(i,2) will be considered.  Usually,
C             NAB(i,j)=N(AB(i,j)), from a previous call to DLAEBZ with
C             IJOB=1.
C             On output, NAB(i,j) will contain
C             max(na(k),min(nb(k),N(AB(i,j)))), where k is the index of
C             the input interval that the output interval
C             (AB(j,1),AB(j,2)] came from, and na(k) and nb(k) are the
C             the input values of NAB(k,1) and NAB(k,2).
C          If IJOB=3, then on output, NAB(i,j) contains N(AB(i,j)),
C             unless N(w) > NVAL(i) for all search points  w , in which
C             case NAB(i,1) will not be modified, i.e., the output
C             value will be the same as the input value (modulo
C             reorderings -- see NVAL and AB), or unless N(w) < NVAL(i)
C             for all search points  w , in which case NAB(i,2) will
C             not be modified.  Normally, NAB should be set to some
C             distinctive value(s) before DLAEBZ is called.
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension (MMAX)
C          Workspace.
C
C  IWORK   (workspace) INTEGER array, dimension (MMAX)
C          Workspace.
C
C  INFO    (output) INTEGER
C          = 0:       All intervals converged.
C          = 1--MMAX: The last INFO intervals did not converge.
C          = MMAX+1:  More than MMAX intervals were generated.
C
C  Further Details
C  ===============
C
C      This routine is intended to be called only by other LAPACK
C  routines, thus the interface is less user-friendly.  It is intended
C  for two purposes:
C
C  (a) finding eigenvalues.  In this case, DLAEBZ should have one or
C      more initial intervals set up in AB, and DLAEBZ should be called
C      with IJOB=1.  This sets up NAB, and also counts the eigenvalues.
C      Intervals with no eigenvalues would usually be thrown out at
C      this point.  Also, if not all the eigenvalues in an interval i
C      are desired, NAB(i,1) can be increased or NAB(i,2) decreased.
C      For example, set NAB(i,1)=NAB(i,2)-1 to get the largest
C      eigenvalue.  DLAEBZ is then called with IJOB=2 and MMAX
C      no smaller than the value of MOUT returned by the call with
C      IJOB=1.  After this (IJOB=2) call, eigenvalues NAB(i,1)+1
C      through NAB(i,2) are approximately AB(i,1) (or AB(i,2)) to the
C      tolerance specified by ABSTOL and RELTOL.
C
C  (b) finding an interval (a',b'] containing eigenvalues w(f),...,w(l).
C      In this case, start with a Gershgorin interval  (a,b).  Set up
C      AB to contain 2 search intervals, both initially (a,b).  One
C      NVAL entry should contain  f-1  and the other should contain  l
C      , while C should contain a and b, resp.  NAB(i,1) should be -1
C      and NAB(i,2) should be N+1, to flag an error if the desired
C      interval does not lie in (a,b).  DLAEBZ is then called with
C      IJOB=3.  On exit, if w(f-1) < w(f), then one of the intervals --
C      j -- will have AB(j,1)=AB(j,2) and NAB(j,1)=NAB(j,2)=f-1, while
C      if, to the specified tolerance, w(f-k)=...=w(f+r), k > 0 and r
C      >= 0, then the interval will have  N(AB(j,1))=NAB(j,1)=f-k and
C      N(AB(j,2))=NAB(j,2)=f+r.  The cases w(l) < w(l+1) and
C      w(l-r)=...=w(l+k) are handled similarly.
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, TWO, HALF
      PARAMETER         (ZERO=0.0D0,TWO=2.0D0,HALF=1.0D0/TWO)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ABSTOL, PIVMIN, RELTOL
      INTEGER           IJOB, INFO, MINP, MMAX, MOUT, N, NBMIN, NITMAX
C     .. Array Arguments ..
      DOUBLE PRECISION  AB(MMAX,*), C(*), D(*), E(*), E2(*), WORK(*)
      INTEGER           IWORK(*), NAB(MMAX,*), NVAL(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  TMP1, TMP2
      INTEGER           ITMP1, ITMP2, J, JI, JIT, JP, KF, KFNEW, KL,
     *                  KLNEW
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Executable Statements ..
C
C     Check for Errors
C
      INFO = 0
      IF (IJOB.LT.1 .OR. IJOB.GT.3) THEN
         INFO = -1
         RETURN
      END IF
C
C     Initialize NAB
C
      IF (IJOB.EQ.1) THEN
C
C        Compute the number of eigenvalues in the initial intervals.
C
         MOUT = 0
         DO 60 JI = 1, MINP
            DO 40 JP = 1, 2
               TMP1 = D(1) - AB(JI,JP)
               IF (ABS(TMP1).LT.PIVMIN) TMP1 = -PIVMIN
               NAB(JI,JP) = 0
               IF (TMP1.LE.ZERO) NAB(JI,JP) = 1
C
               DO 20 J = 2, N
                  TMP1 = D(J) - E2(J-1)/TMP1 - AB(JI,JP)
                  IF (ABS(TMP1).LT.PIVMIN) TMP1 = -PIVMIN
                  IF (TMP1.LE.ZERO) NAB(JI,JP) = NAB(JI,JP) + 1
   20          CONTINUE
   40       CONTINUE
            MOUT = MOUT + NAB(JI,2) - NAB(JI,1)
   60    CONTINUE
         RETURN
      END IF
C
C     Initialize for loop
C
C     KF and KL have the following meaning:
C        Intervals 1,...,KF-1 have converged.
C        Intervals KF,...,KL  still need to be refined.
C
      KF = 1
      KL = MINP
C
C     If IJOB=2, initialize C.
C     If IJOB=3, use the user-supplied starting point.
C
      IF (IJOB.EQ.2) THEN
         DO 80 JI = 1, MINP
            C(JI) = HALF*(AB(JI,1)+AB(JI,2))
   80    CONTINUE
      END IF
C
C     Iteration loop
C
      DO 260 JIT = 1, NITMAX
C
C        Loop over intervals
C
         IF (KL-KF+1.GE.NBMIN .AND. NBMIN.GT.0) THEN
C
C           Begin of Parallel Version of the loop
C
            DO 120 JI = KF, KL
C
C              Compute N(c), the number of eigenvalues less than c
C
               WORK(JI) = D(1) - C(JI)
               IF (ABS(WORK(JI)).LE.PIVMIN) WORK(JI) = -PIVMIN
               IWORK(JI) = 0
               IF (WORK(JI).LE.ZERO) IWORK(JI) = 1
C
               DO 100 J = 2, N
                  WORK(JI) = D(J) - E2(J-1)/WORK(JI) - C(JI)
                  IF (ABS(WORK(JI)).LE.PIVMIN) WORK(JI) = -PIVMIN
                  IF (WORK(JI).LE.ZERO) IWORK(JI) = IWORK(JI) + 1
  100          CONTINUE
  120       CONTINUE
C
            IF (IJOB.LE.2) THEN
C
C              IJOB=2: Choose all intervals containing eigenvalues.
C
               KLNEW = KL
               DO 140 JI = KF, KL
C
C                 Insure that N(w) is monotone
C
                  IWORK(JI) = MIN(NAB(JI,2),MAX(NAB(JI,1),IWORK(JI)))
C
C                 Update the Queue -- add intervals if both halves
C                 contain eigenvalues.
C
                  IF (IWORK(JI).EQ.NAB(JI,2)) THEN
C
C                    No eigenvalue in the upper interval:
C                    just use the lower interval.
C
                     AB(JI,2) = C(JI)
C
                  ELSE IF (IWORK(JI).EQ.NAB(JI,1)) THEN
C
C                    No eigenvalue in the lower interval:
C                    just use the upper interval.
C
                     AB(JI,1) = C(JI)
                  ELSE
                     KLNEW = KLNEW + 1
                     IF (KLNEW.LE.MMAX) THEN
C
C                       Eigenvalue in both intervals -- add upper to
C                       queue.
C
                        AB(KLNEW,2) = AB(JI,2)
                        NAB(KLNEW,2) = NAB(JI,2)
                        AB(KLNEW,1) = C(JI)
                        NAB(KLNEW,1) = IWORK(JI)
                        AB(JI,2) = C(JI)
                        NAB(JI,2) = IWORK(JI)
                     ELSE
                        INFO = MMAX + 1
                     END IF
                  END IF
  140          CONTINUE
               IF (INFO.NE.0) RETURN
               KL = KLNEW
            ELSE
C
C              IJOB=3: Binary search.  Keep only the interval containing
C                      w   s.t. N(w) = NVAL
C
               DO 160 JI = KF, KL
                  IF (IWORK(JI).LE.NVAL(JI)) THEN
                     AB(JI,1) = C(JI)
                     NAB(JI,1) = IWORK(JI)
                  END IF
                  IF (IWORK(JI).GE.NVAL(JI)) THEN
                     AB(JI,2) = C(JI)
                     NAB(JI,2) = IWORK(JI)
                  END IF
  160          CONTINUE
            END IF
C
         ELSE
C
C           End of Parallel Version of the loop
C
C           Begin of Serial Version of the loop
C
            KLNEW = KL
            DO 200 JI = KF, KL
C
C              Compute N(w), the number of eigenvalues less than w
C
               TMP1 = C(JI)
               TMP2 = D(1) - TMP1
               IF (ABS(TMP2).LE.PIVMIN) TMP2 = -PIVMIN
               ITMP1 = 0
               IF (TMP2.LE.ZERO) ITMP1 = 1
C
               DO 180 J = 2, N
                  TMP2 = D(J) - E2(J-1)/TMP2 - TMP1
                  IF (ABS(TMP2).LE.PIVMIN) TMP2 = -PIVMIN
                  IF (TMP2.LE.ZERO) ITMP1 = ITMP1 + 1
  180          CONTINUE
C
               IF (IJOB.LE.2) THEN
C
C                 IJOB=2: Choose all intervals containing eigenvalues.
C
C                 Insure that N(w) is monotone
C
                  ITMP1 = MIN(NAB(JI,2),MAX(NAB(JI,1),ITMP1))
C
C                 Update the Queue -- add intervals if both halves
C                 contain eigenvalues.
C
                  IF (ITMP1.EQ.NAB(JI,2)) THEN
C
C                    No eigenvalue in the upper interval:
C                    just use the lower interval.
C
                     AB(JI,2) = TMP1
C
                  ELSE IF (ITMP1.EQ.NAB(JI,1)) THEN
C
C                    No eigenvalue in the lower interval:
C                    just use the upper interval.
C
                     AB(JI,1) = TMP1
                  ELSE IF (KLNEW.LT.MMAX) THEN
C
C                    Eigenvalue in both intervals -- add upper to queue.
C
                     KLNEW = KLNEW + 1
                     AB(KLNEW,2) = AB(JI,2)
                     NAB(KLNEW,2) = NAB(JI,2)
                     AB(KLNEW,1) = TMP1
                     NAB(KLNEW,1) = ITMP1
                     AB(JI,2) = TMP1
                     NAB(JI,2) = ITMP1
                  ELSE
                     INFO = MMAX + 1
                     RETURN
                  END IF
               ELSE
C
C                 IJOB=3: Binary search.  Keep only the interval
C                         containing  w  s.t. N(w) = NVAL
C
                  IF (ITMP1.LE.NVAL(JI)) THEN
                     AB(JI,1) = TMP1
                     NAB(JI,1) = ITMP1
                  END IF
                  IF (ITMP1.GE.NVAL(JI)) THEN
                     AB(JI,2) = TMP1
                     NAB(JI,2) = ITMP1
                  END IF
               END IF
  200       CONTINUE
            KL = KLNEW
C
C           End of Serial Version of the loop
C
         END IF
C
C        Check for convergence
C
         KFNEW = KF
         DO 220 JI = KF, KL
            TMP1 = ABS(AB(JI,2)-AB(JI,1))
            TMP2 = MAX(ABS(AB(JI,2)),ABS(AB(JI,1)))
            IF (TMP1.LT.MAX(ABSTOL,PIVMIN,RELTOL*TMP2) .OR. NAB(JI,1)
     *          .GE.NAB(JI,2)) THEN
C
C              Converged -- Swap with position KFNEW,
C                           then increment KFNEW
C
               IF (JI.GT.KFNEW) THEN
                  TMP1 = AB(JI,1)
                  TMP2 = AB(JI,2)
                  ITMP1 = NAB(JI,1)
                  ITMP2 = NAB(JI,2)
                  AB(JI,1) = AB(KFNEW,1)
                  AB(JI,2) = AB(KFNEW,2)
                  NAB(JI,1) = NAB(KFNEW,1)
                  NAB(JI,2) = NAB(KFNEW,2)
                  AB(KFNEW,1) = TMP1
                  AB(KFNEW,2) = TMP2
                  NAB(KFNEW,1) = ITMP1
                  NAB(KFNEW,2) = ITMP2
                  IF (IJOB.EQ.3) THEN
                     ITMP1 = NVAL(JI)
                     NVAL(JI) = NVAL(KFNEW)
                     NVAL(KFNEW) = ITMP1
                  END IF
               END IF
               KFNEW = KFNEW + 1
            END IF
  220    CONTINUE
         KF = KFNEW
C
C        Choose Midpoints
C
         DO 240 JI = KF, KL
            C(JI) = HALF*(AB(JI,1)+AB(JI,2))
  240    CONTINUE
C
C        If no more intervals to refine, quit.
C
         IF (KF.GT.KL) GO TO 280
  260 CONTINUE
C
C     Converged
C
  280 CONTINUE
      INFO = MAX(KL+1-KF,0)
      MOUT = KL
C
      RETURN
C
C     End of F08JJZ (DLAEBZ)
C
      END
