      SUBROUTINE E04NCK(PRBTYP,MSGLVL,N,NACTIV,NFREE,LDA,LDT,NUMINF,NZ,
     *                  NRZ,ISTATE,KACTIV,KX,DINKY,JSMLST,KSMLST,JINF,
     *                  JTINY,JBIGST,KBIGST,TRULAM,A,ANORMS,GQ,RLAMDA,T,
     *                  WTINF)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1063 (JUL 1993).
C     MARK 17 REVISED. IER-1575 (JUN 1995).
C
C     ******************************************************************
C     E04NCK  first computes the Lagrange multiplier estimates for the
C     given working set.  It then determines the values and indices of
C     certain significant multipliers.  In this process, the multipliers
C     for inequalities at their upper bounds are adjusted so that a
C     negative multiplier for an inequality constraint indicates non-
C     optimality.  All adjusted multipliers are scaled by the 2-norm
C     of the associated constraint row.  In the following, the term
C     minimum refers to the ordering of numbers on the real line,  and
C     not to their magnitude.
C
C     JSMLST  is the index of the minimum of the set of adjusted
C             multipliers with values less than  - DINKY.  A negative
C             JSMLST defines the index in Q'g of the artificial
C             constraint to be deleted.
C     KSMLST  marks the position of general constraint JSMLST in KACTIV.
C
C     JBIGST  is the index of the largest of the set of adjusted
C             multipliers with values greater than (1 + DINKY).
C     KBIGST  marks its position in KACTIV.
C
C     On exit,  elements 1 thru NACTIV of RLAMDA contain the unadjusted
C     multipliers for the general constraints.  Elements NACTIV onwards
C     of RLAMDA contain the unadjusted multipliers for the bounds.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 31-October-1984.
C     This version of E04NCK dated 14-Sep-92.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DINKY, TRULAM
      INTEGER           JBIGST, JINF, JSMLST, JTINY, KBIGST, KSMLST,
     *                  LDA, LDT, MSGLVL, N, NACTIV, NFREE, NRZ, NUMINF,
     *                  NZ
      CHARACTER*2       PRBTYP
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), ANORMS(*), GQ(N), RLAMDA(N), T(LDT,*),
     *                  WTINF(*)
      INTEGER           ISTATE(*), KACTIV(N), KX(N)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  ANORMJ, BIGGST, BLAM, RLAM, SCDLAM, SMLLST,
     *                  TINYLM
      INTEGER           I, IS, J, K, L, NFIXED
C     .. Local Arrays ..
      CHARACTER*80      REC(80)
C     .. External Subroutines ..
      EXTERNAL          DCOPY, E04NBT, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
C     .. Executable Statements ..
      NFIXED = N - NFREE
C
      JSMLST = 0
      KSMLST = 0
      SMLLST = -DINKY
C
      TINYLM = DINKY
      JTINY = 0
C
      JBIGST = 0
      KBIGST = 0
      BIGGST = ONE + DINKY
C
      IF (NRZ.LT.NZ) THEN
C        ---------------------------------------------------------------
C        Compute JSMLST for the artificial constraints.
C        ---------------------------------------------------------------
         DO 20 J = NRZ + 1, NZ
            RLAM = -ABS(GQ(J))
            IF (RLAM.LT.SMLLST) THEN
               SMLLST = RLAM
               JSMLST = -J
            ELSE IF (RLAM.LT.TINYLM) THEN
               TINYLM = RLAM
               JTINY = J
            END IF
   20    CONTINUE
C
         IF (MSGLVL.GE.20) THEN
            IF (ISUMM.GE.0) THEN
               WRITE (REC,FMT=99999)
               CALL X04BAY(ISUMM,2,REC)
               DO 40 J = NRZ + 1, NZ, 4
                  WRITE (REC,FMT=99996) (GQ(K),K=J,MIN(J+3,NZ))
                  CALL X04BAF(ISUMM,REC(1))
   40          CONTINUE
            END IF
         END IF
C
      END IF
C
C     ---------------------------------------------------------------
C     Compute JSMLST for regular constraints and temporary bounds.
C     ---------------------------------------------------------------
C     First, compute the Lagrange multipliers for the general
C     constraints in the working set, by solving  T'*lamda = Y'g.
C
      IF (N.GT.NZ) CALL DCOPY(N-NZ,GQ(NZ+1),1,RLAMDA,1)
      IF (NACTIV.GT.0) CALL E04NBT(2,LDT,NACTIV,T(1,NZ+1),RLAMDA)
C
C     --------------------------------------------------------------
C     Now set elements NACTIV, NACTIV+1,... of  RLAMDA  equal to
C     the multipliers for the bound constraints.
C     --------------------------------------------------------------
      DO 80 L = 1, NFIXED
         J = KX(NFREE+L)
         BLAM = RLAMDA(NACTIV+L)
         DO 60 K = 1, NACTIV
            I = KACTIV(K)
            BLAM = BLAM - A(I,J)*RLAMDA(K)
   60    CONTINUE
         RLAMDA(NACTIV+L) = BLAM
   80 CONTINUE
C
C     --------------------------------------------------------------
C     Find JSMLST and KSMLST.
C     --------------------------------------------------------------
      DO 100 K = 1, N - NZ
         IF (K.GT.NACTIV) THEN
            J = KX(NZ+K)
         ELSE
            J = KACTIV(K) + N
         END IF
C
         IS = ISTATE(J)
C
         I = J - N
         IF (J.LE.N) ANORMJ = ONE
         IF (J.GT.N) ANORMJ = ANORMS(I)
C
         RLAM = RLAMDA(K)
C
C        Change the sign of the estimate if the constraint is in
C        the working set at its upper bound.
C
         IF (IS.EQ.2) RLAM = -RLAM
         IF (IS.EQ.3) RLAM = ABS(RLAM)
         IF (IS.EQ.4) RLAM = -ABS(RLAM)
C
         IF (IS.NE.3) THEN
            SCDLAM = RLAM*ANORMJ
            IF (SCDLAM.LT.SMLLST) THEN
               SMLLST = SCDLAM
               JSMLST = J
               KSMLST = K
            ELSE IF (SCDLAM.LT.TINYLM) THEN
               TINYLM = SCDLAM
               JTINY = J
            END IF
         END IF
C
         IF (NUMINF.GT.0 .AND. J.GT.JINF) THEN
            SCDLAM = RLAM/WTINF(J)
            IF (SCDLAM.GT.BIGGST) THEN
               BIGGST = SCDLAM
               TRULAM = RLAMDA(K)
               JBIGST = J
               KBIGST = K
            END IF
         END IF
  100 CONTINUE
C
C     --------------------------------------------------------------
C     If required, print the multipliers.
C     --------------------------------------------------------------
      IF (MSGLVL.GE.20) THEN
         IF (ISUMM.GE.0) THEN
            IF (NFIXED.GT.0) THEN
               WRITE (REC,FMT=99998) PRBTYP
               CALL X04BAY(ISUMM,2,REC)
               DO 120 J = 1, NFIXED, 4
                  WRITE (REC,FMT=99995) (KX(NFREE+K),RLAMDA(NACTIV+K),
     *              K=J,MIN(J+3,NFIXED))
                  CALL X04BAF(ISUMM,REC(1))
  120          CONTINUE
            END IF
            IF (NACTIV.GT.0) THEN
               WRITE (REC,FMT=99997) PRBTYP
               CALL X04BAY(ISUMM,2,REC)
               DO 140 J = 1, NACTIV, 4
                  WRITE (REC,FMT=99995) (KACTIV(K),RLAMDA(K),K=J,
     *              MIN(J+3,NACTIV))
                  CALL X04BAF(ISUMM,REC(1))
  140          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C
C     End of  E04NCK. (LSMULS)
C
99999 FORMAT (/' Multipliers for the artificial constraints        ')
99998 FORMAT (/' Multipliers for the ',A2,' bound  constraints   ')
99997 FORMAT (/' Multipliers for the ',A2,' linear constraints   ')
99996 FORMAT (4(5X,1P,D11.2))
99995 FORMAT (4(I5,1P,D11.2))
      END
