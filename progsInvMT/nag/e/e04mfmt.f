      SUBROUTINE E04MFM(PRBTYP,MSGLVL,N,LDA,LDT,NACTIV,NFREE,NZ,ISTATE,
     *                  KACTIV,KX,ZEROLM,NOTOPT,NUMINF,TRUSML,SMLLST,
     *                  JSMLST,KSMLST,TINYST,JTINY,JINF,TRUBIG,BIGGST,
     *                  JBIGST,KBIGST,A,ANORMS,GQ,RLAMDA,T,WTINF)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1559 (JUN 1995).
C
C     ******************************************************************
C     E04MFM  first computes the Lagrange multiplier estimates for the
C     given working set.  It then determines the values and indices of
C     certain significant multipliers.  In this process, the multipliers
C     for inequalities at their upper bounds are adjusted so that a
C     negative multiplier for an inequality constraint indicates non-
C     optimality.  All adjusted multipliers are scaled by the 2-norm
C     of the associated constraint row.  In the following, the term
C     minimum refers to the ordering of numbers on the real line,  and
C     not to their magnitude.
C
C     JSMLST          is the index of the constraint whose multiplier is
C                     the minimum of the set of adjusted multipliers
C                     with values less than  small.
C     RLAMDA(KSMLST)  is the associated multiplier.
C
C     JBIGST          is the index of the constraint whose multiplier is
C                     the largest of the set of adjusted multipliers
C                     with values greater than (1 + small).
C     RLAMDA(KBIGST)  is the associated multiplier.
C
C     On exit,  elements  1  thru  NACTIV  of  RLAMDA  contain the
C     unadjusted multipliers for the general constraints.  Elements
C     NACTIV  onwards of  RLAMDA  contain the unadjusted multipliers
C     for the bounds.
C
C     Original version written 31-October-1984.
C     Based on a version of  LSMULS  dated 30-June-1986.
C     This version of  E04MFM  dated 14-Sep-92.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGGST, SMLLST, TINYST, TRUBIG, TRUSML, ZEROLM
      INTEGER           JBIGST, JINF, JSMLST, JTINY, KBIGST, KSMLST,
     *                  LDA, LDT, MSGLVL, N, NACTIV, NFREE, NOTOPT,
     *                  NUMINF, NZ
      CHARACTER*2       PRBTYP
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), ANORMS(*), GQ(N), RLAMDA(N), T(LDT,*),
     *                  WTINF(*)
      INTEGER           ISTATE(*), KACTIV(N), KX(N)
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  ANORMJ, BLAM, RLAM, SCDLAM
      INTEGER           I, IS, J, K, KK, L, NFIXED
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DTRSV, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
C     .. Executable Statements ..
C
      NFIXED = N - NFREE
C
      JTINY = 0
      JSMLST = 0
      KSMLST = 0
C
      JBIGST = 0
      KBIGST = 0
C
C     ------------------------------------------------------------------
C     Compute  JSMLST  for regular constraints and temporary bounds.
C     ------------------------------------------------------------------
C     First, compute the Lagrange multipliers for the general
C     constraints in the working set, by solving  T'*lamda = Y'g.
C
      IF (N.GT.NZ) CALL DCOPY(N-NZ,GQ(NZ+1),1,RLAMDA,1)
      IF (NACTIV.GT.0) CALL DTRSV('U','T','N',NACTIV,T(1,NZ+1),LDT,
     *                            RLAMDA,1)
C
C     -----------------------------------------------------------------
C     Now set elements  NACTIV, NACTIV+1,... of  RLAMDA  equal to
C     the multipliers for the bound constraints.
C     -----------------------------------------------------------------
      DO 40 L = 1, NFIXED
         J = KX(NFREE+L)
         BLAM = RLAMDA(NACTIV+L)
         DO 20 K = 1, NACTIV
            I = KACTIV(K)
            BLAM = BLAM - A(I,J)*RLAMDA(NACTIV-K+1)
   20    CONTINUE
         RLAMDA(NACTIV+L) = BLAM
   40 CONTINUE
C
C     -----------------------------------------------------------------
C     Find  JSMLST  and  KSMLST.
C     -----------------------------------------------------------------
      DO 60 K = 1, N - NZ
         IF (K.GT.NACTIV) THEN
            J = KX(NZ+K)
         ELSE
            J = KACTIV(NACTIV-K+1) + N
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
C
            IF (SCDLAM.LT.ZEROLM) THEN
               IF (NUMINF.EQ.0) NOTOPT = NOTOPT + 1
C
               IF (SCDLAM.LT.SMLLST) THEN
                  SMLLST = SCDLAM
                  TRUSML = RLAMDA(K)
                  JSMLST = J
                  KSMLST = K
               END IF
            ELSE IF (SCDLAM.LT.TINYST) THEN
               TINYST = SCDLAM
               JTINY = J
            END IF
         END IF
C
         SCDLAM = RLAM/WTINF(J)
         IF (SCDLAM.GT.BIGGST .AND. J.GT.JINF) THEN
            BIGGST = SCDLAM
            TRUBIG = RLAMDA(K)
            JBIGST = J
            KBIGST = K
         END IF
   60 CONTINUE
C
C     -----------------------------------------------------------------
C     If required, print the multipliers.
C     -----------------------------------------------------------------
      IF (MSGLVL.GE.20) THEN
         IF (ISUMM.GE.0) THEN
            IF (NFIXED.GT.0) THEN
               WRITE (REC,FMT=99999) PRBTYP
               CALL X04BAY(ISUMM,2,REC)
               DO 80 K = 1, NFIXED, 4
                  WRITE (REC,FMT=99998) (KX(NFREE+KK),RLAMDA(NACTIV+KK),
     *              KK=K,MIN(K+3,NFIXED))
                  CALL X04BAF(ISUMM,REC(1))
   80          CONTINUE
            END IF
            IF (NACTIV.GT.0) THEN
               WRITE (REC,FMT=99997) PRBTYP
               CALL X04BAY(ISUMM,2,REC)
               DO 100 K = 1, NACTIV, 4
                  WRITE (REC,FMT=99998) (KACTIV(KK),RLAMDA(NACTIV-KK+1),
     *              KK=K,MIN(K+3,NACTIV))
                  CALL X04BAF(ISUMM,REC(1))
  100          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of  E04MFM.  (CMMUL1)
C
99999 FORMAT (/' Multipliers for the ',A2,' bound  constraints   ')
99998 FORMAT (4(I5,1P,D11.2))
99997 FORMAT (/' Multipliers for the ',A2,' linear constraints   ')
      END
