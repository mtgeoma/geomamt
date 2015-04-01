      SUBROUTINE E04MFZ(PRBTYP,MSG,CSET,RSET,UNITQ,ITER,ITMAX,JINF,
     *                  NVIOL,N,NCLIN,LDA,NACTIV,NFREE,NRZ,NZ,ISTATE,
     *                  KACTIV,KX,E04MFU,OBJ,NUMINF,XNORM,A,AX,BL,BU,
     *                  CVEC,FEATOL,FEATLU,X,W)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1569 (JUN 1995).
C
C     ******************************************************************
C     E04MFZ  is a subroutine for linear programming.
C     On entry, it is assumed that an initial working set of
C     linear constraints and bounds is available.  The arrays  ISTATE,
C     KACTIV  and  KX  will have been set accordingly
C     and the arrays  T  and  Q  will contain the TQ factorization of
C     the matrix whose rows are the gradients of the active linear
C     constraints with the columns corresponding to the active bounds
C     removed.  The TQ factorization of the resulting (NACTIV by NFREE)
C     matrix is  A(free)*Q = (0 T),  where Q is (NFREE by NFREE) and T
C     is upper-triangular.
C
C     KACTIV holds the general constraint indices in the order in which
C     they were added. The reverse ordering is used for T since new rows
C     are added at the front of T.
C
C     Over a cycle of iterations, the feasibility tolerance FEATOL
C     increases slightly (from TOLX0 to TOLX1 in steps of TOLINC).
C     this ensures that all steps taken will be positive.
C
C     After IDEGEN consecutive iterations, variables within FEATOL of
C     their bounds are set exactly on their bounds and iterative
C     refinement is used to satisfy the constraints in the working set.
C     FEATOL is then reduced to TOLX0 for the next cycle of iterations.
C
C
C     Values of ISTATE(j) for the linear constraints.......
C
C     ISTATE(j)
C     ---------
C          0    constraint j is not in the working set.
C          1    constraint j is in the working set at its lower bound.
C          2    constraint j is in the working set at its upper bound.
C          3    constraint j is in the working set as an equality.
C
C     Constraint j may be violated by as much as FEATOL(j).
C
C     This version of  E04MFZ  dated  10-Apr-94.
C
C     Copyright  1988/1994  Optimates.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LENLC
      PARAMETER         (LENLC=20)
      INTEGER           MXPARM
      PARAMETER         (MXPARM=30)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      CHARACTER*6       EMPTY
      PARAMETER         (EMPTY='      ')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  OBJ, XNORM
      INTEGER           ITER, ITMAX, JINF, LDA, N, NACTIV, NCLIN, NFREE,
     *                  NRZ, NUMINF, NVIOL, NZ
      LOGICAL           CSET, RSET, UNITQ
      CHARACTER*2       PRBTYP
      CHARACTER*6       MSG
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), AX(*), BL(N+NCLIN), BU(N+NCLIN),
     *                  CVEC(*), FEATLU(N+NCLIN), FEATOL(N+NCLIN), W(*),
     *                  X(N)
      INTEGER           ISTATE(N+NCLIN), KACTIV(N), KX(N)
C     .. Subroutine Arguments ..
      EXTERNAL          E04MFU
C     .. Scalars in Common ..
      DOUBLE PRECISION  ALFA, ASIZE, BIGBND, BIGDX, BNDLOW, BNDUPP,
     *                  DTMAX, DTMIN, EPSPT3, EPSPT5, EPSPT8, EPSPT9,
     *                  TOLACT, TOLFEA, TOLINC, TOLRNK, TOLX0, TRULAM
      INTEGER           IDEGEN, IPRINT, IPRNT, ISDEL, ISUMM, ISUMRY,
     *                  ITMAX1, ITMAX2, ITNFIX, JADD, JDEL, KCHK,
     *                  KCYCLE, KDEGEN, LCRASH, LDQ, LDT, LENNAM,
     *                  LINES1, LINES2, LPROB, MAXACT, MAXNZ, MINSUM,
     *                  MM, MSGLC, MXFREE, NCOLT, NDEGEN, NN, NNCLIN,
     *                  NOUT, NPROB
      LOGICAL           HEADER, PRNT
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADLC(23), RPSVLC(MXPARM), WMACH(15)
      INTEGER           IPADLC(15), IPSVLC(MXPARM), LOCLC(LENLC),
     *                  NFIX(2)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALFAP, ALFHIT, BIGALF, BIGGST, CONDMX, CONDRZ,
     *                  CONDT, DINKY, DNORM, DZZ, ERRMAX, FLMAX, GFNORM,
     *                  GRZNRM, GZNORM, OBJSIZ, SMLLST, SUMINF, TINYST,
     *                  TRUBIG, TRUSML, WSSIZE, ZEROLM
      INTEGER           IADD, IFIX, INFORM, IS, IT, J, JBIGST, JMAX,
     *                  JSMLST, JTINY, KBIGST, KDEL, KSMLST, LAD,
     *                  LANORM, LCQ, LD, LDR, LGQ, LQ, LR, LRLAM, LT,
     *                  LWRK, LWTINF, MSGLVL, NCTOTL, NFIXED, NGQ,
     *                  NMOVED, NOTOPT, NTFIXD
      LOGICAL           FIRSTV, FP, HITLOW, LP, MOVE, ONBND, OVERFL,
     *                  UNBNDD
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMLC(MXPARM)
      INTEGER           IPRMLC(MXPARM)
      CHARACTER*80      REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2, F06BLF
      EXTERNAL          DDOT, DNRM2, F06BLF
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DSCAL, E04MFH, E04MFL,
     *                  E04MFM, E04MFQ, E04MFR, E04MFS, E04MFY, E04NBW,
     *                  E04NFP, E04NFR, F06FBF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MOD
C     .. Common blocks ..
      COMMON            /AE04MF/LOCLC
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
      COMMON            /BE04NB/LENNAM, LDT, NCOLT, LDQ
      COMMON            /CE04MF/TOLX0, TOLINC, IDEGEN, KDEGEN, NDEGEN,
     *                  ITNFIX, NFIX
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /DE04MF/ALFA, TRULAM, ISDEL, JDEL, JADD, HEADER,
     *                  PRNT
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
      COMMON            /FE04MF/IPSVLC, IPRNT, ISUMRY, ITMAX1, ITMAX2,
     *                  KCHK, KCYCLE, LCRASH, LPROB, MAXACT, MXFREE,
     *                  MAXNZ, MM, MINSUM, MSGLC, NN, NNCLIN, NPROB,
     *                  IPADLC
      COMMON            /GE04MF/RPSVLC, BIGBND, BIGDX, BNDLOW, BNDUPP,
     *                  TOLACT, TOLFEA, TOLRNK, RPADLC
C     .. Equivalences ..
      EQUIVALENCE       (IPRMLC(1),IPRNT), (RPRMLC(1),BIGBND)
      EQUIVALENCE       (MSGLC,MSGLVL)
C     .. Save statement ..
      SAVE              /AX02ZA/, /AE04MF/, /FE04MF/, /GE04MF/, FIRSTV
C     .. Executable Statements ..
C
C     Specify the machine-dependent parameters.
C
      FLMAX = WMACH(7)
C
      IF (CSET) THEN
         NGQ = 2
      ELSE
         NGQ = 1
      END IF
C
      LP = PRBTYP .EQ. 'lp' .OR. PRBTYP .EQ. 'LP'
      FP = .NOT. LP
C
      LDR = LDT
      IT = 1
C
      LANORM = LOCLC(4)
      LAD = LOCLC(5)
C
      LD = LOCLC(7)
      LGQ = LOCLC(8)
      LCQ = LOCLC(9)
      LRLAM = LOCLC(10)
C
      LR = LOCLC(11)
      LT = LOCLC(12)
      LQ = LOCLC(13)
      LWTINF = LOCLC(14)
      LWRK = LOCLC(15)
C
      IF (ITER.EQ.0) THEN
C        -------------------------
C        First entry.  Initialize.
C        -------------------------
         JADD = 0
         JDEL = 0
         ISDEL = 0
         FIRSTV = .FALSE.
C
         ALFA = ZERO
         DZZ = ONE
      END IF
C
      NCTOTL = N + NCLIN
      NVIOL = 0
C
      CONDMX = FLMAX
C
      CALL E04MFH(N,NCLIN,LDA,ISTATE,BIGBND,NUMINF,SUMINF,BL,BU,A,
     *            FEATOL,W(LGQ),X,W(LWTINF))
C
      IF (NUMINF.GT.0) THEN
         CALL E04NBW(6,N,NZ,NFREE,LDQ,UNITQ,KX,W(LGQ),W(LQ),W(LWRK))
      ELSE IF (LP) THEN
         CALL DCOPY(N,W(LCQ),1,W(LGQ),1)
      END IF
C
      IF (NUMINF.EQ.0 .AND. LP) THEN
         OBJ = DDOT(N,CVEC,1,X,1)
      ELSE
         OBJ = SUMINF
      END IF
C
      MSG = EMPTY
C
C*    ======================Start of main loop==========================
C     +    do while (msg .eq. empty)
   20 IF (MSG.EQ.EMPTY) THEN
C
         GZNORM = ZERO
         IF (NZ.GT.0) GZNORM = DNRM2(NZ,W(LGQ),1)
C
         IF (NRZ.EQ.NZ) THEN
            GRZNRM = GZNORM
         ELSE
            GRZNRM = ZERO
            IF (NRZ.GT.0) GRZNRM = DNRM2(NRZ,W(LGQ),1)
         END IF
C
         GFNORM = GZNORM
         IF (NFREE.GT.0 .AND. NACTIV.GT.0) GFNORM = DNRM2(NFREE,W(LGQ),
     *       1)
C
C        ---------------------------------------------------------------
C        Print the details of this iteration.
C        ---------------------------------------------------------------
C        Define small quantities that reflect the size of x, R and
C        the constraints in the working set.
C
         IF (PRNT) THEN
            CONDT = ONE
            IF (NACTIV.GT.0) CONDT = F06BLF(DTMAX,DTMIN,OVERFL)
C
            CALL E04MFU(PRBTYP,HEADER,RSET,MSGLVL,ITER,ISDEL,JDEL,JADD,
     *                  N,NCLIN,NACTIV,NFREE,NZ,NRZ,LDR,LDT,ISTATE,ALFA,
     *                  CONDRZ,CONDT,DZZ,GZNORM,NUMINF,SUMINF,NOTOPT,
     *                  OBJ,TRULAM,AX,W(LR),W(LT),X,W(LWRK))
            JDEL = 0
            JADD = 0
            ALFA = ZERO
         END IF
C
         IF (NUMINF.GT.0) THEN
            DINKY = EPSPT8*ABS(SUMINF)
         ELSE
            OBJSIZ = ONE + ABS(OBJ)
            WSSIZE = ZERO
            IF (NACTIV.GT.0) WSSIZE = DTMAX
            DINKY = EPSPT8*MAX(WSSIZE,OBJSIZ,GFNORM)
         END IF
C
C        If the reduced gradient Z'g is small enough,
C        Lagrange multipliers will be computed.
C
         IF (NUMINF.EQ.0 .AND. FP) THEN
            MSG = 'feasbl'
            NFIXED = N - NFREE
            CALL F06FBF(NACTIV+NFIXED,ZERO,W(LRLAM),1)
            GO TO 20
         END IF
C
         IF (GRZNRM.LE.DINKY) THEN
C           =========================================================
C           The point  x  is a constrained stationary point.
C           Compute Lagrange multipliers.
C           =========================================================
C           Define what we mean by 'tiny' and non-optimal multipliers.
C
            NOTOPT = 0
            JDEL = 0
            ZEROLM = -DINKY
            SMLLST = -DINKY
            BIGGST = DINKY + ONE
            TINYST = DINKY
C
            CALL E04MFM(PRBTYP,MSGLVL,N,LDA,LDT,NACTIV,NFREE,NZ,ISTATE,
     *                  KACTIV,KX,ZEROLM,NOTOPT,NUMINF,TRUSML,SMLLST,
     *                  JSMLST,KSMLST,TINYST,JTINY,JINF,TRUBIG,BIGGST,
     *                  JBIGST,KBIGST,A,W(LANORM),W(LGQ),W(LRLAM),W(LT),
     *                  W(LWTINF))
C
            IF (NRZ.LT.NZ) THEN
               CALL E04MFL(MSGLVL,N,NRZ,NZ,ZEROLM,NOTOPT,NUMINF,TRUSML,
     *                     SMLLST,JSMLST,TINYST,JTINY,W(LGQ))
            END IF
C
            IF (ABS(JSMLST).GT.0) THEN
C              ------------------------------------------------------
C              Delete a constraint.
C              ------------------------------------------------------
C              E04MFM  or  E04MFL  found a non-optimal multiplier.
C
               TRULAM = TRUSML
               JDEL = JSMLST
C
               IF (JSMLST.GT.0) THEN
C
C                 Regular constraint.
C
                  KDEL = KSMLST
                  ISDEL = ISTATE(JDEL)
                  ISTATE(JDEL) = 0
               END IF
            ELSE IF (MINSUM.GT.0) THEN
               IF (NUMINF.GT.0 .AND. JBIGST.GT.0) THEN
C
C                 No feasible point exists for the constraints but
C                 the sum of the constraint violations can be reduced
C                 by moving off constraints with multipliers greater
C                 than 1.
C
                  JDEL = JBIGST
                  KDEL = KBIGST
                  ISDEL = ISTATE(JDEL)
                  IF (TRUBIG.LE.ZERO) IS = -1
                  IF (TRUBIG.GT.ZERO) IS = -2
                  ISTATE(JDEL) = IS
                  TRULAM = TRUBIG
                  FIRSTV = .TRUE.
                  NUMINF = NUMINF + 1
               END IF
            END IF
C
            IF (JDEL.EQ.0) THEN
               IF (NUMINF.GT.0) THEN
                  MSG = 'infeas'
               ELSE
                  MSG = 'optiml'
               END IF
               GO TO 20
            END IF
C
C           Constraint  JDEL  has been deleted.
C           Update the  TQ  factorization.
C
            CALL E04NFP(UNITQ,IT,N,NACTIV,NFREE,NGQ,NZ,NRZ,LDA,LDQ,LDT,
     *                  JDEL,KDEL,KACTIV,KX,A,W(LT),W(LGQ),W(LQ),W(LD),
     *                  W(LRLAM))
            IF (RSET) CALL E04MFY(NRZ,LDR,W(LR),ONE)
C
            PRNT = .FALSE.
         ELSE
C           ============================================================
C           Compute a search direction.
C           ============================================================
            IF (ITER.GE.ITMAX) THEN
               MSG = 'itnlim'
               GO TO 20
            END IF
C
            PRNT = .TRUE.
            ITER = ITER + 1
C
            CALL DCOPY(NRZ,W(LGQ),1,W(LD),1)
            CALL DSCAL(NRZ,(-ONE),W(LD),1)
C
            DNORM = DNRM2(NRZ,W(LD),1)
C
            CALL E04NBW(1,N,NRZ,NFREE,LDQ,UNITQ,KX,W(LD),W(LQ),W(LWRK))
            CALL DGEMV('No transpose',NCLIN,N,ONE,A,LDA,W(LD),1,ZERO,
     *                 W(LAD),1)
C
C           ---------------------------------------------------------
C           Find the constraint we bump into along d.
C           Update  x  and  Ax  if the step alfa is nonzero.
C           ---------------------------------------------------------
C           ALFHIT is initialized to BIGALF. If it remains that value
C           after the call to  E04MFS, it is regarded as infinite.
C
            BIGALF = F06BLF(BIGDX,DNORM,OVERFL)
C
            CALL E04MFS(FIRSTV,N,NCLIN,ISTATE,BIGALF,BIGBND,DNORM,
     *                  HITLOW,MOVE,ONBND,UNBNDD,ALFHIT,ALFAP,JADD,
     *                  W(LANORM),W(LAD),AX,BL,BU,FEATOL,FEATLU,W(LD),X)
C
            IF (UNBNDD) THEN
               MSG = 'unbndd'
               GO TO 20
            END IF
C
            ALFA = ALFHIT
            CALL DAXPY(N,ALFA,W(LD),1,X,1)
C
            IF (NCLIN.GT.0) CALL DAXPY(NCLIN,ALFA,W(LAD),1,AX,1)
            XNORM = DNRM2(N,X,1)
C
C           ---------------------------------------------------------
C           Add a constraint to the working set.
C           Update the  TQ  factors of the working set.
C           Use  d  as temporary work space.
C           ---------------------------------------------------------
            IF (BL(JADD).EQ.BU(JADD)) THEN
               ISTATE(JADD) = 3
            ELSE IF (HITLOW) THEN
               ISTATE(JADD) = 1
            ELSE
               ISTATE(JADD) = 2
            END IF
C
            IF (JADD.GT.N) THEN
               IADD = JADD - N
            ELSE
               IF (ALFA.GE.ZERO) THEN
                  IF (HITLOW) THEN
                     X(JADD) = BL(JADD)
                  ELSE
                     X(JADD) = BU(JADD)
                  END IF
               END IF
               DO 40 IFIX = 1, NFREE
                  IF (KX(IFIX).EQ.JADD) GO TO 60
   40          CONTINUE
   60       END IF
C
            CALL E04NFR(UNITQ,RSET,INFORM,IFIX,IADD,JADD,IT,NACTIV,NZ,
     *                  NFREE,NRZ,NGQ,N,LDA,LDQ,LDR,LDT,KX,CONDMX,DZZ,A,
     *                  W(LR),W(LT),W(LGQ),W(LQ),W(LWRK),W(LRLAM),W(LD),
     *                  MSGLVL)
C
            NZ = NZ - 1
            NRZ = NRZ - 1
C
            IF (JADD.LE.N) THEN
C
C              A simple bound has been added.
C
               NFREE = NFREE - 1
            ELSE
C
C              A general constraint has been added.
C
               NACTIV = NACTIV + 1
               KACTIV(NACTIV) = IADD
            END IF
C
C           Increment FEATOL.
C
            CALL DAXPY(NCTOTL,TOLINC,FEATLU,1,FEATOL,1)
C
            IF (MOD(ITER,KCHK).EQ.0) THEN
C              ------------------------------------------------------
C              Check the feasibility of constraints with non-
C              negative ISTATE values.  If some violations have
C              occurred.  Set INFORM to force iterative
C              refinement and a switch to phase 1.
C              ------------------------------------------------------
               CALL E04MFQ(N,NCLIN,ISTATE,BIGBND,NVIOL,JMAX,ERRMAX,AX,
     *                     BL,BU,FEATOL,X)
C
               IF (NVIOL.GT.0) THEN
                  IF (MSGLVL.GT.0) THEN
                     WRITE (REC,FMT=99999) ERRMAX, JMAX
                     CALL X04BAF(IPRINT,REC(1))
                  END IF
               END IF
            END IF
C
            IF (MOD(ITER,IDEGEN).EQ.0) THEN
C
C              Every  IDEGEN  iterations, reset  FEATOL  and
C              move  x  on to the working set if it is close.
C
               CALL E04MFR('End of cycle',MSGLVL,N,NCLIN,NMOVED,ITER,
     *                     NUMINF,ISTATE,BL,BU,FEATOL,FEATLU,X)
C
               NVIOL = NVIOL + NMOVED
            END IF
C
            IF (NVIOL.GT.0) THEN
               MSG = 'resetx'
               GO TO 20
            END IF
C
            IF (NUMINF.NE.0) THEN
               CALL E04MFH(N,NCLIN,LDA,ISTATE,BIGBND,NUMINF,SUMINF,BL,
     *                     BU,A,FEATOL,W(LGQ),X,W(LWTINF))
C
               IF (NUMINF.GT.0) THEN
                  CALL E04NBW(6,N,NZ,NFREE,LDQ,UNITQ,KX,W(LGQ),W(LQ),
     *                        W(LWRK))
               ELSE IF (LP) THEN
                  CALL DCOPY(N,W(LCQ),1,W(LGQ),1)
               END IF
            END IF
C
            IF (NUMINF.EQ.0 .AND. LP) THEN
               OBJ = DDOT(N,CVEC,1,X,1)
            ELSE
               OBJ = SUMINF
            END IF
         END IF
         GO TO 20
C        +    end while
      END IF
C     ======================end of main loop============================
C
      IF (MSG.EQ.'optiml') THEN
         IF (LP) THEN
            IF (NRZ.LT.NZ) THEN
               MSG = 'weak  '
            ELSE
               NTFIXD = 0
               DO 80 J = 1, N
                  IF (ISTATE(J).EQ.4) NTFIXD = NTFIXD + 1
   80          CONTINUE
               IF (NTFIXD.GT.0) MSG = 'weak  '
            END IF
            IF (ABS(JTINY).GT.0) MSG = 'weak  '
         END IF
      ELSE IF (MSG.EQ.'unbndd' .AND. NUMINF.GT.0) THEN
         MSG = 'infeas'
      END IF
C
      RETURN
C
C
C     End of  E04MFZ.  (LPCORE)
C
99999 FORMAT (' XXX  Iterative refinement.  The max violation is ',1P,
     *       D10.2,' in constraint',I5)
      END
