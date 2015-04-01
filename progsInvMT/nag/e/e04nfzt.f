      SUBROUTINE E04NFZ(PRBTYP,MSG,CSET,UNITQ,ITER,ITMAX,NVIOL,N,NCLIN,
     *                  LDA,LDH,NACTIV,NFREE,NRZ,NZ,ISTATE,KACTIV,KX,
     *                  QPHESS,E04NFS,OBJQP,XNORM,HSIZE,A,AX,BL,BU,CVEC,
     *                  FEATOL,FEATLU,H,X,W)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1598 (JUN 1995).
C
C     ******************************************************************
C     E04NFZ  is a subroutine for general quadratic programming.
C     On entry, it is assumed that an initial working set of
C     linear constraints and bounds is available.
C     The arrays  ISTATE, KACTIV and KX  will have been set accordingly
C     and the arrays  T  and  Q  will contain the TQ factorization of
C     the matrix whose rows are the gradients of the active linear
C     constraints with the columns corresponding to the active bounds
C     removed.  The TQ factorization of the resulting (NACTIV by NFREE)
C     matrix is  A(free)*Q = (0 T),  where Q is (NFREE by NFREE) and T
C     is upper-triangular.
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
C     This version of  E04NFZ  dated 11-Nov-92.
C     Copyright  1988--1994  Optimates.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LENLC
      PARAMETER         (LENLC=20)
      INTEGER           MXPARM
      PARAMETER         (MXPARM=30)
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER         (ZERO=0.0D+0,HALF=0.5D+0,ONE=1.0D+0)
      CHARACTER*6       EMPTY
      PARAMETER         (EMPTY='      ')
      INTEGER           MREFN
      PARAMETER         (MREFN=1)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HSIZE, OBJQP, XNORM
      INTEGER           ITER, ITMAX, LDA, LDH, N, NACTIV, NCLIN, NFREE,
     *                  NRZ, NVIOL, NZ
      LOGICAL           CSET, UNITQ
      CHARACTER*2       PRBTYP
      CHARACTER*6       MSG
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), AX(*), BL(N+NCLIN), BU(N+NCLIN),
     *                  CVEC(*), FEATLU(N+NCLIN), FEATOL(N+NCLIN),
     *                  H(LDH,*), W(*), X(N)
      INTEGER           ISTATE(N+NCLIN), KACTIV(N), KX(N)
C     .. Subroutine Arguments ..
      EXTERNAL          E04NFS, QPHESS
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
     *                  CONDT, DINKY, DNORM, DRZMAX, DRZMIN, DRZZ,
     *                  ERRMAX, FLMAX, GFNORM, GRZNRM, GZDZ, GZNORM,
     *                  OBJCHG, OBJSIZ, SMLLST, SUMINF, TINYST, TRUBIG,
     *                  TRUSML, WSSIZE, ZEROLM
      INTEGER           IADD, IFIX, INFORM, IREFN, ISSAVE, IT, JBIGST,
     *                  JDSAVE, JINF, JMAX, JSMLST, JTHCOL, JTINY,
     *                  KBIGST, KDEL, KSMLST, LAD, LANORM, LCQ, LD, LDR,
     *                  LGQ, LHX, LQ, LR, LRLAM, LT, LWRK, LWTINF,
     *                  MSGLVL, NCTOTL, NGQ, NMOVED, NOTOPT, NUMINF
      LOGICAL           DEADPT, DELREG, FIRSTV, GIVEUP, HITCON, HITLOW,
     *                  MINMZR, MOVE, ONBND, OVERFL, POSDEF, RENEWR,
     *                  RSET, SINGLR, STATPT, UNBNDD, UNCON, UNITGZ
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMLC(MXPARM)
      INTEGER           IPRMLC(MXPARM)
      CHARACTER*80      REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2, F06BLF
      EXTERNAL          DDOT, DNRM2, F06BLF
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, E04MFL, E04MFM, E04MFQ, E04MFR,
     *                  E04MFS, E04NBW, E04NFP, E04NFR, E04NFV, E04NFY,
     *                  F06FLF, X04BAF
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
      SAVE              /AX02ZA/, /AE04MF/, /FE04MF/, /GE04MF/
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
      LDR = LDT
      IT = 1
C
      LANORM = LOCLC(4)
      LAD = LOCLC(5)
      LHX = LOCLC(6)
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
C     Initialize.
C
      IREFN = 0
      JINF = 0
      NCTOTL = N + NCLIN
      NVIOL = 0
      NUMINF = 0
      SUMINF = ZERO
      CONDMX = FLMAX
C
      DELREG = .FALSE.
      FIRSTV = .FALSE.
      POSDEF = .TRUE.
      RENEWR = .FALSE.
      RSET = .TRUE.
      SINGLR = .FALSE.
      UNCON = .FALSE.
      UNITGZ = .FALSE.
C
      NOTOPT = 0
      DRZZ = ONE
C
      MSG = EMPTY
C
C*    ======================Start of main loop==========================
C     +    do while (msg .eq. empty)
   20 IF (MSG.EQ.EMPTY) THEN
C
         IF (NZ.GT.0) THEN
            GZNORM = DNRM2(NZ,W(LGQ),1)
         ELSE
            GZNORM = ZERO
         END IF
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
         OBJSIZ = ONE + ABS(OBJQP)
         WSSIZE = ZERO
         IF (NACTIV.GT.0) WSSIZE = DTMAX
         DINKY = EPSPT8*MAX(WSSIZE,OBJSIZ,GFNORM)
         IF (UNCON) THEN
            UNITGZ = GRZNRM .LE. DINKY
         END IF
C
C        If the reduced gradient Z'g is small and Hz is positive
C        definite,  x is a minimizer on the working set.
C        A maximum number of unconstrained steps is imposed to
C        allow for  DINKY  being too large because of bad scaling.
C
         STATPT = GRZNRM .LE. DINKY
         GIVEUP = IREFN .GT. MREFN
C
         MINMZR = STATPT .AND. POSDEF
         DEADPT = STATPT .AND. SINGLR
C
C        ------------------------------------------------------------
C        Print the details of this iteration.
C        ------------------------------------------------------------
C        Define small quantities that reflect the size of x, R and
C        the constraints in the working set.
C
         IF (PRNT) THEN
            IF (NRZ.GT.0) THEN
               CALL F06FLF(NRZ,W(LR),LDR+1,DRZMAX,DRZMIN)
               CONDRZ = F06BLF(DRZMAX,DRZMIN,OVERFL)
            ELSE
               CONDRZ = ONE
            END IF
C
            IF (NACTIV.GT.0) THEN
               CONDT = F06BLF(DTMAX,DTMIN,OVERFL)
            ELSE
               CONDT = ONE
            END IF
C
            CALL E04NFS(PRBTYP,HEADER,RSET,MSGLVL,ITER,ISDEL,JDEL,JADD,
     *                  N,NCLIN,NACTIV,NFREE,NZ,NRZ,LDR,LDT,ISTATE,ALFA,
     *                  CONDRZ,CONDT,DRZZ,GRZNRM,NUMINF,SUMINF,NOTOPT,
     *                  OBJQP,TRULAM,AX,W(LR),W(LT),X,W(LWRK))
         END IF
C
         IF (MINMZR .OR. GIVEUP) THEN
C           =========================================================
C           The point  x  is a constrained stationary point.
C           Compute Lagrange multipliers.
C           =========================================================
C           Define what we mean by ``non-optimal'' multipliers.
C
            NOTOPT = 0
            JDEL = 0
            ZEROLM = DINKY
            SMLLST = DINKY
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
            IF (NOTOPT.EQ.0 .AND. POSDEF) THEN
               MSG = 'optiml'
               GO TO 20
            END IF
C
C           ---------------------------------------------------------
C           Delete one of three types of constraint
C           (1) regular           JSMLST > 0   ISTATE(JSMLST) = 1, 2
C           (2) temporary bound   JSMLST > 0,  ISTATE(JSMLST) = 4
C           (3) artificial        JSMLST < 0
C           ---------------------------------------------------------
            TRULAM = TRUSML
            JDEL = JSMLST
            DELREG = .FALSE.
C
            IF (NRZ+1.GT.MAXNZ .AND. JDEL.NE.0) THEN
               MSG = 'Rz2big'
               GO TO 20
            END IF
C
            IF (JDEL.GT.0) THEN
C
C              Regular constraint or temporary bound.
C              DELREG  says that a regular constraint was deleted.
C              JDSAVE, ISSAVE are only defined if  DELREG  is true.
C
               KDEL = KSMLST
               ISDEL = ISTATE(JDEL)
               ISTATE(JDEL) = 0
               DELREG = ISDEL .NE. 4
               IF (DELREG) THEN
                  JDSAVE = JDEL
                  ISSAVE = ISDEL
               END IF
            END IF
C
C           Update the factorizations.
C
            CALL E04NFP(UNITQ,IT,N,NACTIV,NFREE,NGQ,NZ,NRZ,LDA,LDQ,LDT,
     *                  JDEL,KDEL,KACTIV,KX,A,W(LT),W(LGQ),W(LQ),W(LD),
     *                  W(LRLAM))
C
            RENEWR = .TRUE.
            CALL E04NFY(SINGLR,POSDEF,RENEWR,UNITQ,N,NRZ,NFREE,LDQ,LDH,
     *                  LDR,KX,HSIZE,DRZZ,TOLRNK,QPHESS,H,W(LR),W(LQ),
     *                  W(LWRK),W(LD))
C
            IREFN = 0
            PRNT = .FALSE.
            UNCON = .FALSE.
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
            CALL E04NFV(DELREG,POSDEF,STATPT,UNITGZ,UNITQ,N,NCLIN,NFREE,
     *                  LDA,LDQ,LDR,NRZ,ISSAVE,JDSAVE,KX,DNORM,GZDZ,A,
     *                  W(LAD),W(LD),W(LGQ),W(LR),W(LQ),W(LWRK))
C
C           ---------------------------------------------------------
C           Find the constraint we bump into along  d.
C           Update  x  and  Ax  if the step  ALFA  is nonzero.
C           ---------------------------------------------------------
C           E04MFS initializes  ALFHIT  to BIGALF. If it is still
C           that value on exit,  it is regarded as infinite.
C
            BIGALF = F06BLF(BIGDX,DNORM,OVERFL)
C
            CALL E04MFS(FIRSTV,N,NCLIN,ISTATE,BIGALF,BIGBND,DNORM,
     *                  HITLOW,MOVE,ONBND,UNBNDD,ALFHIT,ALFAP,JADD,
     *                  W(LANORM),W(LAD),AX,BL,BU,FEATOL,FEATLU,W(LD),X)
C
C           ---------------------------------------------------------
C           If Hz is positive definite,  ALFA = 1.0  will be the step
C           to the minimizer of the quadratic on the current working
C           set.  If the unit step does not violate the nearest
C           constraint by more than FEATOL,  the constraint is not
C           added to the working set.
C           ---------------------------------------------------------
            UNCON = ALFAP .GT. ONE .AND. POSDEF
            HITCON = .NOT. UNCON
C
            IF (HITCON) THEN
               ALFA = ALFHIT
               IREFN = 0
            ELSE
               IREFN = IREFN + 1
               JADD = 0
               ALFA = ONE
            END IF
C
            IF (HITCON .AND. UNBNDD) THEN
               MSG = 'unbndd'
               GO TO 20
            END IF
C
C           Predict the change in the QP objective function.
C
            IF (POSDEF) THEN
               OBJCHG = ALFA*GZDZ*(ONE-HALF*ALFA)
            ELSE
               OBJCHG = ALFA*GZDZ + HALF*ALFA**2*DRZZ
            END IF
C
C           Check for a dead point or unbounded solution.
C
            IF (OBJCHG.GE.-EPSPT9*OBJSIZ .AND. DEADPT) THEN
               MSG = 'deadpt'
               GO TO 20
            END IF
C
            IF (OBJCHG.GE.EPSPT9*OBJSIZ) THEN
               MSG = 'resetx'
               GO TO 20
            END IF
C
            CALL DAXPY(N,ALFA,W(LD),1,X,1)
            IF (NCLIN.GT.0) CALL DAXPY(NCLIN,ALFA,W(LAD),1,AX,1)
            XNORM = DNRM2(N,X,1)
C
            IF (HITCON) THEN
C              ------------------------------------------------------
C              Add a constraint to the working set.
C              Update the TQ factors of the working set.
C              Use  d  as temporary work space.
C              ------------------------------------------------------
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
                  IF (HITLOW) THEN
                     X(JADD) = BL(JADD)
                  ELSE
                     X(JADD) = BU(JADD)
                  END IF
C
                  DO 40 IFIX = 1, NFREE
                     IF (KX(IFIX).EQ.JADD) GO TO 60
   40             CONTINUE
   60          END IF
C
               CALL E04NFR(UNITQ,RSET,INFORM,IFIX,IADD,JADD,IT,NACTIV,
     *                     NZ,NFREE,NRZ,NGQ,N,LDA,LDQ,LDR,LDT,KX,CONDMX,
     *                     DRZZ,A,W(LR),W(LT),W(LGQ),W(LQ),W(LWRK),
     *                     W(LRLAM),W(LD),MSGLVL)
               NRZ = NRZ - 1
               NZ = NZ - 1
C
               IF (JADD.LE.N) THEN
C
C                 A simple bound has been added.
C
                  NFREE = NFREE - 1
               ELSE
C
C                 A general constraint has been added.
C
                  NACTIV = NACTIV + 1
                  KACTIV(NACTIV) = IADD
               END IF
C
C              ------------------------------------------------------
C              Check if  Hz  has become positive definite.
C              Recompute the last column of Rz if unacceptable
C              growth has occurred.
C              -----------------------------------------------------
               IF ( .NOT. POSDEF) THEN
                  CALL E04NFY(SINGLR,POSDEF,RENEWR,UNITQ,N,NRZ,NFREE,
     *                        LDQ,LDH,LDR,KX,HSIZE,DRZZ,TOLRNK,QPHESS,H,
     *                        W(LR),W(LQ),W(LWRK),W(LD))
               END IF
            END IF
C
C           Increment FEATOL.
C
            CALL DAXPY(NCTOTL,TOLINC,FEATLU,1,FEATOL,1)
C
            IF (MOD(ITER,KCHK).EQ.0) THEN
C              ------------------------------------------------------
C              Check the feasibility of constraints with non-
C              negative  ISTATE  values.  If violations have
C              occurred,  force iterative refinement and a switch
C              to phase 1.
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
               NVIOL = NVIOL + NMOVED
            END IF
C
            IF (NVIOL.GT.0) THEN
               MSG = 'resetx'
               GO TO 20
            END IF
C
C           ---------------------------------------------------------
C           Compute the QP objective and transformed gradient.
C           ------------------------------------------------------
            IF (CSET) THEN
               OBJQP = DDOT(N,CVEC,1,X,1)
            ELSE
               OBJQP = ZERO
            END IF
C
            JTHCOL = 0
            CALL QPHESS(N,JTHCOL,H,LDH,X,W(LHX))
            OBJQP = OBJQP + HALF*DDOT(N,W(LHX),1,X,1)
            CALL DCOPY(N,W(LHX),1,W(LGQ),1)
            CALL E04NBW(6,N,NZ,NFREE,LDQ,UNITQ,KX,W(LGQ),W(LQ),W(LWRK))
            IF (CSET) CALL DAXPY(N,ONE,W(LCQ),1,W(LGQ),1)
         END IF
         GO TO 20
C        +    end while
      END IF
C     ======================end of main loop============================
C
      RETURN
C
C
C     End of  E04NFZ.  (QPCORE)
C
99999 FORMAT (' XXX  Iterative refinement.  The max violation is ',1P,
     *       D10.2,' in constraint',I5)
      END
