      SUBROUTINE E04NFF(N,NCLIN,A,LDA,BL,BU,CVEC,H,LDH,QPHESS,ISTATE,X,
     *                  ITER,OBJ,AX,CLAMDA,IW,LENIW,W,LENW,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1587 (JUN 1995).
C
C     ******************************************************************
C     E04NFF  solves problems of the form
C
C              minimize               f(x)
C                 x
C                                    (  x )
C              subject to    bl  .le.(    ).ge.  bu,
C                                    ( Ax )
C
C     where  '  denotes the transpose of a column vector,  x  denotes
C     the n-vector of parameters and  f(x) is one of the following...
C
C     FP            =              none    (find a feasible point)
C     LP            =    c'x
C     QP1           =          1/2 x'Hx     H n x n symmetric
C     QP2 (default) =    c'x + 1/2 x'Hx     H n x n symmetric
C     QP3           =          1/2 x'H'Hx   H m x n upper trapezoidal
C     QP4           =    c'x + 1/2 x'H'Hx   H m x n upper trapezoidal
C
C     The matrix  H  is stored in the two-dimensional array  H  of
C     row dimension  LDH.  H  can be entered explicitly as the matrix
C     H,  or implicitly via a user-supplied version of the
C     subroutine QPHESS.  If  LDH = 0,  H is not touched.
C
C     The vector  c  is entered in the one-dimensional array  CVEC.
C
C     nclin  is the number of general linear constraints (rows of  A).
C     (nclin may be zero.)
C
C     The first  n  elements of  bl  and   bu  are lower and upper
C     bounds on the variables.  The next  nclin  elements are
C     lower and upper bounds on the general linear constraints.
C
C     The matrix  A  of coefficients in the general linear constraints
C     is entered as the two-dimensional array  A  (of dimension
C     LDA by n).  If nclin = 0, A is not referenced.
C
C     The vector  x  must contain an initial estimate of the solution,
C     and will contain the computed solution on output.
C
C     Documentation for QPOPT is coming real soon now.
C     Wait for the release of users guide for QPOPT (Version 1.0), by
C     P. E. Gill, W. Murray and M. A. Saunders.
C
C     Version 1.0-6  Jun 30, 1991. (Nag Mk 16 version).
C     Version 1.0-7  Mar 21, 1993. Summary file added.
C     Version 1.0-8  Apr 10, 1994. Sum of infeas. added as an option.
C     Version 1.0-9  Jul 15, 1994. Debug output eliminated.
C
C     This version of  E04NFF  dated 15-Jul-94.
C     Copyright  1988/1994  Optimates.
C     ******************************************************************
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E04NFF')
      INTEGER           LENLC
      PARAMETER         (LENLC=20)
      INTEGER           MXPARM
      PARAMETER         (MXPARM=30)
      DOUBLE PRECISION  ZERO, POINT3, HALF
      PARAMETER         (ZERO=0.0D+0,POINT3=3.3D-1,HALF=0.5D+0)
      DOUBLE PRECISION  POINT8, POINT9, ONE
      PARAMETER         (POINT8=0.8D+0,POINT9=0.9D+0,ONE=1.0D+0)
      DOUBLE PRECISION  HUNDRD
      PARAMETER         (HUNDRD=1.0D+2)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  OBJ
      INTEGER           IFAIL, ITER, LDA, LDH, LENIW, LENW, N, NCLIN
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), AX(*), BL(N+NCLIN), BU(N+NCLIN),
     *                  CLAMDA(N+NCLIN), CVEC(*), H(LDH,*), W(LENW),
     *                  X(N)
      INTEGER           ISTATE(N+NCLIN), IW(LENIW)
C     .. Subroutine Arguments ..
      EXTERNAL          QPHESS
C     .. Scalars in Common ..
      DOUBLE PRECISION  ALFA, ASIZE, BIGBND, BIGDX, BNDLOW, BNDUPP,
     *                  DTMAX, DTMIN, EPSPT3, EPSPT5, EPSPT8, EPSPT9,
     *                  TOLACT, TOLFEA, TOLRNK, TRULAM, TOLINC, TOLX0
      INTEGER           IPRINT, IPRNT, ISDEL, ISUMM, ISUMRY, ITMAX1,
     *                  ITMAX2, JADD, JDEL, KCHK, KCYCLE, LCRASH, LDQ,
     *                  LDT, LENNAM, LINES1, LINES2, LPROB, LQPTYP, M,
     *                  MAXACT, MAXNZ, MINSUM, MM, MSGLC, MXFREE, NCOLT,
     *                  NN, NNCLIN, NOUT, NPROB, IDEGEN, KDEGEN, NDEGEN,
     *                  ITNFIX
      LOGICAL           HEADER, PRNT
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADLC(23), RPSVLC(MXPARM), WMACH(15)
      INTEGER           IPADLC(15), IPSVLC(MXPARM), LOCLC(LENLC),
     *                  NFIX(2)
C     .. Local Scalars ..
      DOUBLE PRECISION  AMIN, CONDMX, EPSMCH, ERRMAX, FEAMAX, FEAMIN,
     *                  HSIZE, RTEPS, XNORM
      INTEGER           IANRMJ, INFORM, IT, ITMAX, J, JINF, JMAX,
     *                  JTHCOL, LANORM, LCQ, LD, LDR, LFEATU, LGQ, LHX,
     *                  LITOTL, LKACTV, LKX, LQ, LR, LRLAM, LT, LWRK,
     *                  LWTINF, LWTOTL, MINACT, MINFXD, MSGLVL, NACT1,
     *                  NACTIV, NARTIF, NCNLN, NCTOTL, NERR, NERROR,
     *                  NFREE, NGQ, NMOVED, NREJTD, NRZ, NUMINF, NVIOL,
     *                  NZ
      LOGICAL           COLD, CSET, DONE, FOUND, HALTED, HOT, NAMED,
     *                  ROWERR, RSET, UNITQ, VERTEX, WARM
      CHARACTER*2       PRBTYP
      CHARACTER*4       START
      CHARACTER*6       MSG
      CHARACTER*11      TITLE
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMLC(MXPARM)
      INTEGER           IPRMLC(MXPARM)
      CHARACTER*8       NAMES(1)
      CHARACTER*80      ERRREC(2), REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2
      INTEGER           P01ABF
      EXTERNAL          DDOT, DNRM2, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DSCAL, E04MFG, E04MFJ, E04MFK,
     *                  E04MFP, E04MFR, E04MFT, E04MFZ, E04NBW, E04NFQ,
     *                  E04NFS, E04NFT, E04NFW, E04NFX, E04NFZ, F06DFF,
     *                  F06FBF, F06FLF, X02ZAZ, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Common blocks ..
      COMMON            /AE04MF/LOCLC
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AE04NF/LQPTYP, M
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
C     .. Data statements ..
      DATA              TITLE/' *** E04NFF'/
C     .. Executable Statements ..
C
C     Set the machine-dependent constants.
C
      CALL X02ZAZ
C
      EPSMCH = WMACH(3)
      RTEPS = WMACH(4)
      NOUT = WMACH(11)
      NERR = WMACH(12)
C
      EPSPT3 = EPSMCH**POINT3
      EPSPT5 = RTEPS
      EPSPT8 = EPSMCH**POINT8
      EPSPT9 = EPSMCH**POINT9
C
      NAMED = .FALSE.
C
      ITER = 0
      HEADER = .TRUE.
      PRNT = .TRUE.
C
C     Set the default values of the parameters.
C
      CALL E04NFW(N,NCLIN,TITLE)
C
      CONDMX = MAX(ONE/EPSPT5,HUNDRD)
C
      LQPTYP = LPROB
      M = MM
      NCTOTL = N + NCLIN
C
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Set all parameters determined by the problem type.
C        key  LQPTYP        objective
C        ---  ------        ---------
C        FP      1              none    (find a feasible point)
C        LP      2    c'x
C        QP1     3          1/2 x'Ax     A n x n symmetric
C        QP2     4    c'x + 1/2 x'Ax     A n x n symmetric
C        QP3     5          1/2 x'A'Ax   A m x n upper trapezoidal
C        QP4     6    c'x + 1/2 x'A'Ax   A m x n upper trapezoidal
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF (LQPTYP.EQ.1) THEN
         PRBTYP = 'FP'
         CSET = .FALSE.
C
      ELSE IF (LQPTYP.EQ.2) THEN
         PRBTYP = 'LP'
         CSET = .TRUE.
C
      ELSE IF (LQPTYP.GE.3 .AND. LQPTYP.LE.6) THEN
         PRBTYP = 'QP'
         CSET = .TRUE.
         IF (LQPTYP.EQ.3 .OR. LQPTYP.EQ.5) CSET = .FALSE.
      ELSE
         PRBTYP = 'illegal'
         MSG = 'noprob'
         GO TO 60
      END IF
C
C     Assign the dimensions of arrays in the parameter list of E04NFZ.
C     Economies of storage are possible if the minimum number of active
C     constraints and the minimum number of fixed variables are known in
C     advance.  The expert user should alter MINACT and MINFXD
C     accordingly.
C     If a linear program is being solved and the matrix of general
C     constraints has fewer rows than columns, i.e.,  nclin .lt. n,  a
C     non-zero value is known for MINFXD.  In this case, VERTEX must
C     be set  .true..
C
      VERTEX = PRBTYP .NE. 'QP' .AND. NCLIN .LT. N
C
      MINFXD = N - MXFREE
      MINACT = MXFREE - MAXNZ
C
      LDT = MAX(MAXNZ,MAXACT)
      NCOLT = MXFREE
      LDR = LDT
      IF (NCLIN.EQ.0) THEN
         LDQ = 1
      ELSE
         LDQ = MAX(1,MXFREE)
      END IF
C
      NCNLN = 0
      LENNAM = 1
C
C     ==================================================================
C     Cold start:  Only  x  is provided.
C     Warm start:  Initial working set is specified in  ISTATE.
C     Hot  start:  The work arrays  IW  and  W  are assumed to have been
C                  initialized during a previous run.
C                  The first four elements of  IW  contain details
C                  on the dimension of the initial working set.
C     ==================================================================
      IF (LCRASH.EQ.0) THEN
         START = 'COLD'
      ELSE IF (LCRASH.EQ.1) THEN
         START = 'WARM'
      ELSE IF (LCRASH.EQ.2) THEN
         START = 'HOT '
      END IF
C
      COLD = LCRASH .EQ. 0
      WARM = LCRASH .EQ. 1
      HOT = LCRASH .EQ. 2
C
C     Allocate remaining work arrays.
C
      LITOTL = 3
      LWTOTL = 0
      CALL E04NFT(CSET,N,NCLIN,LITOTL,LWTOTL)
C
C     Check input parameters and storage limits.
C
      CALL E04MFP(NERROR,MSGLVL,START,LENIW,LENW,LITOTL,LWTOTL,N,NCLIN,
     *            NCNLN,ISTATE,NAMED,NAMES,BIGBND,BL,BU,LDA,LDH,MM,NERR,
     *            LQPTYP,IFAIL)
C
      IF (NERROR.GT.0) THEN
         MSG = 'errors'
         GO TO 60
      END IF
C
      LKACTV = LOCLC(1)
      LKX = LOCLC(2)
C
      LFEATU = LOCLC(3)
      LANORM = LOCLC(4)
      LHX = LOCLC(6)
      LD = LOCLC(7)
      LGQ = LOCLC(8)
      LCQ = LOCLC(9)
      LRLAM = LOCLC(10)
      LR = LOCLC(11)
      LT = LOCLC(12)
      LQ = LOCLC(13)
      LWTINF = LOCLC(14)
      LWRK = LOCLC(15)
C
C     ------------------------------------------------------------------
C     Define the initial feasibility tolerances in CLAMDA.
C     ------------------------------------------------------------------
      IF (TOLFEA.GT.ZERO) CALL F06FBF(N+NCLIN,TOLFEA,W(LFEATU),1)
C
      CALL E04MFR('Initialize anti-cycling variables',MSGLVL,N,NCLIN,
     *            NMOVED,ITER,NUMINF,ISTATE,BL,BU,CLAMDA,W(LFEATU),X)
C
      IF (COLD .OR. WARM) THEN
C        ---------------------------------------------------------------
C        Cold or warm start.  Just about everything must be initialized.
C        The only exception is ISTATE during a warm start.
C        ---------------------------------------------------------------
         IANRMJ = LANORM
         DO 20 J = 1, NCLIN
            W(IANRMJ) = DNRM2(N,A(J,1),LDA)
            IANRMJ = IANRMJ + 1
   20    CONTINUE
         IF (NCLIN.GT.0) CALL F06FLF(NCLIN,W(LANORM),1,ASIZE,AMIN)
C
         CALL F06FLF(NCTOTL,W(LFEATU),1,FEAMAX,FEAMIN)
         CALL DCOPY(NCTOTL,W(LFEATU),1,W(LWTINF),1)
         CALL DSCAL(NCTOTL,(ONE/FEAMIN),W(LWTINF),1)
C
C        ---------------------------------------------------------------
C        Define the initial working set.
C               NFREE ,  NACTIV,  KACTIV, KX,
C               ISTATE (if START  = 'COLD')
C               NARTIF (if VERTEX = 'TRUE')
C        ---------------------------------------------------------------
         CALL E04MFT(START,VERTEX,NCLIN,NCTOTL,NACTIV,NARTIF,NFREE,N,
     *               LDA,ISTATE,IW(LKACTV),IW(LKX),BIGBND,TOLACT,A,AX,
     *               BL,BU,CLAMDA,X,W(LGQ),W(LWRK))
C
C        ---------------------------------------------------------------
C        Compute the TQ factorization of the working set matrix.
C        ---------------------------------------------------------------
         UNITQ = .TRUE.
         NZ = NFREE
         IF (NACTIV.GT.0) THEN
            IT = NACTIV + 1
            NACT1 = NACTIV
            NACTIV = 0
            NGQ = 0
C
            CALL E04NFQ(UNITQ,VERTEX,1,NACT1,IT,NACTIV,NARTIF,NZ,NFREE,
     *                  NREJTD,NGQ,N,LDQ,LDA,LDT,ISTATE,IW(LKACTV),
     *                  IW(LKX),CONDMX,A,W(LT),W(LGQ),W(LQ),W(LWRK),
     *                  W(LD),W(LRLAM),MSGLVL)
         END IF
C
      ELSE IF (HOT) THEN
C        ---------------------------------------------------------------
C        Arrays  IW  and  W  have been defined in a previous run.
C        The first three elements of  IW  are  UNITQ,  NFREE and NACTIV.
C        ---------------------------------------------------------------
         UNITQ = IW(1) .EQ. 1
         NFREE = IW(2)
         NACTIV = IW(3)
C
         NZ = NFREE - NACTIV
      END IF
C
      IF (CSET) THEN
C
C        Install the transformed linear term in CQ.
C
         CALL DCOPY(N,CVEC,1,W(LCQ),1)
         CALL E04NBW(6,N,NZ,NFREE,LDQ,UNITQ,IW(LKX),W(LCQ),W(LQ),W(LWRK)
     *               )
C
      END IF
C
      RSET = .FALSE.
      IF (PRBTYP.EQ.'LP') THEN
         ITMAX = MAX(ITMAX1,ITMAX2)
      ELSE
         ITMAX = ITMAX1
      END IF
C
      JINF = 0
C
C     +    When minimizing the sum of infeasibilities,
C     +    NRZ    =  NZ  implies steepest-descent in the two-norm.
C     +    NRZ    =  0   implies steepest-descent in the infinity norm.
      NRZ = 0
C
C     ==================================================================
C     repeat               (until working set residuals are acceptable)
C     ---------------------------------------------------------------
C     Move x onto the constraints in the working set.
C     ---------------------------------------------------------------
   40 CALL E04MFJ(ROWERR,UNITQ,NCLIN,NACTIV,NFREE,NZ,N,LDQ,LDA,LDT,
     *            ISTATE,IW(LKACTV),IW(LKX),JMAX,ERRMAX,XNORM,A,AX,BL,
     *            BU,CLAMDA,W(LT),X,W(LQ),W(LD),W(LWRK))
C
      IF (ROWERR) THEN
         MSG = 'rowerr'
         NUMINF = 1
         GO TO 60
      END IF
C
      CALL E04MFZ(PRBTYP,MSG,CSET,RSET,UNITQ,ITER,ITMAX,JINF,NVIOL,N,
     *            NCLIN,LDA,NACTIV,NFREE,NRZ,NZ,ISTATE,IW(LKACTV),
     *            IW(LKX),E04NFS,OBJ,NUMINF,XNORM,A,AX,BL,BU,CVEC,
     *            CLAMDA,W(LFEATU),X,W)
C
      IF (PRBTYP.EQ.'QP' .AND. MSG.EQ.'feasbl') THEN
         IF (MSGLVL.EQ.5 .OR. MSGLVL.GE.10) THEN
            WRITE (REC,FMT=99999) ITER
            CALL X04BAF(IPRINT,REC(1))
         END IF
         IF (MSGLVL.GE.5) THEN
            IF (ISUMM.GE.0 .AND. ISUMM.NE.IPRINT) THEN
               WRITE (REC,FMT=99999) ITER
               CALL X04BAF(ISUMM,REC(1))
            END IF
         END IF
         RSET = .TRUE.
         ITMAX = ITER + ITMAX2
C
C        --------------------------------------------------------
C        Compute the first QP objective and transformed gradient.
C        --------------------------------------------------------
         IF (CSET) THEN
            OBJ = DDOT(N,CVEC,1,X,1)
         ELSE
            OBJ = ZERO
         END IF
C
         JTHCOL = 0
         CALL QPHESS(N,JTHCOL,H,LDH,X,W(LHX))
         OBJ = OBJ + HALF*DDOT(N,W(LHX),1,X,1)
         CALL DCOPY(N,W(LHX),1,W(LGQ),1)
         CALL E04NBW(6,N,NZ,NFREE,LDQ,UNITQ,IW(LKX),W(LGQ),W(LQ),W(LWRK)
     *               )
         IF (CSET) CALL DAXPY(N,ONE,W(LCQ),1,W(LGQ),1)
C
C        ------------------------------------------------------------
C        Compute the Cholesky factor R of an initial reduced Hessian.
C        The magnitudes of the diagonals of  R  are nonincreasing.
C        ------------------------------------------------------------
         IF (CSET) THEN
            NGQ = 2
         ELSE
            NGQ = 1
         END IF
         HSIZE = ONE
         CALL E04NFX(UNITQ,QPHESS,MAXNZ,N,NGQ,NRZ,NZ,NFREE,LDQ,LDH,LDR,
     *               IW(LKX),HSIZE,TOLRNK,W(LGQ),H,W(LR),W(LQ),W(LWRK),
     *               W(LRLAM))
C
         CALL E04NFZ(PRBTYP,MSG,CSET,UNITQ,ITER,ITMAX,NVIOL,N,NCLIN,LDA,
     *               LDH,NACTIV,NFREE,NRZ,NZ,ISTATE,IW(LKACTV),IW(LKX),
     *               QPHESS,E04NFS,OBJ,XNORM,HSIZE,A,AX,BL,BU,CVEC,
     *               CLAMDA,W(LFEATU),H,X,W)
      END IF
C
      FOUND = MSG .EQ. 'optiml' .OR. MSG .EQ. 'feasbl' .OR. MSG .EQ.
     *        'deadpt' .OR. MSG .EQ. 'weak  ' .OR. MSG .EQ.
     *        'unbndd' .OR. MSG .EQ. 'infeas'
      HALTED = MSG .EQ. 'itnlim' .OR. MSG .EQ. 'Rz2big'
C
      IF (FOUND) THEN
         CALL E04MFR('Optimal',MSGLVL,N,NCLIN,NMOVED,ITER,NUMINF,ISTATE,
     *               BL,BU,CLAMDA,W(LFEATU),X)
      END IF
C
      DONE = FOUND .AND. NVIOL .EQ. 0 .AND. NMOVED .EQ. 0
C
C     until      done  .or.  halted
      IF ( .NOT. (DONE .OR. HALTED)) GO TO 40
C     ==================================================================
C     Set   CLAMDA.  Print the full solution.
C     Clean up.  Save values for a subsequent hot start.
C     ------------------------------------------------------------------
C
      CALL E04MFG(NFREE,LDA,N,NCLIN,NCTOTL,NACTIV,ISTATE,IW(LKACTV),
     *            IW(LKX),A,BL,BU,X,CLAMDA,W(LFEATU),W(LWRK),W(LRLAM),X)
      CALL E04MFK(MSGLVL,N,NCLIN,NCTOTL,BIGBND,NAMED,NAMES,ISTATE,BL,BU,
     *            CLAMDA,W(LFEATU),W(LWRK))
C
      IW(1) = 0
      IF (UNITQ) IW(1) = 1
      IW(2) = NFREE
      IW(3) = NACTIV
C
C     ==================================================================
C     Print messages if required.
C     Recover the optional parameters set by the user.
C     ==================================================================
   60 IF (MSG.EQ.'feasbl') THEN
         INFORM = 0
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99997)
            CALL X04BAY(IPRINT,2,REC)
         END IF
C
      ELSE IF (MSG.EQ.'optiml') THEN
         INFORM = 0
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99996) PRBTYP
            CALL X04BAY(IPRINT,2,REC)
         END IF
C
      ELSE IF (MSG.EQ.'deadpt' .OR. MSG.EQ.'weak  ') THEN
         INFORM = 1
         IF (MSGLVL.GT.0) THEN
            IF (PRBTYP.EQ.'QP') THEN
               WRITE (REC,FMT=99995)
               CALL X04BAY(IPRINT,3,REC)
               IF (NZ.GT.NRZ) THEN
                  WRITE (REC,FMT=99991) NZ - NRZ
                  CALL X04BAF(IPRINT,REC(1))
               END IF
            ELSE
               WRITE (REC,FMT=99993)
               CALL X04BAY(IPRINT,2,REC)
            END IF
         END IF
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) THEN
            IF (PRBTYP.EQ.'QP') THEN
               WRITE (ERRREC,FMT=99994)
            ELSE
               WRITE (ERRREC,FMT=99992)
            END IF
         END IF
C
      ELSE IF (MSG.EQ.'unbndd') THEN
         INFORM = 2
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99990) PRBTYP
            CALL X04BAY(IPRINT,2,REC)
         END IF
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) WRITE (ERRREC,FMT=99989)
     *       PRBTYP
C
      ELSE IF (MSG.EQ.'infeas') THEN
         INFORM = 3
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99988)
            CALL X04BAY(IPRINT,2,REC)
         END IF
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) WRITE (ERRREC,FMT=99987)
C
      ELSE IF (MSG.EQ.'rowerr') THEN
         INFORM = 3
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99975)
            CALL X04BAY(IPRINT,3,REC)
         END IF
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) WRITE (ERRREC,FMT=99974)
C
      ELSE IF (MSG.EQ.'itnlim') THEN
         INFORM = 4
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99986)
            CALL X04BAY(IPRINT,2,REC)
         END IF
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) WRITE (ERRREC,FMT=99985)
C
      ELSE IF (MSG.EQ.'Rz2big') THEN
         INFORM = 5
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99984) MAXNZ
            CALL X04BAY(IPRINT,3,REC)
         END IF
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) WRITE (ERRREC,FMT=99983)
C
      ELSE IF (MSG.EQ.'errors') THEN
         INFORM = 6
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99982) NERROR
            CALL X04BAY(IPRINT,2,REC)
         END IF
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) WRITE (ERRREC,FMT=99981)
     *       NERROR
C
      ELSE IF (MSG.EQ.'noprob') THEN
         INFORM = 7
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99980)
            CALL X04BAY(IPRINT,2,REC)
         END IF
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) WRITE (ERRREC,FMT=99979)
      END IF
C
      IF (MSGLVL.GT.0) THEN
C
         IF (INFORM.LT.5) THEN
            IF (NUMINF.EQ.0) THEN
               IF (PRBTYP.NE.'FP') THEN
                  WRITE (REC,FMT=99978) PRBTYP, OBJ
                  CALL X04BAY(IPRINT,2,REC)
               END IF
            ELSE IF (INFORM.EQ.3) THEN
               IF (MSG.EQ.'infeas') THEN
                  IF (MINSUM.EQ.0) THEN
                     WRITE (REC,FMT=99973) OBJ
                     CALL X04BAY(IPRINT,2,REC)
                  ELSE
                     WRITE (REC,FMT=99977) OBJ
                     CALL X04BAY(IPRINT,2,REC)
                  END IF
               ELSE IF (MSG.EQ.'rowerr') THEN
                  WRITE (REC,FMT=99972) ERRMAX
                  CALL X04BAY(IPRINT,2,REC)
               END IF
            ELSE
               WRITE (REC,FMT=99976) OBJ
               CALL X04BAY(IPRINT,2,REC)
            END IF
         END IF
      END IF
C
      IF (INFORM.LT.6) THEN
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99998) PRBTYP, ITER
            CALL X04BAY(IPRINT,2,REC)
         END IF
      END IF
C
      CALL F06DFF(MXPARM,IPSVLC,1,IPRMLC,1)
      CALL DCOPY(MXPARM,RPSVLC,1,RPRMLC,1)
C
      IF ((INFORM.GE.1 .AND. INFORM.LE.7)
     *    .AND. (IFAIL.EQ.0 .OR. IFAIL.EQ.-1)) CALL X04BAY(NERR,2,
     *    ERRREC)
      IFAIL = P01ABF(IFAIL,INFORM,SRNAME,0,REC)
      RETURN
C
C
C
C     End of  E04NFF.  (QPOPT)
C
99999 FORMAT (' Itn',I6,' -- Feasible point found.')
99998 FORMAT (/' Exit from ',A2,' problem after ',I5,' iterations.')
99997 FORMAT (/' Exit E04NFF - Feasible point found.')
99996 FORMAT (/' Exit E04NFF - Optimal ',A2,' solution.')
99995 FORMAT (/' Exit E04NFF - Iterations terminated at a dead-point',
     *       /15X,'(check the optimality conditions).     ')
99994 FORMAT (/' ** Iterations terminated at a dead-point.')
99993 FORMAT (/' Exit E04NFF - Optimal solution is not unique.')
99992 FORMAT (/' ** Optimal solution is not unique.')
99991 FORMAT ('             - Artificial constraints in working set = ',
     *       I4,'.')
99990 FORMAT (/' Exit E04NFF - ',A2,' solution is unbounded.')
99989 FORMAT (/' ** ',A2,' solution is unbounded.')
99988 FORMAT (/' Exit E04NFF - No feasible point for the linear constr',
     *       'aints.')
99987 FORMAT (/' ** No feasible point for the linear constraints.')
99986 FORMAT (/' Exit E04NFF - Too many iterations.')
99985 FORMAT (/' ** Too many iterations.')
99984 FORMAT (/' Exit E04NFF - Reduced Hessian exceeds assigned dimens',
     *       'ion.',/'             - Maximum degrees of freedom = ',I4,
     *       '.')
99983 FORMAT (/' ** Reduced Hessian exceeds assigned dimension.')
99982 FORMAT (/' Exit E04NFF - ',I7,' errors found in the input parame',
     *       'ters.  Problem abandoned.')
99981 FORMAT (/' ** ',I7,' errors found in the input parameters.  Prob',
     *       'lem abandoned.')
99980 FORMAT (/' Exit E04NFF - Problem type not recognized.  Problem a',
     *       'bandoned.')
99979 FORMAT (/' ** Problem type not recognized.  Problem abandoned.')
99978 FORMAT (/' Final ',A2,' objective value =',G16.7)
99977 FORMAT (/' Minimum sum of infeasibilities =',G16.7)
99976 FORMAT (/' Final sum of infeasibilities =',G16.7)
99975 FORMAT (/' Exit E04NFF - Cannot satisfy the working set constrai',
     *       'nts to the accuracy',/'requested.')
99974 FORMAT (/' ** Cannot satisfy the working set constraints to the ',
     *       'accuracy requested.')
99973 FORMAT (/' Sum of infeasibilities =',G16.7)
99972 FORMAT (/' Maximum row error =',G16.7)
      END
