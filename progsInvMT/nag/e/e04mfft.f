      SUBROUTINE E04MFF(N,NCLIN,A,LDA,BL,BU,CVEC,ISTATE,X,ITER,OBJ,AX,
     *                  CLAMDA,IW,LENIW,W,LENW,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1556 (JUN 1995).
C
C     ******************************************************************
C     E04MFF  solves the linear programming problem
C
C           minimize               c' x
C              x
C                                 (  x )
C           subject to    bl  .le.(    ).ge.  bu,
C                                 ( Ax )
C
C     where  A  is a constant  nclin by n  matrix.
C     The feasible region is defined by a mixture of linear equality or
C     inequality constraints on  x.
C
C     n  is the number of variables (dimension of x).
C        (n must be positive.)
C
C     nclin  is the number of general linear constraints (rows of  A).
C        (nclin may be zero.)
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
C     Documentation for  E04MFF  is coming real soon now.
C     Wait for the release of  users guide for LPOPT (Version 1.00-6),
C     by P. E. Gill, W. Murray and M. A. Saunders,
C
C     Version 1.0-6    Jun 30, 1991.  (Nag Mk 16 version.)
C     Version 1.0-7    Mar 21, 1993.  Summary file added.
C     Version 1.0-8    Apr 10, 1994.  Sum of infeas. added as an option.
C     Version 1.0-9    Jul 15, 1994.  Debug printing removed.
C
C     Copyright  1989/1994  Optimates.
C     ******************************************************************
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E04MFF')
      INTEGER           LENLC
      PARAMETER         (LENLC=20)
      INTEGER           MXPARM
      PARAMETER         (MXPARM=30)
      DOUBLE PRECISION  ZERO, POINT3
      PARAMETER         (ZERO=0.0D+0,POINT3=3.3D-1)
      DOUBLE PRECISION  POINT8, POINT9, ONE
      PARAMETER         (POINT8=0.8D+0,POINT9=0.9D+0,ONE=1.0D+0)
      DOUBLE PRECISION  HUNDRD
      PARAMETER         (HUNDRD=1.0D+2)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  OBJ
      INTEGER           IFAIL, ITER, LDA, LENIW, LENW, N, NCLIN
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), AX(*), BL(N+NCLIN), BU(N+NCLIN),
     *                  CLAMDA(N+NCLIN), CVEC(*), W(LENW), X(N)
      INTEGER           ISTATE(N+NCLIN), IW(LENIW)
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
      DOUBLE PRECISION  AMIN, CONDMX, EPSMCH, ERRMAX, FEAMAX, FEAMIN,
     *                  RTEPS, XNORM
      INTEGER           IANRMJ, INFORM, IT, ITMAX, J, JINF, JMAX,
     *                  LANORM, LCQ, LD, LDH, LFEATU, LGQ, LITOTL,
     *                  LKACTV, LKX, LLPTYP, LQ, LRLAM, LT, LWRK,
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
      CHARACTER*80      ERRREC(2), REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2
      INTEGER           P01ABF
      EXTERNAL          DNRM2, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DSCAL, E04MFG, E04MFJ, E04MFK, E04MFP,
     *                  E04MFR, E04MFT, E04MFU, E04MFV, E04MFW, E04MFZ,
     *                  E04NBW, E04NFQ, F06DFF, F06FBF, F06FLF, X02ZAZ,
     *                  X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
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
C     .. Data statements ..
      DATA              TITLE/' *** E04MFF'/
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
      INFORM = 0
      ITER = 0
C
      HEADER = .TRUE.
      PRNT = .TRUE.
C
      CONDMX = MAX(ONE/EPSPT5,HUNDRD)
C
C     Set the default values of the parameters.
C
      CALL E04MFW(N,NCLIN,TITLE)
C
      LLPTYP = LPROB
      NCTOTL = N + NCLIN
C
C     Set all parameters determined by the problem type.
C
      IF (LLPTYP.EQ.1) THEN
         PRBTYP = 'FP'
         CSET = .FALSE.
      ELSE IF (LLPTYP.EQ.2) THEN
         PRBTYP = 'LP'
         CSET = .TRUE.
      ELSE
         PRBTYP = 'illegal'
         MSG = 'noprob'
         GO TO 60
      END IF
C
C     Assign the dimensions of arrays in the parameter list of E04MFZ.
C     Economies of storage are possible if the minimum number of active
C     constraints and the minimum number of fixed variables are known in
C     advance.  The expert user should alter MINACT and MINFXD
C     accordingly.
C     If a linear program is being solved and the matrix of general
C     constraints has fewer rows than columns, i.e.,  nclin .lt. n,
C     a non-zero value is known for MINFXD.  Note that in this case,
C     VERTEX must be set  .true..
C
      VERTEX = NCLIN .LT. N
C
      MINFXD = N - MXFREE
      MINACT = MXFREE - MAXNZ
C
      LDT = MAX(MAXNZ,MAXACT)
      NCOLT = MXFREE
      IF (NCLIN.EQ.0) THEN
         LDQ = 1
      ELSE
         LDQ = MAX(1,MXFREE)
      END IF
C
      NCNLN = 0
      LENNAM = 1
      LDH = 1
      MM = 0
C
C     ==================================================================
C     Cold start:  Only  x  is provided.
C     Warm start:  Initial working set is specified in  ISTATE.
C     Hot  start:  The work arrays  IW  and  W  are assumed to have been
C                  initialized during a previous run.
C                  The first three elements of  IW  contain details
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
      CALL E04MFV(CSET,N,NCLIN,LITOTL,LWTOTL)
C
C     Check input parameters and storage limits.
C
      CALL E04MFP(NERROR,MSGLVL,START,LENIW,LENW,LITOTL,LWTOTL,N,NCLIN,
     *            NCNLN,ISTATE,NAMED,NAMES,BIGBND,BL,BU,LDA,LDH,MM,NERR,
     *            LLPTYP,IFAIL)
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
      LD = LOCLC(7)
      LGQ = LOCLC(8)
      LCQ = LOCLC(9)
      LRLAM = LOCLC(10)
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
C
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
      END IF
C
      RSET = .FALSE.
      ITMAX = ITMAX2
      JINF = 0
C
C     +    Take your pick when minimizing the sum of infeasibilities:
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
     *            BU,W(LFEATU),W(LT),X,W(LQ),W(LD),W(LWRK))
C
      IF (ROWERR) THEN
         MSG = 'rowerr'
         NUMINF = 1
         GO TO 60
      END IF
C
      CALL E04MFZ(PRBTYP,MSG,CSET,RSET,UNITQ,ITER,ITMAX,JINF,NVIOL,N,
     *            NCLIN,LDA,NACTIV,NFREE,NRZ,NZ,ISTATE,IW(LKACTV),
     *            IW(LKX),E04MFU,OBJ,NUMINF,XNORM,A,AX,BL,BU,CVEC,
     *            CLAMDA,W(LFEATU),X,W)
C
      FOUND = MSG .EQ. 'feasbl' .OR. MSG .EQ. 'optiml' .OR. MSG .EQ.
     *        'weak  ' .OR. MSG .EQ. 'unbndd' .OR. MSG .EQ. 'infeas'
      HALTED = MSG .EQ. 'itnlim'
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
C     ===========================================================
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
   60 IF (MSG.EQ.'optiml') THEN
         INFORM = 0
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99997) PRBTYP
            CALL X04BAY(IPRINT,2,REC)
         END IF
C
      ELSE IF (MSG.EQ.'feasbl') THEN
         INFORM = 0
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99998)
            CALL X04BAY(IPRINT,2,REC)
         END IF
C
      ELSE IF (MSG.EQ.'weak  ') THEN
         INFORM = 1
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99996) PRBTYP
            CALL X04BAY(IPRINT,2,REC)
         END IF
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) WRITE (ERRREC,FMT=99995)
     *       PRBTYP
C
      ELSE IF (MSG.EQ.'unbndd') THEN
         INFORM = 2
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99994) PRBTYP
            CALL X04BAY(IPRINT,2,REC)
         END IF
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) WRITE (ERRREC,FMT=99993)
     *       PRBTYP
C
      ELSE IF (MSG.EQ.'infeas') THEN
         INFORM = 3
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99992)
            CALL X04BAY(IPRINT,2,REC)
         END IF
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) WRITE (ERRREC,FMT=99991)
C
      ELSE IF (MSG.EQ.'rowerr') THEN
         INFORM = 3
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99981)
            CALL X04BAY(IPRINT,3,REC)
         END IF
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) WRITE (ERRREC,FMT=99980)
C
      ELSE IF (MSG.EQ.'itnlim') THEN
         INFORM = 4
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99990)
            CALL X04BAY(IPRINT,2,REC)
         END IF
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) WRITE (ERRREC,FMT=99989)
C
      ELSE IF (MSG.EQ.'errors') THEN
         INFORM = 6
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99988) NERROR
            CALL X04BAY(IPRINT,2,REC)
         END IF
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) WRITE (ERRREC,FMT=99987)
     *       NERROR
C
      ELSE IF (MSG.EQ.'noprob') THEN
         INFORM = 7
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99986)
            CALL X04BAY(IPRINT,2,REC)
         END IF
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.-1) WRITE (ERRREC,FMT=99985)
      END IF
C
      IF (MSGLVL.GT.0) THEN
C
         IF (INFORM.LT.5) THEN
            IF (NUMINF.EQ.0) THEN
               IF (PRBTYP.NE.'FP') THEN
                  WRITE (REC,FMT=99984) PRBTYP, OBJ
                  CALL X04BAY(IPRINT,2,REC)
               END IF
            ELSE IF (INFORM.EQ.3) THEN
               IF (MSG.EQ.'infeas') THEN
                  IF (MINSUM.EQ.0) THEN
                     WRITE (REC,FMT=99979) OBJ
                     CALL X04BAY(IPRINT,2,REC)
                  ELSE
                     WRITE (REC,FMT=99983) OBJ
                     CALL X04BAY(IPRINT,2,REC)
                  END IF
               ELSE IF (MSG.EQ.'rowerr') THEN
                  WRITE (REC,FMT=99978) ERRMAX
                  CALL X04BAY(IPRINT,2,REC)
               END IF
            ELSE
               WRITE (REC,FMT=99982) OBJ
               CALL X04BAY(IPRINT,2,REC)
            END IF
         END IF
      END IF
C
      IF (INFORM.LT.6) THEN
         IF (MSGLVL.GT.0) THEN
            WRITE (REC,FMT=99999) PRBTYP, ITER
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
C     End of  E04MFF.  (LPOPT)
C
99999 FORMAT (/' Exit from ',A2,' problem after ',I5,' iterations.')
99998 FORMAT (/' Exit E04MFF - Feasible point found.     ')
99997 FORMAT (/' Exit E04MFF - Optimal ',A2,' solution.')
99996 FORMAT (/' Exit E04MFF - Weak ',A2,' solution.')
99995 FORMAT (/' ** Weak ',A2,' solution.')
99994 FORMAT (/' Exit E04MFF - ',A2,' solution is unbounded.')
99993 FORMAT (/' ** ',A2,' solution is unbounded.')
99992 FORMAT (/' Exit E04MFF - No feasible point for the linear constr',
     *       'aints.')
99991 FORMAT (/' ** No feasible point for the linear constraints.')
99990 FORMAT (/' Exit E04MFF - Too many iterations.')
99989 FORMAT (/' ** Too many iterations.')
99988 FORMAT (/' Exit E04MFF - ',I7,' errors found in the input parame',
     *       'ters.  Problem abandoned.')
99987 FORMAT (/' ** ',I7,' errors found in the input parameters.  Prob',
     *       'lem abandoned.')
99986 FORMAT (/' Exit E04MFF - Problem type not recognized.  Problem a',
     *       'bandoned.')
99985 FORMAT (/' ** Problem type not recognized.  Problem abandoned.')
99984 FORMAT (/' Final ',A2,' objective value =',G16.7)
99983 FORMAT (/' Minimum sum of infeasibilities =',G16.7)
99982 FORMAT (/' Final sum of infeasibilities =',G16.7)
99981 FORMAT (/' Exit E04MFF - Cannot satisfy the working set constrai',
     *       'nts to the accuracy',/'requested.')
99980 FORMAT (/' ** Cannot satisfy the working set constraints to the ',
     *       'accuracy requested.')
99979 FORMAT (/' Sum of infeasibilities =',G16.7)
99978 FORMAT (/' Maximum row error =',G16.7)
      END
