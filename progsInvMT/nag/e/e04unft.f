      SUBROUTINE E04UNF(M,N,NCLIN,NCNLN,LDA,LDCJU,LDFJU,LDR,A,BL,BU,Y,
     *                  CONFUN,OBJFUN,ITER,ISTATE,C,CJACU,F,FJACU,
     *                  CLAMDA,OBJF,R,X,IW,LENIW,W,LENW,IUSER,USER,
     *                  IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     ------------------------------------------------------------------
C
C     E04UNF  solves the constrained least-squares problem
C
C            minimize           F(x) = (1/2)(y - f(x)'(y - f(x))
C
C
C                                    (    x  )
C            subject to    bl  .le.  (  A*x  )  .le.  bu
C                                    (  c(x) )
C
C     where  y  is a given constant m-vector, f(x)  is an m-vector of
C     smooth functions,  A  is a constant matrix and  c(x)  is a vector
C     of smooth nonlinear functions. The feasible region is defined
C     by a mixture of linear and nonlinear equality or inequality
C     constraints on  x.
C
C     The dimensions of the problem are...
C
C     M        the number of data points (the dimension of  y),
C
C     N        the number of variables (dimension of  x),
C
C     NCLIN    the number of linear constraints (rows of the matrix  A),
C
C     NCNLN    the number of nonlinear constraints (dimension of  c(x)).
C
C
C     E04UNF  uses a sequential quadratic programming algorithm, with a
C     positive-definite quasi-Newton approximation to the transformed
C     Hessian  Q'HQ  of the Lagrangian function (which will be stored in
C     the array  R).
C
C     Version 4.04, June        9, 1989. Matches NPSOL Version 4.04.
C     Version 4.06, November    5, 1991. Matches NPSOL Version 4.06.
C     Version 5.01, July       21, 1994. Data vector y added.
C                                        Matches NPSOL Version 5.01.
C
C     Copyright  1983--1994  Optimates.
C     This software is not in the public domain. Its use is governed
C     by a license agreement with Stanford University. It is illegal
C     to make copies except as authorised by the license agreement.
C
C     ------------------------------------------------------------------
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E04UNF')
      INTEGER           MXPARM
      PARAMETER         (MXPARM=30)
      INTEGER           LENLS
      PARAMETER         (LENLS=20)
      INTEGER           LENNP
      PARAMETER         (LENNP=35)
      INTEGER           LENNL
      PARAMETER         (LENNL=20)
      DOUBLE PRECISION  ZERO, POINT3, POINT8
      PARAMETER         (ZERO=0.0D+0,POINT3=3.3D-1,POINT8=0.8D+0)
      DOUBLE PRECISION  POINT9, ONE, TEN
      PARAMETER         (POINT9=0.9D+0,ONE=1.0D+0,TEN=10.0D+0)
      DOUBLE PRECISION  HUNDRD, GROWTH
      PARAMETER         (HUNDRD=1.0D+2,GROWTH=1.0D+2)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  OBJF
      INTEGER           IFAIL, ITER, LDA, LDCJU, LDFJU, LDR, LENIW,
     *                  LENW, M, N, NCLIN, NCNLN
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), BL(N+NCLIN+NCNLN), BU(N+NCLIN+NCNLN),
     *                  C(*), CJACU(LDCJU,*), CLAMDA(N+NCLIN+NCNLN),
     *                  F(M), FJACU(LDFJU,*), R(LDR,*), USER(*),
     *                  W(LENW), X(N), Y(M)
      INTEGER           ISTATE(N+NCLIN+NCNLN), IUSER(*), IW(LENIW)
C     .. Subroutine Arguments ..
      EXTERNAL          CONFUN, OBJFUN
C     .. Scalars in Common ..
      DOUBLE PRECISION  ASIZE, BIGBND, BIGDX, BNDLOW, BNDUPP, CDINT,
     *                  CTOL, DRMAX, DRMIN, DTMAX, DTMIN, DXLIM, EPSPT3,
     *                  EPSPT5, EPSPT8, EPSPT9, EPSRF, ETA, FDINT, FTOL,
     *                  HCNDBD, RCNDBD, RFROBN, RHODMP, RHOMAX, RHONRM,
     *                  SCALE, TOLACT, TOLFEA, TOLRNK
      INTEGER           IPRINT, IPRNT, ISUMM, ISUMRY, ITMAX1, ITMAX2,
     *                  ITMXNP, JVRFY1, JVRFY2, JVRFY3, JVRFY4, KSAVE,
     *                  LCRASH, LDQ, LDT, LENNAM, LFDSET, LFORMH,
     *                  LINES1, LINES2, LPROB, LTYPEH, LVERFY, LVLDER,
     *                  LVLDIF, LVRFYC, MSGLS, MSGNP, NACTIV, NCDIFF,
     *                  NCOLT, NFDIFF, NFREE, NLNF, NLNJ, NLNX, NLOAD,
     *                  NN, NNCLIN, NNCNLN, NOUT, NPROB, NRESET, NSAVE,
     *                  NZ
      LOGICAL           INCRUN, UNITQ
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADLS(23), RPADNL(30), RPADNP(22),
     *                  RPSVLS(MXPARM), RPSVNL(MXPARM), RPSVNP(MXPARM),
     *                  WMACH(15)
      INTEGER           IPADLS(19), IPADNL(28), IPADNP(15),
     *                  IPSVLS(MXPARM), IPSVNL(MXPARM), IPSVNP(MXPARM),
     *                  JVERFY(4), LOCLS(LENLS), LOCNL(LENNL),
     *                  LOCNP(LENNP)
C     .. Local Scalars ..
      DOUBLE PRECISION  AMIN, COND, CONDMX, CTX, DXNORM, EPSMCH, ERRMAX,
     *                  FDCHK, FDNORM, FEAMAX, FEAMIN, OBJ, ROOTN,
     *                  RTEPS, SSQ1, SUMINF, XNORM
      INTEGER           I, IANRMJ, IKX, INFO, INFORM, ITMXSV, ITNS, J,
     *                  JINF, JMAX, LANORM, LAQP, LAX, LCJAC, LCJDX,
     *                  LCLAM, LCMUL, LDAQP, LDCJ, LDFJ, LDX, LF,
     *                  LFEATL, LFJAC, LFJDX, LGQ, LGRAD, LHCTRL,
     *                  LHFRWD, LIPERM, LITOTL, LKACTV, LKX, LNEEDC, LQ,
     *                  LRES, LRES0, LRLAM, LT, LWRK1, LWRK2, LWRK3,
     *                  LWTINF, LWTOTL, LX1, LYF, MAXACT, MAXNZ, MINACT,
     *                  MINFXD, MODE, MSGQP, MXFREE, NACT1, NARTIF,
     *                  NCTOTL, NERR, NERROR, NFUN, NGQ, NGRAD, NLPERR,
     *                  NMAJOR, NMINOR, NPLIN, NRANK, NREJTD, NRES,
     *                  NSTATE, NUMINF, NZ1
      LOGICAL           COLD, LINOBJ, NAMED, NEEDFD, OVERFL, ROWERR,
     *                  VERTEX
      CHARACTER*11      TITLE
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMLS(MXPARM), RPRMNL(MXPARM), RPRMNP(MXPARM)
      INTEGER           IPRMLS(MXPARM), IPRMNL(MXPARM), IPRMNP(MXPARM)
      CHARACTER*8       NAMES(1)
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2, F06BLF, F06RJF
      INTEGER           P01ABF
      EXTERNAL          DNRM2, F06BLF, F06RJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DSCAL, E04NBW, E04NBZ, E04NCG,
     *                  E04NCH, E04NCU, E04NCX, E04NCY, E04NCZ, E04UCP,
     *                  E04UCS, E04UDR, E04UNP, E04UNV, E04UNZ, E04UPS,
     *                  E04UPX, E04UPY, F06DFF, F06FBF, F06FLF, F06QFF,
     *                  F06QHF, X02ZAZ, X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX, MIN, SQRT
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AE04NC/LOCLS
      COMMON            /AE04UC/LOCNP
      COMMON            /AE04UP/LOCNL
      COMMON            /AX02ZA/WMACH
      COMMON            /BE04NB/LENNAM, LDT, NCOLT, LDQ
      COMMON            /BE04UC/LVLDIF, NCDIFF, NFDIFF, LFDSET
      COMMON            /BE04UP/IPSVNL, LTYPEH, NRESET, IPADNL
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /CE04UC/LVRFYC, JVERFY
      COMMON            /CE04UP/RPSVNL, RPADNL
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
      COMMON            /DE04NC/IPSVLS, IPRNT, ISUMRY, ITMAX1, ITMAX2,
     *                  LCRASH, LFORMH, LPROB, MSGLS, NN, NNCLIN, NPROB,
     *                  IPADLS
      COMMON            /DE04UC/RHOMAX, RHONRM, RHODMP, SCALE, INCRUN
      COMMON            /EE04NB/RCNDBD, RFROBN, DRMAX, DRMIN
      COMMON            /EE04NC/RPSVLS, BIGBND, BIGDX, BNDLOW, BNDUPP,
     *                  TOLACT, TOLFEA, TOLRNK, RPADLS
      COMMON            /FE04NC/NACTIV, NFREE, NZ, UNITQ
      COMMON            /GE04UC/IPSVNP, ITMXNP, JVRFY1, JVRFY2, JVRFY3,
     *                  JVRFY4, LVLDER, LVERFY, MSGNP, NLNF, NLNJ, NLNX,
     *                  NNCNLN, NSAVE, NLOAD, KSAVE, IPADNP
      COMMON            /HE04UC/RPSVNP, CDINT, CTOL, DXLIM, EPSRF, ETA,
     *                  FDINT, FTOL, HCNDBD, RPADNP
C     .. Equivalences ..
      EQUIVALENCE       (IPRMLS(1),IPRNT), (RPRMLS(1),BIGBND)
      EQUIVALENCE       (IPRMNP(1),ITMXNP), (RPRMNP(1),CDINT)
      EQUIVALENCE       (IPRMNL(1),LTYPEH), (RPRMNL(1),RPADNL(1))
      EQUIVALENCE       (ITMXNP,NMAJOR), (ITMAX2,NMINOR), (MSGLS,MSGQP)
C     .. Save statement ..
      SAVE              /DE04NC/, /EE04NC/, /GE04UC/, /HE04UC/,
     *                  /BE04UP/, /CE04UP/, /AX02ZA/, /FE04NC/
C     .. Data statements ..
      DATA              TITLE/' *** E04UNF'/
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
      RHOMAX = ONE/EPSMCH
      ROOTN = SQRT(DBLE(N))
C
C     Default names will be provided for variables during printing.
C
      NAMED = .FALSE.
      INFORM = 0
      ITER = 0
C
C     Set the default values for the parameters.
C
      CALL E04UPX(M,N,NCLIN,NCNLN,TITLE)
C
      NEEDFD = LVLDER .EQ. 0 .OR. LVLDER .EQ. 2 .OR.
     *         (LVLDER.EQ.1 .AND. NCNLN.GT.0)
      COLD = LCRASH .EQ. 0
      LVLDIF = 0
      IF (NEEDFD) LVLDIF = 1
C
      NPLIN = N + NCLIN
      NCTOTL = NPLIN + NCNLN
C
C     Assign the dimensions of arrays in the parameter list of E04UNZ.
C     Economies of storage are possible if the minimum number of active
C     constraints and the minimum number of fixed variables are known in
C     advance.  The expert user should alter MINACT and MINFXD
C     accordingly.
C
      MINACT = 0
      MINFXD = 0
C
      MXFREE = N - MINFXD
      MAXACT = MAX(1,MIN(N,NCLIN))
      MAXNZ = N - (MINFXD+MINACT)
C
      IF (NCLIN+NCNLN.EQ.0) THEN
         LDQ = 1
         LDT = 1
         NCOLT = 1
      ELSE
         LDQ = MAX(1,MXFREE)
         LDT = MAX(MAXNZ,MAXACT)
         NCOLT = MXFREE
      END IF
C
      LENNAM = 1
C
      LDAQP = MAX(NCLIN+NCNLN,1)
C
C     E04UNV defines the arrays that contain the locations of various
C     work arrays within  W  and  IW.
C
      LITOTL = 0
      LWTOTL = 0
      CALL E04UCP(N,NCLIN,NCNLN,NCTOTL,LITOTL,LWTOTL)
      CALL E04UNV(M,N,LITOTL,LWTOTL)
C
C     Allocate addresses not allocated in E04UCP or E04UNV.
C
      LAX = LWTOTL + 1
      LWTOTL = LAX + NCLIN - 1
      LAX = MIN(LAX,LWTOTL)
C
C     Check input parameters and storage limits.
C
      CALL E04NBZ(NERROR,MSGNP,LCRASH,LENIW,LENW,LITOTL,LWTOTL,N,NCLIN,
     *            NCNLN,ISTATE,NAMED,NAMES,BIGBND,BL,BU,CLAMDA,M,LDA,
     *            LDR,LDCJU,LDFJU,NERR,IFAIL)
C
      IF (NERROR.GT.0) THEN
         INFORM = 9
         GO TO 60
      END IF
C
      LKACTV = LOCLS(1)
      LANORM = LOCLS(2)
      LCJDX = LOCLS(3)
      LRES = LOCLS(5)
      LRES0 = LOCLS(6)
      LGQ = LOCLS(9)
      LRLAM = LOCLS(10)
      LT = LOCLS(11)
      LQ = LOCLS(12)
      LWTINF = LOCLS(13)
      LWRK1 = LOCLS(14)
C
      LKX = LOCNP(1)
      LIPERM = LOCNP(2)
      LAQP = LOCNP(3)
      LDX = LOCNP(7)
      LFEATL = LOCNP(10)
      LX1 = LOCNP(11)
      LWRK2 = LOCNP(12)
C
      LCMUL = LOCNP(16)
      LWRK3 = LOCNP(21)
      LNEEDC = LOCNP(24)
      LHFRWD = LOCNP(25)
      LHCTRL = LOCNP(26)
      LCJAC = LOCNP(27)
      LGRAD = LOCNP(28)
C
      LF = LOCNL(1)
      LFJAC = LOCNL(2)
      LFJDX = LOCNL(3)
      LYF = LOCNL(4)
C
      LDCJ = MAX(NCNLN,1)
      LDFJ = MAX(M,1)
C
      TOLRNK = ONE/HCNDBD
      RCNDBD = SQRT(HCNDBD)
C
C     ------------------------------------------------------------------
C     Load the arrays of feasibility tolerances.
C     ------------------------------------------------------------------
      IF (TOLFEA.GT.ZERO) CALL F06FBF(NPLIN,TOLFEA,W(LFEATL),1)
C
      IF (NCNLN.GT.0 .AND. CTOL.GT.ZERO) CALL F06FBF(NCNLN,CTOL,
     *    W(LFEATL+NPLIN),1)
C
      IF (LFDSET.EQ.0) THEN
         FDCHK = SQRT(EPSRF)
      ELSE IF (LFDSET.EQ.1) THEN
         FDCHK = FDINT
      ELSE
         FDCHK = W(LHFRWD)
      END IF
C
      NFUN = 0
      NGRAD = 0
      NSTATE = 1
C
      XNORM = DNRM2(N,X,1)
      CALL DCOPY(N,X,1,W(LX1),1)
C
C     ------------------------------------------------------------------
C     If required,  compute the problem functions.
C     If the constraints are nonlinear,  the first call of CONFUN
C     sets up any constant elements in the Jacobian matrix.  A copy of
C     the Jacobian (with constant elements set) is placed in  CJACU.
C     ------------------------------------------------------------------
      IF (LVERFY.GE.10) THEN
         LVRFYC = LVERFY - 10
C
         CALL E04UPY(INFO,MSGNP,NSTATE,LVLDER,NFUN,NGRAD,LDCJ,LDCJU,
     *               LDFJ,LDFJU,M,N,NCNLN,CONFUN,OBJFUN,IW(LNEEDC),
     *               BIGBND,EPSRF,CDINT,FDINT,FDCHK,FDNORM,XNORM,BL,BU,
     *               C,W(LWRK3),W(LCJAC),CJACU,W(LCJDX),F,W(LF),W(LFJAC)
     *               ,FJACU,W(LFJDX),W(LDX),W(LHFRWD),W(LHCTRL),X,
     *               W(LWRK1),W(LWRK2),W(LYF),IUSER,USER)
C
         IF (INFO.NE.0) THEN
            IF (INFO.GT.0) INFORM = 7
            IF (INFO.LT.0) INFORM = INFO
            GO TO 60
         END IF
         NSTATE = 0
      END IF
C
      IF (LCRASH.LT.2) THEN
C        ===============================================================
C        Cold or warm start.  Use  E04NCZ  to obtain a point that
C        satisfies the linear constraints.
C        ===============================================================
C
         IF (NCLIN.GT.0) THEN
            IANRMJ = LANORM
            DO 20 J = 1, NCLIN
               W(IANRMJ) = DNRM2(N,A(J,1),LDA)
               IANRMJ = IANRMJ + 1
   20       CONTINUE
            CALL F06FLF(NCLIN,W(LANORM),1,ASIZE,AMIN)
         END IF
C
         CALL F06FLF(NPLIN,W(LFEATL),1,FEAMAX,FEAMIN)
         CALL DCOPY(NPLIN,W(LFEATL),1,W(LWTINF),1)
         CALL DSCAL(NPLIN,(ONE/FEAMIN),W(LWTINF),1)
C
C        ===============================================================
C        The input values of X and (optionally)  ISTATE are used by
C        E04NCU  to define an initial working set.
C        ===============================================================
         VERTEX = .FALSE.
         CALL E04NCU(COLD,VERTEX,NCLIN,NPLIN,NACTIV,NARTIF,NFREE,N,LDA,
     *               ISTATE,IW(LKACTV),BIGBND,TOLACT,A,W(LAX),BL,BU,X,
     *               W(LWRK1),W(LWRK2))
C
         UNITQ = .TRUE.
         NRES = 0
         NGQ = 0
         CONDMX = MAX(ONE/EPSPT5,HUNDRD)
C
         IKX = LKX
         DO 40 I = 1, N
            IW(IKX) = I
            IKX = IKX + 1
   40    CONTINUE
C
         IF (COLD) THEN
            NRANK = 0
         ELSE
            NRANK = NLNX
            CALL F06FBF(NLNX,(ZERO),W(LRES0),1)
         END IF
C
C        ---------------------------------------------------------------
C        Re-order KX so that the free variables come first.
C        If a warm start is required, NRANK will be nonzero and the
C        factor R will be updated.
C        ---------------------------------------------------------------
         CALL E04NCX(UNITQ,INFORM,NZ,NFREE,NRANK,NRES,NGQ,N,LDQ,LDA,LDR,
     *               LDT,ISTATE,IW(LKX),CONDMX,A,R,W(LT),W(LRES0),W(LGQ)
     *               ,W(LQ),W(LWRK1),W(LWRK2),W(LRLAM),MSGNP)
C
C        ---------------------------------------------------------------
C        Factorize the linear constraints in the initial working set.
C        ---------------------------------------------------------------
         IF (NACTIV.GT.0) THEN
            NACT1 = NACTIV
            NACTIV = 0
C
            CALL E04NCY(UNITQ,VERTEX,INFORM,1,NACT1,NACTIV,NARTIF,NZ,
     *                  NFREE,NRANK,NREJTD,NRES,NGQ,N,LDQ,LDA,LDR,LDT,
     *                  ISTATE,IW(LKACTV),IW(LKX),CONDMX,A,R,W(LT),
     *                  W(LRES0),W(LGQ),W(LQ),W(LWRK1),W(LWRK2),W(LRLAM)
     *                  ,MSGNP)
         END IF
C
         SSQ1 = ZERO
         LINOBJ = .FALSE.
         CALL E04NCH(LINOBJ,ROWERR,UNITQ,NCLIN,NACTIV,NFREE,NRANK,NZ,N,
     *               NPLIN,LDQ,LDA,LDR,LDT,ISTATE,IW(LKACTV),IW(LKX),
     *               JMAX,ERRMAX,CTX,XNORM,A,W(LAX),BL,BU,W(LGQ),W(LRES)
     *               ,W(LRES0),W(LFEATL),R,W(LT),X,W(LQ),W(LWRK1),
     *               W(LWRK2))
C
C        ---------------------------------------------------------------
C        Call  E04NCZ  to find a feasible  x.
C        ---------------------------------------------------------------
C        Use  WORK2  as the multiplier vector.
C
         JINF = 0
         LCLAM = LWRK2
C
         ITMXSV = ITMAX1
         ITMAX1 = NMINOR
C
         CALL E04NCZ('FP problem',NAMED,NAMES,LINOBJ,UNITQ,NLPERR,ITNS,
     *               JINF,NCLIN,NPLIN,NACTIV,NFREE,NRANK,NZ,NZ1,N,LDA,
     *               LDR,ISTATE,IW(LKACTV),IW(LKX),CTX,OBJ,SSQ1,SUMINF,
     *               NUMINF,XNORM,BL,BU,A,W(LCLAM),W(LAX),W(LFEATL),R,X,
     *               W)
C
         ITMAX1 = ITMXSV
C
         IF (NLPERR.GT.0) THEN
            INFORM = 2
            GO TO 60
         ELSE IF (MSGQP.GT.0) THEN
            WRITE (REC,FMT=99987)
            CALL X04BAY(IPRINT,2,REC)
         END IF
      END IF
C
C     ==================================================================
C     Check the gradients at this first feasible x.
C     ==================================================================
      CALL DAXPY(N,(-ONE),X,1,W(LX1),1)
      DXNORM = DNRM2(N,W(LX1),1)
C
      IF (LVERFY.GE.10 .AND. DXNORM.LE.TEN*EPSMCH) THEN
C        Relax, we already have everything at this x.
      ELSE
         LVRFYC = LVERFY
         IF (LVERFY.GE.10) LVRFYC = -1
C
         CALL E04UPY(INFO,MSGNP,NSTATE,LVLDER,NFUN,NGRAD,LDCJ,LDCJU,
     *               LDFJ,LDFJU,M,N,NCNLN,CONFUN,OBJFUN,IW(LNEEDC),
     *               BIGBND,EPSRF,CDINT,FDINT,FDCHK,FDNORM,XNORM,BL,BU,
     *               C,W(LWRK3),W(LCJAC),CJACU,W(LCJDX),F,W(LF),W(LFJAC)
     *               ,FJACU,W(LFJDX),W(LDX),W(LHFRWD),W(LHCTRL),X,
     *               W(LWRK1),W(LWRK2),W(LYF),IUSER,USER)
C
         IF (INFO.NE.0) THEN
            IF (INFO.GT.0) INFORM = 7
            IF (INFO.LT.0) INFORM = INFO
            GO TO 60
         END IF
      END IF
C
      MODE = 2
      CALL E04UNP(MODE,N,M,Y,F,W(LFJAC),LDFJ,OBJF,W(LGQ),W(LYF))
      CALL E04NBW(6,N,NZ,NFREE,LDQ,UNITQ,IW(LKX),W(LGQ),W(LQ),W(LWRK1))
C
      IF (COLD) THEN
C        ---------------------------------------------------------------
C        Cold start.
C        ---------------------------------------------------------------
         IF (LTYPEH.EQ.0) THEN
C
C           Initialize  R'R  as J'J.
C
            CALL E04UPS(N,LDQ,NFREE,UNITQ,IW(LKX),M,W(LFJAC),LDFJ,R,LDR,
     *                  W(LQ),W(LWRK1))
            RFROBN = F06RJF('Frobenius norm','Upper',
     *               'Non-unit diagonal',N,N,R,LDR,W)
         ELSE
C
C           Initialize  R'R  as the identity.
C
            CALL F06QHF('Upper-triangular',N,N,ZERO,ONE,R,LDR)
            RFROBN = ROOTN
         END IF
         IF (NCNLN.GT.0) CALL F06FBF(NCNLN,(ZERO),W(LCMUL),1)
      ELSE
C        ---------------------------------------------------------------
C        Warm start.
C        Set the multipliers for the nonlinear constraints.
C        Check the condition of the initial factor R.
C        ---------------------------------------------------------------
         IF (NCNLN.GT.0) CALL DCOPY(NCNLN,CLAMDA(NPLIN+1),1,W(LCMUL),1)
         RFROBN = F06RJF('Frobenius norm','Upper','Non-unit diagonal',N,
     *            N,R,LDR,W)
      END IF
C
C     ------------------------------------------------------------------
C     Check the condition of the initial factor R.
C     Compute the diagonal condition estimator of R.
C     ------------------------------------------------------------------
      CALL F06FLF(N,R,LDR+1,DRMAX,DRMIN)
      COND = F06BLF(DRMAX,DRMIN,OVERFL)
C
      IF (COND.GT.RCNDBD .OR. RFROBN.GT.ROOTN*GROWTH*DRMAX) THEN
C        ---------------------------------------------------------------
C        Refactorize the Hessian and bound the condition estimator.
C        ---------------------------------------------------------------
         IF ( .NOT. COLD) THEN
            IF (MSGNP.GT.0) THEN
               WRITE (REC,FMT=99986)
               CALL X04BAF(IPRINT,REC(1))
            END IF
         END IF
         CALL E04UDR(UNITQ,N,NFREE,NZ,LDQ,LDR,IW(LIPERM),IW(LKX),W(LGQ),
     *               R,W(LQ),W(LWRK1),W(LRES0))
      END IF
C
C     ==================================================================
C     Solve the problem.
C     ==================================================================
      IF (NCNLN.EQ.0) THEN
C        ---------------------------------------------------------------
C        The problem has only linear constraints and bounds.
C        ---------------------------------------------------------------
         CALL E04UNZ(NAMED,NAMES,UNITQ,INFORM,ITER,M,N,NCLIN,NCNLN,
     *               NACTIV,NFREE,NZ,LDCJ,LDCJU,LDFJ,LDFJU,LDA,LDR,NFUN,
     *               NGRAD,ISTATE,IW(LKACTV),IW(LKX),OBJF,FDNORM,XNORM,
     *               CONFUN,OBJFUN,A,W(LAX),BL,BU,C,W(LCJAC),CJACU,
     *               CLAMDA,Y,F,W(LFJAC),FJACU,W(LFEATL),W(LGRAD),R,X,
     *               IW,W,LENW,IUSER,USER)
      ELSE
C        ---------------------------------------------------------------
C        The problem has some nonlinear constraints.
C        ---------------------------------------------------------------
         IF (NCLIN.GT.0) CALL F06QFF('General',NCLIN,N,A,LDA,W(LAQP),
     *                               LDAQP)
C
C        Try and add some nonlinear constraint indices to KACTIV.
C
         CALL E04UCS(COLD,N,NCLIN,NCNLN,NCTOTL,NACTIV,NFREE,NZ,ISTATE,
     *               IW(LKACTV),BIGBND,TOLACT,BL,BU,C)
C
         CALL E04UNZ(NAMED,NAMES,UNITQ,INFORM,ITER,M,N,NCLIN,NCNLN,
     *               NACTIV,NFREE,NZ,LDCJ,LDCJU,LDFJ,LDFJU,LDAQP,LDR,
     *               NFUN,NGRAD,ISTATE,IW(LKACTV),IW(LKX),OBJF,FDNORM,
     *               XNORM,CONFUN,OBJFUN,W(LAQP),W(LAX),BL,BU,C,W(LCJAC)
     *               ,CJACU,CLAMDA,Y,F,W(LFJAC),FJACU,W(LFEATL),W(LGRAD)
     *               ,R,X,IW,W,LENW,IUSER,USER)
C
      END IF
C
C     ------------------------------------------------------------------
C     If required, overwrite R with the factor of the Hessian.
C     ------------------------------------------------------------------
C
      IF (LFORMH.GT.0) THEN
         CALL E04NCG('Hessian',UNITQ,NFREE,N,N,LDQ,LDR,IW(LKX),R,W(LQ),
     *               W(LWRK1),W(LWRK2))
      END IF
C
C     Print messages if required.
C
   60 IF (MSGNP.GT.0) THEN
         IF (INFORM.LT.0) WRITE (REC,FMT=99999)
         IF (INFORM.EQ.0) WRITE (REC,FMT=99998)
         IF (INFORM.EQ.1) WRITE (REC,FMT=99997)
         IF (INFORM.EQ.2) WRITE (REC,FMT=99996)
         IF (INFORM.EQ.3) WRITE (REC,FMT=99995)
         IF (INFORM.EQ.4) WRITE (REC,FMT=99994)
         IF (INFORM.EQ.6) WRITE (REC,FMT=99993)
         IF (INFORM.EQ.7) WRITE (REC,FMT=99992)
         IF (INFORM.EQ.9) WRITE (REC,FMT=99991) NERROR
         CALL X04BAY(IPRINT,2,REC)
C
         IF (INFORM.GE.0 .AND. INFORM.LT.7) THEN
            IF (NLPERR.EQ.0) THEN
               WRITE (REC,FMT=99990) OBJF
               CALL X04BAY(IPRINT,2,REC)
            ELSE
               IF (NLPERR.EQ.3) THEN
                  WRITE (REC,FMT=99989) SUMINF
                  CALL X04BAY(IPRINT,2,REC)
               ELSE
                  WRITE (REC,FMT=99988) SUMINF
                  CALL X04BAY(IPRINT,2,REC)
               END IF
            END IF
         END IF
      END IF
C
C     Recover the optional parameters set by the user.
C
      CALL F06DFF(MXPARM,IPSVLS,1,IPRMLS,1)
      CALL DCOPY(MXPARM,RPSVLS,1,RPRMLS,1)
      CALL F06DFF(MXPARM,IPSVNP,1,IPRMNP,1)
      CALL DCOPY(MXPARM,RPSVNP,1,RPRMNP,1)
      CALL F06DFF(MXPARM,IPSVNL,1,IPRMNL,1)
      CALL DCOPY(MXPARM,RPSVNL,1,RPRMNL,1)
C
      IF (INFORM.LT.9) THEN
         IF (NCNLN.GT.0) CALL F06QFF('General',NCNLN,N,W(LCJAC),LDCJ,
     *                               CJACU,LDCJU)
         CALL F06QFF('General',M,N,W(LFJAC),LDFJ,FJACU,LDFJU)
      END IF
C
      IF (INFORM.NE.0 .AND. (IFAIL.EQ.0 .OR. IFAIL.EQ.-1)) THEN
         IF (INFORM.LT.0) WRITE (REC,FMT=99985)
         IF (INFORM.EQ.1) WRITE (REC,FMT=99984)
         IF (INFORM.EQ.2) WRITE (REC,FMT=99983)
         IF (INFORM.EQ.3) WRITE (REC,FMT=99982)
         IF (INFORM.EQ.4) WRITE (REC,FMT=99981)
         IF (INFORM.EQ.6) WRITE (REC,FMT=99980)
         IF (INFORM.EQ.7) WRITE (REC,FMT=99979)
         IF (INFORM.EQ.9) WRITE (REC,FMT=99978) NERROR
         CALL X04BAY(NERR,2,REC)
      END IF
      IFAIL = P01ABF(IFAIL,INFORM,SRNAME,0,REC)
      RETURN
C
C
C
C     End of  E04UNF.  (NLSSOL)
C
99999 FORMAT (/' Exit E04UNF - User requested termination.')
99998 FORMAT (/' Exit E04UNF - Optimal solution found.')
99997 FORMAT (/' Exit E04UNF - Optimal solution found, but requested a',
     *       'ccuracy not achieved.')
99996 FORMAT (/' Exit E04UNF - No feasible point for the linear constr',
     *       'aints.')
99995 FORMAT (/' Exit E04UNF - No feasible point for the nonlinear con',
     *       'straints.')
99994 FORMAT (/' Exit E04UNF - Too many major iterations.             ')
99993 FORMAT (/' Exit E04UNF - Current point cannot be improved upon. ')
99992 FORMAT (/' Exit E04UNF - Large errors found in the derivatives. ')
99991 FORMAT (/' Exit E04UNF - ',I7,' errors found in the input parame',
     *       'ters.  Problem abandoned.')
99990 FORMAT (/' Final objective value =',G16.7)
99989 FORMAT (/' Minimum sum of infeasibilities =',G16.7)
99988 FORMAT (/' Final sum of infeasibilities =',G16.7)
99987 FORMAT (/' The linear constraints are feasible.')
99986 FORMAT (' XXX  Bad initial Hessian,   R  refactorized.')
99985 FORMAT (/' ** User requested termination.')
99984 FORMAT (/' ** Optimal solution found, but requested accuracy not',
     *       ' achieved.')
99983 FORMAT (/' ** No feasible point for the linear constraints.')
99982 FORMAT (/' ** No feasible point for the nonlinear constraints.')
99981 FORMAT (/' ** Too many major iterations.             ')
99980 FORMAT (/' ** Current point cannot be improved upon. ')
99979 FORMAT (/' ** Large errors found in the derivatives. ')
99978 FORMAT (/' ** ',I7,' errors found in the input parameters.  Prob',
     *       'lem abandoned.')
      END
