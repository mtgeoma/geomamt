      SUBROUTINE E04UCZ(NAMED,NAMES,UNITQ,INFORM,MAJITS,N,NCLIN,NCNLN,
     *                  NCTOTL,NACTIV,NFREE,NZ,LDCJ,LDCJU,LDAQP,LDR,
     *                  NFUN,NGRAD,ISTATE,KACTIV,KX,OBJF,FDNORM,XNORM,
     *                  CONFUN,OBJFUN,AQP,AX,BL,BU,C,CJAC,CJACU,CLAMDA,
     *                  FEATOL,GRAD,GRADU,R,X,IW,W,LENW,IUSER,USER)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1093 (JUL 1993).
C     MARK 17 REVISED. IER-1613 (JUN 1995).
C
C     ******************************************************************
C     E04UCZ  is the core routine for  E04UCF,  a sequential quadratic
C     programming (SQP) method for nonlinearly constrained optimization.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version      February-1982.
C     This version of E04UCZ dated 23-Dec-92.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LENLS
      PARAMETER         (LENLS=20)
      INTEGER           LENNP
      PARAMETER         (LENNP=35)
      INTEGER           MXPARM
      PARAMETER         (MXPARM=30)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      DOUBLE PRECISION  GROWTH
      PARAMETER         (GROWTH=1.0D+2)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  FDNORM, OBJF, XNORM
      INTEGER           INFORM, LDAQP, LDCJ, LDCJU, LDR, LENW, MAJITS,
     *                  N, NACTIV, NCLIN, NCNLN, NCTOTL, NFREE, NFUN,
     *                  NGRAD, NZ
      LOGICAL           NAMED, UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  AQP(LDAQP,*), AX(*), BL(NCTOTL), BU(NCTOTL),
     *                  C(*), CJAC(LDCJ,*), CJACU(LDCJU,*),
     *                  CLAMDA(NCTOTL), FEATOL(NCTOTL), GRAD(N),
     *                  GRADU(N), R(LDR,*), USER(*), W(LENW), X(N)
      INTEGER           ISTATE(*), IUSER(*), IW(*), KACTIV(N), KX(N)
      CHARACTER*8       NAMES(*)
C     .. Subroutine Arguments ..
      EXTERNAL          CONFUN, OBJFUN
C     .. Scalars in Common ..
      DOUBLE PRECISION  ASIZE, BIGBND, BIGDX, BNDLOW, BNDUPP, CDINT,
     *                  CTOL, DRMAX, DRMIN, DTMAX, DTMIN, DXLIM, EPSPT3,
     *                  EPSPT5, EPSPT8, EPSPT9, EPSRF, ETA, FDINT, FTOL,
     *                  HCNDBD, RCNDBD, RFROBN, TOLACT, TOLFEA, TOLRNK
      INTEGER           IPRINT, IPRNT, ISUMM, ISUMRY, ITMAX1, ITMAX2,
     *                  ITMXNP, JVRFY1, JVRFY2, JVRFY3, JVRFY4, KSAVE,
     *                  LCRASH, LDT, LDZY, LENNAM, LFDSET, LFORMH,
     *                  LINES1, LINES2, LPROB, LVERFY, LVLDER, LVLDIF,
     *                  MSGLS, MSGNP, NCDIFF, NCOLT, NFDIFF, NLNF, NLNJ,
     *                  NLNX, NLOAD, NN, NNCLIN, NNCNLN, NOUT, NPROB,
     *                  NSAVE
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADLS(23), RPADNP(22), RPSVLS(MXPARM),
     *                  RPSVNP(MXPARM), WMACH(15)
      INTEGER           IPADLS(19), IPADNP(15), IPSVLS(MXPARM),
     *                  IPSVNP(MXPARM), LOCLS(LENLS), LOCNP(LENNP)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALFA, ALFBND, ALFDX, ALFLIM, ALFMAX, ALFMIN,
     *                  ALFSML, CNORM, COND, CONDH, CONDHZ, CONDT,
     *                  CVNORM, DINKY, DRZMAX, DRZMIN, DXNORM, ERRMAX,
     *                  FLMAX, GDX, GFNORM, GLF1, GLF2, GLNORM, GLTEST,
     *                  GRDALF, GTEST, GZNORM, OBJ, OBJALF, OBJSIZ,
     *                  QPCURV, ROOTN, RTFTOL, RTMAX, XSIZE
      INTEGER           INFO, JMAX, LADX, LANORM, LAQP, LBL, LBU,
     *                  LC1MUL, LCJAC1, LCJDX, LCJDX1, LCMUL, LCS1,
     *                  LCS2, LDCJ1, LDLAM, LDSLK, LDX, LGQ, LGQ1,
     *                  LHCTRL, LHFRWD, LHPQ, LINACT, LIPERM, LNEEDC,
     *                  LQPTOL, LQRWRK, LRHO, LRLAM, LRPQ, LSLK, LSLK1,
     *                  LT, LVIOLN, LWRK1, LWRK2, LWRK3, LWTINF, LX1,
     *                  LZY, MAJIT0, MINITS, MNR, MNRSUM, MODE, MSGQP,
     *                  NCQP, NL, NLNACT, NLSERR, NMAJOR, NMINOR, NPLIN,
     *                  NQPERR, NQPINF, NSTATE, NUMINF, NVIOL
      LOGICAL           CENTRL, CONVPT, CONVRG, DONE, ERROR, FEASQP,
     *                  GOODGQ, INFEAS, NEEDFD, NEWGQ, OPTIML, OVERFL
      CHARACTER*5       MJRMSG
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMLS(MXPARM), RPRMNP(MXPARM)
      INTEGER           IPRMLS(MXPARM), IPRMNP(MXPARM)
      LOGICAL           KTCOND(2)
      CHARACTER*80      REC(5)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2, F06BLF
      EXTERNAL          DDOT, DNRM2, F06BLF
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, E04MFK, E04NBW, E04UCL,
     *                  E04UCN, E04UCR, E04UCT, E04UCU, E04UCW, E04UDR,
     *                  E04UDS, E04UDT, E04UPP, F06DBF, F06FLF, F06QFF,
     *                  X04BAF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX, MIN, SQRT
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AE04NC/LOCLS
      COMMON            /AE04UC/LOCNP
      COMMON            /AX02ZA/WMACH
      COMMON            /BE04NB/LENNAM, LDT, NCOLT, LDZY
      COMMON            /BE04UC/LVLDIF, NCDIFF, NFDIFF, LFDSET
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
      COMMON            /DE04NC/IPSVLS, IPRNT, ISUMRY, ITMAX1, ITMAX2,
     *                  LCRASH, LFORMH, LPROB, MSGLS, NN, NNCLIN, NPROB,
     *                  IPADLS
      COMMON            /EE04NB/RCNDBD, RFROBN, DRMAX, DRMIN
      COMMON            /EE04NC/RPSVLS, BIGBND, BIGDX, BNDLOW, BNDUPP,
     *                  TOLACT, TOLFEA, TOLRNK, RPADLS
      COMMON            /GE04UC/IPSVNP, ITMXNP, JVRFY1, JVRFY2, JVRFY3,
     *                  JVRFY4, LVLDER, LVERFY, MSGNP, NLNF, NLNJ, NLNX,
     *                  NNCNLN, NSAVE, NLOAD, KSAVE, IPADNP
      COMMON            /HE04UC/RPSVNP, CDINT, CTOL, DXLIM, EPSRF, ETA,
     *                  FDINT, FTOL, HCNDBD, RPADNP
C     .. Equivalences ..
      EQUIVALENCE       (IPRMLS(1),IPRNT), (RPRMLS(1),BIGBND)
      EQUIVALENCE       (IPRMNP(1),ITMXNP), (RPRMNP(1),CDINT)
      EQUIVALENCE       (ITMXNP,NMAJOR), (ITMAX2,NMINOR), (MSGLS,MSGQP)
C     .. Save statement ..
      SAVE              /AX02ZA/, /DE04NC/, /EE04NC/, /GE04UC/, /HE04UC/
C     .. Executable Statements ..
C
C     Specify machine-dependent parameters.
C
      FLMAX = WMACH(7)
      RTMAX = WMACH(8)
C
      LANORM = LOCLS(2)
      LRPQ = LOCLS(5)
      LQRWRK = LOCLS(6)
      LHPQ = LOCLS(8)
      LGQ = LOCLS(9)
      LRLAM = LOCLS(10)
      LT = LOCLS(11)
      LZY = LOCLS(12)
      LWTINF = LOCLS(13)
      LWRK1 = LOCLS(14)
      LQPTOL = LOCLS(15)
C
      LIPERM = LOCNP(2)
      LAQP = LOCNP(3)
      LADX = LOCNP(4)
      LBL = LOCNP(5)
      LBU = LOCNP(6)
      LDX = LOCNP(7)
      LGQ1 = LOCNP(8)
      LX1 = LOCNP(11)
      LWRK2 = LOCNP(12)
      LCS1 = LOCNP(13)
      LCS2 = LOCNP(14)
      LC1MUL = LOCNP(15)
      LCMUL = LOCNP(16)
      LCJDX1 = LOCNP(17)
      LDLAM = LOCNP(18)
      LDSLK = LOCNP(19)
      LRHO = LOCNP(20)
      LWRK3 = LOCNP(21)
      LSLK1 = LOCNP(22)
      LSLK = LOCNP(23)
      LNEEDC = LOCNP(24)
      LHFRWD = LOCNP(25)
      LHCTRL = LOCNP(26)
C
      LCJAC1 = LAQP + NCLIN
      LCJDX = LADX + NCLIN
      LVIOLN = LWRK3
C
C     Initialize
C
      MJRMSG = '     '
      NQPINF = 0
      MNRSUM = 0
C
      MAJIT0 = MAJITS
      NPLIN = N + NCLIN
      NCQP = NCLIN + NCNLN
      NL = MIN(NPLIN+1,NCTOTL)
C
      LDCJ1 = MAX(NCQP,1)
C
      NEEDFD = LVLDER .EQ. 0 .OR. LVLDER .EQ. 2 .OR.
     *         (LVLDER.EQ.1 .AND. NCNLN.GT.0)
C
      ALFA = ZERO
      ALFDX = ZERO
      RTFTOL = SQRT(FTOL)
      ROOTN = SQRT(DBLE(N))
C
C     ------------------------------------------------------------------
C     Information from the feasibility phase will be used to generate a
C     hot start for the first QP subproblem.
C     ------------------------------------------------------------------
      CALL DCOPY(NCTOTL,FEATOL,1,W(LQPTOL),1)
C
      NSTATE = 0
C
      OBJALF = OBJF
      IF (NCNLN.GT.0) THEN
         OBJALF = OBJALF - DDOT(NCNLN,W(LCMUL),1,C,1)
      END IF
C
      NEWGQ = .FALSE.
C
C*    ==================================================================
C+    repeat                             (until converged or error exit)
C
C     ===============================================================
C     See if we want to save the details of this iteration.
C     ===============================================================
C     20 IF (MOD(MAJITS,KSAVE).EQ.0 .AND. MAJITS.NE.MAJIT0) THEN
C        CALL NPSAVR(UNITQ,N,NCLIN,NCNLN,LDR,LDQ,NFREE,NSAVE,MAJITS,
C    *               ISTATE,KX,W(LHFRWD),W(LHCTRL),W(LCMUL),R,W(LRHO),X,
C    *               X,W(LQ))
C     END IF
C
C   *    ===============================================================
C   +    repeat                         (Until a good gradient is found)
C
   20 MINITS = 0
C
   40 CENTRL = LVLDIF .EQ. 2
C
      IF (NEWGQ) THEN
         IF (NEEDFD) THEN
C           ------------------------------------------------------
C           Compute any missing gradient elements and the
C           transformed gradient of the objective.
C           ------------------------------------------------------
            CALL E04UDS(CENTRL,MODE,LDCJ,LDCJU,N,NCNLN,BIGBND,CDINT,
     *                  FDINT,FDNORM,OBJF,CONFUN,OBJFUN,IW(LNEEDC),BL,
     *                  BU,C,W(LWRK2),W(LWRK3),CJAC,CJACU,GRAD,GRADU,
     *                  W(LHFRWD),W(LHCTRL),X,IUSER,USER)
            INFORM = MODE
            IF (MODE.LT.0) GO TO 60
C
         END IF
C
         CALL DCOPY(N,GRAD,1,W(LGQ),1)
         CALL E04NBW(6,N,NZ,NFREE,LDZY,UNITQ,KX,W(LGQ),W(LZY),W(LWRK1))
         NEWGQ = .FALSE.
      END IF
C
C     ============================================================
C     (1) Solve an inequality quadratic program (IQP) for the
C         search direction and multiplier estimates.
C     (2) For each nonlinear inequality constraint,  compute
C         the slack variable for which the merit function is
C         minimized.
C     (3) Compute the search direction for the slack variables
C         and multipliers.
C
C     Note that the array VIOLN is WRK3.
C     ============================================================
      CALL E04UCU(FEASQP,UNITQ,NQPERR,MAJITS,MNR,N,NCLIN,NCNLN,LDCJ,
     *            LDAQP,LDR,LINACT,NLNACT,NACTIV,NFREE,NZ,NUMINF,ISTATE,
     *            KACTIV,KX,DXNORM,GDX,QPCURV,AQP,W(LADX),W(LANORM),AX,
     *            BL,BU,C,CJAC,CLAMDA,W(LCMUL),W(LCS1),W(LDLAM),W(LDSLK)
     *            ,W(LDX),W(LBL),W(LBU),W(LQPTOL),R,W(LRHO),W(LSLK),
     *            W(LVIOLN),X,W(LWTINF),W)
C
      MINITS = MINITS + MNR
      MNRSUM = MNRSUM + MNR
C
      IF (FEASQP) THEN
         NQPINF = 0
      ELSE
         NQPINF = NQPINF + 1
         MJRMSG(2:2) = 'Infeasible subproblem'
      END IF
C
C     ============================================================
C     Compute quantities needed for the convergence test.
C     ============================================================
C     Compute the norms of the projected gradient and the
C     gradient with respect to the free variables.
C
      GZNORM = ZERO
      IF (NZ.GT.0) GZNORM = DNRM2(NZ,W(LGQ),1)
      GFNORM = GZNORM
      IF (NFREE.GT.0 .AND. NACTIV.GT.0) GFNORM = DNRM2(NFREE,W(LGQ),1)
C
C     If the forward-difference estimate of the transformed
C     gradient of the Lagrangian function is small,  switch to
C     central differences, recompute the derivatives and re-solve
C     the QP.
C
      GOODGQ = .TRUE.
      IF (NEEDFD .AND. .NOT. CENTRL) THEN
         GLNORM = DNRM2(N,W(LHPQ),1)
         IF (NCNLN.EQ.0) THEN
            CNORM = ZERO
         ELSE
            CNORM = DNRM2(NCNLN,C,1)
         END IF
C
         GLTEST = (ONE+ABS(OBJF)+ABS(CNORM))*EPSRF/FDNORM
         IF (GLNORM.LE.GLTEST) THEN
            GOODGQ = .FALSE.
            MJRMSG(3:3) = 'Central differences'
            LVLDIF = 2
            NEWGQ = .TRUE.
            IF (MSGNP.GE.5) THEN
               IF (MINITS.GT.0) THEN
                  WRITE (REC,FMT=99999) MINITS
                  CALL X04BAF(IPRINT,REC(1))
               END IF
            END IF
         END IF
C
      END IF
C
C     +       UNTIL     (GOODGQ)
      IF ( .NOT. GOODGQ) GO TO 40
C
C     ===============================================================
C     (1) Compute the number of constraints that are violated by more
C         than FEATOL.
C     (2) Compute the 2-norm of the residuals of the constraints in
C         the QP working set.
C     ===============================================================
      CALL E04UCW(N,NCLIN,NCNLN,ISTATE,BIGBND,CVNORM,ERRMAX,JMAX,NVIOL,
     *            AX,BL,BU,C,FEATOL,X,W(LWRK2))
C
C     Define small quantities that reflect the magnitude of OBJF and
C     the norm of GRAD(free).
C
      OBJSIZ = ONE + ABS(OBJF)
      XSIZE = ONE + XNORM
      GTEST = MAX(OBJSIZ,GFNORM)
      DINKY = RTFTOL*GTEST
C
      IF (NACTIV.EQ.0) THEN
         CONDT = ZERO
      ELSE IF (NACTIV.EQ.1) THEN
         CONDT = DTMIN
      ELSE
         CONDT = F06BLF(DTMAX,DTMIN,OVERFL)
      END IF
C
      CALL F06FLF(N,R,LDR+1,DRMAX,DRMIN)
C
      CONDH = F06BLF(DRMAX,DRMIN,OVERFL)
      IF (CONDH.LT.RTMAX) THEN
         CONDH = CONDH*CONDH
      ELSE
         CONDH = FLMAX
      END IF
C
      IF (NZ.EQ.0) THEN
         CONDHZ = ONE
      ELSE IF (NZ.EQ.N) THEN
         CONDHZ = CONDH
      ELSE
         CALL F06FLF(NZ,R,LDR+1,DRZMAX,DRZMIN)
         CONDHZ = F06BLF(DRZMAX,DRZMIN,OVERFL)
         IF (CONDHZ.LT.RTMAX) THEN
            CONDHZ = CONDHZ*CONDHZ
         ELSE
            CONDHZ = FLMAX
         END IF
      END IF
C
C     ---------------------------------------------------------------
C     Test for convergence.
C     The point test CONVPT checks for a K-T point at the initial
C     point or after a large change in X.
C     ---------------------------------------------------------------
      CONVPT = DXNORM .LE. EPSPT8*GTEST .AND. NVIOL .EQ. 0 .AND.
     *         NQPERR .LE. 1
C
      KTCOND(1) = GZNORM .LT. DINKY
      KTCOND(2) = NVIOL .EQ. 0
      OPTIML = KTCOND(1) .AND. KTCOND(2)
C
      CONVRG = MAJITS .GT. 0 .AND. ALFDX .LE. RTFTOL*XSIZE
C
      INFEAS = CONVRG .AND. .NOT. FEASQP .OR. NQPINF .GT. 7
C
      DONE = CONVPT .OR. (CONVRG .AND. OPTIML) .OR. INFEAS
C
      OBJALF = OBJF
      GRDALF = GDX
      GLF1 = GDX
      IF (NCNLN.GT.0) THEN
         GLF1 = GLF1 - DDOT(NCNLN,W(LCJDX),1,CLAMDA(NL),1)
C
C        Compute the value and directional derivative of the
C        augmented Lagrangian merit function.
C        The penalty parameters may be increased or decreased.
C
         CALL E04UCN(FEASQP,N,NCLIN,NCNLN,OBJALF,GRDALF,QPCURV,ISTATE,
     *               W(LCJDX),W(LCMUL),W(LCS1),W(LDLAM),W(LRHO),
     *               W(LVIOLN),W(LWRK1),W(LWRK2))
      END IF
C
C     ===============================================================
C     Print the details of this iteration.
C     ===============================================================
      CALL E04UCT(KTCOND,CONVRG,MJRMSG,MSGNP,MSGQP,LDR,LDT,N,NCLIN,
     *            NCNLN,NCTOTL,NACTIV,LINACT,NLNACT,NZ,NFREE,MAJIT0,
     *            MAJITS,MINITS,ISTATE,ALFA,NFUN,CONDHZ,CONDH,CONDT,
     *            OBJALF,OBJF,GZNORM,CVNORM,AX,C,R,W(LT),W(LVIOLN),X,
     *            W(LWRK1))
C
      ALFA = ZERO
      ERROR = MAJITS .GE. NMAJOR
C
      IF ( .NOT. (DONE .OR. ERROR)) THEN
         MAJITS = MAJITS + 1
C
C        Make copies of information needed for the BFGS update.
C
         CALL DCOPY(N,X,1,W(LX1),1)
         CALL DCOPY(N,W(LGQ),1,W(LGQ1),1)
C
         IF (NCNLN.GT.0) THEN
            CALL DCOPY(NCNLN,W(LCJDX),1,W(LCJDX1),1)
            CALL DCOPY(NCNLN,W(LCMUL),1,W(LC1MUL),1)
            CALL DCOPY(NCNLN,W(LSLK),1,W(LSLK1),1)
         END IF
C
C        ============================================================
C        Compute the parameters for the linesearch.
C        ============================================================
C        ALFMIN is the smallest allowable step predicted by the QP
C        subproblem.
C
         ALFMIN = ONE
         IF ( .NOT. FEASQP) ALFMIN = ZERO
C
C        ------------------------------------------------------------
C        ALFMAX is the largest feasible steplength subject to a user-
C        defined limit ALFLIM on the change in X.
C        ------------------------------------------------------------
         IF (NCNLN.GT.0 .AND. NEEDFD) THEN
            ALFMAX = ONE
         ELSE
            ALFMAX = F06BLF(BIGDX,DXNORM,OVERFL)
            CALL E04UDT(INFO,N,NCLIN,NCNLN,ALFA,ALFMIN,ALFMAX,BIGBND,
     *                  DXNORM,W(LANORM),W(LADX),AX,BL,BU,W(LDSLK),
     *                  W(LDX),W(LSLK),X)
            ALFMAX = ALFA
            IF (ALFMAX.LT.ONE+EPSPT3 .AND. FEASQP) ALFMAX = ONE
         END IF
C
C        ------------------------------------------------------------
C        ALFBND is a tentative upper bound on the steplength.  If the
C        merit function is decreasing at ALFBND and certain
C        conditions hold,  ALFBND will be increased in multiples of
C        two (subject to not being greater than ALFMAX).
C        ------------------------------------------------------------
         IF (NCNLN.EQ.0) THEN
            ALFBND = ALFMAX
         ELSE
            ALFBND = MIN(ONE,ALFMAX)
         END IF
C
C        ------------------------------------------------------------
C        ALFSML trips the computation of central differences.  If a
C        trial steplength falls below ALFSML, the linesearch is
C        terminated.
C        ------------------------------------------------------------
         ALFSML = ZERO
         IF (NEEDFD .AND. .NOT. CENTRL) THEN
            ALFSML = F06BLF(FDNORM,DXNORM,OVERFL)
            ALFSML = MIN(ALFSML,ALFMAX)
         END IF
C
C        ============================================================
C        Compute the steplength using safeguarded interpolation.
C        ============================================================
         ALFLIM = F06BLF((ONE+XNORM)*DXLIM,DXNORM,OVERFL)
         ALFA = MIN(ALFLIM,ONE)
C
         CALL E04UCR(NEEDFD,NLSERR,N,NCNLN,LDCJ,LDCJU,NFUN,NGRAD,
     *               IW(LNEEDC),CONFUN,OBJFUN,ALFA,ALFBND,ALFMAX,ALFSML,
     *               DXNORM,EPSRF,ETA,GDX,GRDALF,GLF1,GLF2,OBJF,OBJALF,
     *               QPCURV,XNORM,C,W(LWRK1),CJAC,CJACU,W(LCJDX),
     *               W(LWRK3),W(LC1MUL),W(LCMUL),W(LCS1),W(LCS2),W(LDX),
     *               W(LDLAM),W(LDSLK),GRAD,GRADU,CLAMDA(NL),W(LRHO),
     *               W(LSLK1),W(LSLK),W(LX1),X,W(LWRK2),IUSER,USER)
C
C           ------------------------------------------------------------
C           E04UCR  sets NLSERR to the following values...
C
C           < 0  if the user wants to stop.
C             1  if the search is successful and ALFA < ALFMAX.
C             2  if the search is successful and ALFA = ALFMAX.
C             3  if a better point was found but too many functions
C                were needed (not sufficient decrease).
C
C           Values of NLSERR occurring with a nonzero value of ALFA.
C             4  if ALFMAX < TOLABS (too small to do a search).
C             5  if ALFA  < ALFSML (E04UCJ only -- maybe want to switch
C                to central differences to get a better direction).
C             6  if the search found that there is no useful step.
C                The interval of uncertainty is less than 2*TOLABS.
C                The minimizer is very close to ALFA = zero
C                or the gradients are not sufficiently accurate.
C             7  if there were too many function calls.
C             8  if the input parameters were bad
C                (ALFMAX le TOLTNY  or  uphill).
C           ------------------------------------------------------------
         IF (NLSERR.LT.0) THEN
            INFORM = NLSERR
            GO TO 60
         END IF
C
         IF (ALFA.GT.ALFLIM) MJRMSG(4:4) = 'L'
C
         ERROR = NLSERR .GE. 4
         IF (ERROR) THEN
C           ---------------------------------------------------------
C           The linesearch failed to find a better point.
C           If exact gradients or central differences are being used,
C           or the KT conditions are satisfied, stop.  Otherwise,
C           switch to central differences and solve the QP again.
C           ---------------------------------------------------------
            IF (NEEDFD .AND. .NOT. CENTRL) THEN
               IF ( .NOT. OPTIML) THEN
                  ERROR = .FALSE.
                  MJRMSG(3:3) = 'Central differences'
                  LVLDIF = 2
                  NEWGQ = .TRUE.
                  IF (MSGNP.GE.5) THEN
                     IF (MINITS.GT.0) THEN
                        WRITE (REC,FMT=99999) MINITS
                        CALL X04BAF(IPRINT,REC(1))
                     END IF
                  END IF
               END IF
            END IF
         ELSE
            IF (NEEDFD) THEN
C              ======================================================
C              Compute the missing gradients.
C              ======================================================
               MODE = 1
               NGRAD = NGRAD + 1
C
               IF (NCNLN.GT.0) THEN
                  CALL F06DBF(NCNLN,(1),IW(LNEEDC),1)
C
                  CALL CONFUN(MODE,NCNLN,N,LDCJU,IW(LNEEDC),X,W(LWRK1),
     *                        CJACU,NSTATE,IUSER,USER)
                  INFORM = MODE
                  IF (MODE.LT.0) GO TO 60
C
                  CALL F06QFF('General',NCNLN,N,CJACU,LDCJU,CJAC,LDCJ)
               END IF
C
               CALL OBJFUN(MODE,N,X,OBJ,GRADU,NSTATE,IUSER,USER)
               INFORM = MODE
               IF (MODE.LT.0) GO TO 60
C
               CALL DCOPY(N,GRADU,1,GRAD,1)
C
               CALL E04UDS(CENTRL,MODE,LDCJ,LDCJU,N,NCNLN,BIGBND,CDINT,
     *                     FDINT,FDNORM,OBJF,CONFUN,OBJFUN,IW(LNEEDC),
     *                     BL,BU,C,W(LWRK2),W(LWRK3),CJAC,CJACU,GRAD,
     *                     GRADU,W(LHFRWD),W(LHCTRL),X,IUSER,USER)
C
               INFORM = MODE
               IF (MODE.LT.0) GO TO 60
C
               GDX = DDOT(N,GRAD,1,W(LDX),1)
               GLF2 = GDX
               IF (NCNLN.GT.0) THEN
                  CALL DGEMV('N',NCNLN,N,ONE,CJAC,LDCJ,W(LDX),1,ZERO,
     *                       W(LCJDX),1)
                  GLF2 = GLF2 - DDOT(NCNLN,W(LCJDX),1,CLAMDA(NL),1)
               END IF
            END IF
C
            CALL DCOPY(N,GRAD,1,W(LGQ),1)
            CALL E04NBW(6,N,NZ,NFREE,LDZY,UNITQ,KX,W(LGQ),W(LZY),
     *                  W(LWRK1))
C
            XNORM = DNRM2(N,X,1)
C
            IF (NCNLN.GT.0 .AND. ALFA.GE.ONE) CALL DCOPY(NCNLN,
     *          CLAMDA(NL),1,W(LCMUL),1)
C
            IF (NCLIN.GT.0) CALL DAXPY(NCLIN,ALFA,W(LADX),1,AX,1)
            ALFDX = ALFA*DXNORM
C
C           =========================================================
C           Update the factors of the approximate Hessian of the
C           Lagrangian function.
C           =========================================================
            CALL E04UCL(MJRMSG,UNITQ,N,NCNLN,NFREE,NZ,LDCJ1,LDCJ,LDZY,
     *                  LDR,KX,ALFA,GLF1,GLF2,QPCURV,W(LCJAC1),CJAC,
     *                  W(LCJDX1),W(LCJDX),W(LCS1),W(LCS2),W(LGQ1),
     *                  W(LGQ),W(LHPQ),W(LRPQ),CLAMDA(NL),R,W(LWRK3),
     *                  W(LZY),W(LWRK2),W(LWRK1))
C
            CALL F06FLF(N,R,LDR+1,DRMAX,DRMIN)
            COND = F06BLF(DRMAX,DRMIN,OVERFL)
C
            IF (COND.GT.RCNDBD .OR. RFROBN.GT.ROOTN*GROWTH*DRMAX) THEN
C              ------------------------------------------------------
C              Reset the condition estimator and range-space
C              partition of Q'HQ.
C              ------------------------------------------------------
               MJRMSG(5:5) = 'Refactorize Hessian'
C
               CALL E04UDR(UNITQ,N,NFREE,NZ,LDZY,LDR,IW(LIPERM),KX,
     *                     W(LGQ),R,W(LZY),W(LWRK1),W(LQRWRK))
            END IF
         END IF
      END IF
C
C     +    UNTIL     (DONE  .OR.  ERROR)
      IF ( .NOT. (DONE .OR. ERROR)) GO TO 20
C
C     ======================end of main loop============================
C
      IF (DONE) THEN
         IF (CONVRG .AND. OPTIML) THEN
            INFORM = 0
         ELSE IF (CONVPT) THEN
            INFORM = 1
         ELSE IF (INFEAS) THEN
            INFORM = 3
         END IF
      ELSE IF (ERROR) THEN
         IF (MAJITS.GE.NMAJOR) THEN
            INFORM = 4
         ELSE IF (OPTIML) THEN
            INFORM = 1
         ELSE
            INFORM = 6
         END IF
      END IF
C
C     ------------------------------------------------------------------
C     Set  CLAMDA.  Print the full solution.
C     ------------------------------------------------------------------
   60 IF (MSGNP.GT.0) THEN
         WRITE (REC,FMT=99998) MAJITS, MNRSUM
         CALL X04BAY(IPRINT,3,REC)
      END IF
C
      CALL E04UPP(NFREE,LDAQP,N,NCLIN,NCTOTL,NACTIV,ISTATE,KACTIV,KX,
     *            AQP,BL,BU,C,CLAMDA,FEATOL,W(LWRK1),W(LRLAM),X)
      CALL E04MFK(MSGNP,N,NCLIN,NCTOTL,BIGBND,NAMED,NAMES,ISTATE,BL,BU,
     *            CLAMDA,FEATOL,W(LWRK1))
C
      IF (NCNLN.GT.0) CALL DCOPY(NCNLN,W(LCMUL),1,CLAMDA(N+NCLIN+1),1)
C
      RETURN
C
C
C
C     End of  E04UCZ. (NPCORE)
C
99999 FORMAT (' Mnr itn ',I4,' -- Re-solve QP subproblem.')
99998 FORMAT (/' Exit from NP problem after ',I5,' major iterations,',
     *       /'                            ',I5,' minor iterations.')
      END
