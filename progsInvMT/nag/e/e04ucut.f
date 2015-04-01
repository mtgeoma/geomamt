      SUBROUTINE E04UCU(FEASQP,UNITQ,NQPERR,MAJITS,MINITS,N,NCLIN,NCNLN,
     *                  LDCJ,LDAQP,LDR,LINACT,NLNACT,NACTIV,NFREE,NZ,
     *                  NUMINF,ISTATE,KACTIV,KX,DXNORM,GDX,QPCURV,AQP,
     *                  ADX,ANORM,AX,BL,BU,C,CJAC,CLAMDA,CMUL,CS,DLAM,
     *                  DSLK,DX,QPBL,QPBU,QPTOL,R,RHO,SLK,VIOLN,X,WTINF,
     *                  W)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 14C REVISED. IER-881 (NOV 1990).
C     MARK 16 REVISED. IER-1089 (JUL 1993).
C     MARK 17 REVISED. IER-1609 (JUN 1995).
C
C     ******************************************************************
C     E04UCU   does the following:
C
C     (1)  Generate the upper and lower bounds for the QP  subproblem.
C
C     (2)  Compute the  TQ  factors of the rows of  AQP  specified by
C          the array  ISTATE.  The part of the factorization defined by
C          the first contiguous group of linear constraints does not
C          need to be recomputed.  The remaining rows (which could be
C          comprised of both linear and nonlinear constraints) are
C          included as new rows of the  TQ  factorization stored in
C          T and ZY.  Note that if there are no nonlinear constraints,
C          no factorization is required.
C
C     (3)  Solve the  QP  subproblem.
C                 minimize     1/2 (W p - d)'(Wp - d) + g'p
C
C                 subject to   qpbl .le. (  p ) .le. qpbu,
C                                        ( Ap )
C
C          where  W  is a matrix (not stored) such that  W'W = H  and
C          WQ = R,  d  is the zero vector,  and  g  is the gradient.
C          If the subproblem is infeasible, compute the point which
C          minimizes the sum of infeasibilities.
C
C     (4)   Find the value of each slack variable for which the merit
C          function is minimized.
C
C     (5)   Compute  DSLK,  DLAM  and  DX,  the search directions for
C          the slack variables, the multipliers and the variables.
C
C     Systems Optimization Laboratory, Stanford University.
C     Fortran 66 version written 10-January-1983.
C     Level-2 matrix routines added 18-May-1988.
C     This version of E04UCU dated 23-Dec-92.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LENLS
      PARAMETER         (LENLS=20)
      INTEGER           MXPARM
      PARAMETER         (MXPARM=30)
      LOGICAL           QPNAMD, VERTEX
      PARAMETER         (QPNAMD=.FALSE.,VERTEX=.FALSE.)
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0,TWO=2.0D+0)
      DOUBLE PRECISION  HUNDRD
      PARAMETER         (HUNDRD=1.0D+2)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DXNORM, GDX, QPCURV
      INTEGER           LDAQP, LDCJ, LDR, LINACT, MAJITS, MINITS, N,
     *                  NACTIV, NCLIN, NCNLN, NFREE, NLNACT, NQPERR,
     *                  NUMINF, NZ
      LOGICAL           FEASQP, UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  ADX(*), ANORM(*), AQP(LDAQP,*), AX(*), BL(*),
     *                  BU(*), C(*), CJAC(LDCJ,*), CLAMDA(*), CMUL(*),
     *                  CS(*), DLAM(*), DSLK(*), DX(N), QPBL(*),
     *                  QPBU(*), QPTOL(*), R(LDR,*), RHO(*), SLK(*),
     *                  VIOLN(*), W(*), WTINF(*), X(N)
      INTEGER           ISTATE(*), KACTIV(N), KX(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  ASIZE, BIGBND, BIGDX, BNDLOW, BNDUPP, CDINT,
     *                  CTOL, DTMAX, DTMIN, DXLIM, EPSPT3, EPSPT5,
     *                  EPSPT8, EPSPT9, EPSRF, ETA, FDINT, FTOL, HCNDBD,
     *                  RHODMP, RHOMAX, RHONRM, SCALE, TOLACT, TOLFEA,
     *                  TOLRNK
      INTEGER           IPRINT, IPRNT, ISUMM, ISUMRY, ITMAX1, ITMAX2,
     *                  ITMXNP, JVRFY1, JVRFY2, JVRFY3, JVRFY4, KSAVE,
     *                  LCRASH, LDT, LDZY, LENNAM, LFORMH, LINES1,
     *                  LINES2, LPROB, LVERFY, LVLDER, MSGLS, MSGNP,
     *                  NCOLT, NLNF, NLNJ, NLNX, NLOAD, NN, NNCLIN,
     *                  NNCNLN, NOUT, NPROB, NSAVE
      LOGICAL           INCRUN
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADLS(23), RPADNP(22), RPSVLS(MXPARM),
     *                  RPSVNP(MXPARM)
      INTEGER           IPADLS(19), IPADNP(15), IPSVLS(MXPARM),
     *                  IPSVNP(MXPARM), LOCLS(LENLS)
C     .. Local Scalars ..
      DOUBLE PRECISION  AMIN, BIGLOW, BIGUPP, BLJ, BUJ, CON, CONDMX,
     *                  QUOTNT, SSQ, SSQ1, SUMINF, VIOL, WEIGHT, WSCALE,
     *                  WTMAX, WTMIN
      INTEGER           I, INFORM, ISWAP, J, JINF, K, K1, K2, KVIOL, L,
     *                  LGQ, LHPQ, LRLAM, LRPQ, LRPQ0, LT, LWRK1, LZY,
     *                  MSGQP, NARTIF, NCQP, NCTOTL, NGQ, NMAJOR,
     *                  NMINOR, NPLIN, NRANK, NREJTD, NRPQ, NTRY, NVIOL,
     *                  NZ1
      LOGICAL           LINOBJ, OVERFL
C     .. Local Arrays ..
      DOUBLE PRECISION  RPRMLS(MXPARM), RPRMNP(MXPARM)
      INTEGER           IPRMLS(MXPARM), IPRMNP(MXPARM)
      CHARACTER*8       NAMES(1)
      CHARACTER*80      REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2, F06BLF
      EXTERNAL          DDOT, DNRM2, F06BLF
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DSCAL, DTRMV, DTRSV,
     *                  E04NBW, E04NCY, E04NCZ, E04UCM, F06DBF, F06FBF,
     *                  F06FLF, F06QFF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AE04NC/LOCLS
      COMMON            /BE04NB/LENNAM, LDT, NCOLT, LDZY
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /DE04NB/ASIZE, DTMAX, DTMIN
      COMMON            /DE04NC/IPSVLS, IPRNT, ISUMRY, ITMAX1, ITMAX2,
     *                  LCRASH, LFORMH, LPROB, MSGLS, NN, NNCLIN, NPROB,
     *                  IPADLS
      COMMON            /DE04UC/RHOMAX, RHONRM, RHODMP, SCALE, INCRUN
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
      SAVE              /DE04NC/, /EE04NC/, /GE04UC/, /HE04UC/
C     .. Executable Statements ..
C
      LRPQ = LOCLS(5)
      LRPQ0 = LOCLS(6)
      LHPQ = LOCLS(8)
      LGQ = LOCLS(9)
      LRLAM = LOCLS(10)
      LT = LOCLS(11)
      LZY = LOCLS(12)
      LWRK1 = LOCLS(14)
C
      NRPQ = 0
      NGQ = 1
C
      FEASQP = .TRUE.
      LINOBJ = .TRUE.
C
      BIGLOW = -BIGBND
      BIGUPP = BIGBND
      SSQ1 = ZERO
C
      NPLIN = N + NCLIN
      NCTOTL = NPLIN + NCNLN
      NCQP = NCLIN + NCNLN
      NRANK = N
      NREJTD = 0
C
      IF (MSGQP.GT.0) THEN
         WRITE (REC,FMT=99999) MAJITS
         CALL X04BAY(IPRINT,3,REC)
      END IF
C
C     ==================================================================
C     Generate the upper and lower bounds upon the search direction, the
C     weights on the sum of infeasibilities and the nonlinear constraint
C     violations.
C     ==================================================================
      WSCALE = -ONE
      DO 40 J = 1, NCTOTL
C
         IF (J.LE.N) THEN
            CON = X(J)
         ELSE IF (J.LE.NPLIN) THEN
            CON = AX(J-N)
         ELSE
            CON = C(J-NPLIN)
         END IF
C
         BLJ = BL(J)
         BUJ = BU(J)
         IF (BLJ.GT.BIGLOW) BLJ = BLJ - CON
         IF (BUJ.LT.BIGUPP) BUJ = BUJ - CON
C
         WEIGHT = ONE
         IF (J.LE.NPLIN) THEN
            IF (ABS(BLJ).LE.QPTOL(J)) BLJ = ZERO
            IF (ABS(BUJ).LE.QPTOL(J)) BUJ = ZERO
         ELSE
            I = J - NPLIN
            VIOL = ZERO
            IF (BL(J).GT.BIGLOW) THEN
               IF (BLJ.GT.ZERO) THEN
                  VIOL = BLJ
                  IF (RHO(I).GT.ZERO) THEN
                     WEIGHT = VIOL*RHO(I)
                  ELSE
                     WEIGHT = VIOL
                  END IF
                  WSCALE = MAX(WSCALE,WEIGHT)
                  GO TO 20
               END IF
            END IF
C
            IF (BU(J).LT.BIGUPP) THEN
               IF (BUJ.LT.ZERO) THEN
                  VIOL = BUJ
                  IF (RHO(I).GT.ZERO) THEN
                     WEIGHT = -VIOL*RHO(I)
                  ELSE
                     WEIGHT = -VIOL
                  END IF
                  WSCALE = MAX(WSCALE,WEIGHT)
               END IF
            END IF
C
C           Set the vector of nonlinear constraint violations.
C
   20       VIOLN(I) = VIOL
         END IF
C
         WTINF(J) = WEIGHT
         QPBL(J) = BLJ
         QPBU(J) = BUJ
C
   40 CONTINUE
C
      IF (WSCALE.GT.ZERO) THEN
         WSCALE = ONE/WSCALE
         CALL DSCAL(NCTOTL,(WSCALE),WTINF,1)
      END IF
C
      CALL F06FLF(NCTOTL,WTINF,1,WTMAX,WTMIN)
      WTMIN = EPSPT9*WTMAX
      DO 60 J = 1, NCTOTL
         WTINF(J) = MAX(WTINF(J),WTMIN)
   60 CONTINUE
C
C     Set the maximum allowable condition estimator of the constraints
C     in the working set.  Note that a relatively well-conditioned
C     working set is used to start the QP iterations.
C
      CONDMX = MAX(ONE/EPSPT3,HUNDRD)
C
      IF (NCNLN.GT.0) THEN
C        ===============================================================
C        Refactorize part of the  QP  constraint matrix.
C        ===============================================================
C        Load the new Jacobian into the  QP  matrix  A.  Compute the
C        2-norms of the rows of the Jacobian.
C
         CALL F06QFF('General',NCNLN,N,CJAC,LDCJ,AQP(NCLIN+1,1),LDAQP)
C
         DO 80 J = NCLIN + 1, NCQP
            ANORM(J) = DNRM2(N,AQP(J,1),LDAQP)
   80    CONTINUE
C
C        Count the number of linear constraints in the working set and
C        move them to the front of KACTIV.  Compute the norm of the
C        matrix of constraints in the working set.
C        Let K1  point to the first nonlinear constraint.  Constraints
C        with indices KACTIV(K1),..., KACTIV(NACTIV)  must be
C        refactorized.
C
         ASIZE = ZERO
         LINACT = 0
         K1 = NACTIV + 1
         DO 100 K = 1, NACTIV
            I = KACTIV(K)
            ASIZE = MAX(ASIZE,ANORM(I))
C
            IF (I.LE.NCLIN) THEN
               LINACT = LINACT + 1
               IF (LINACT.NE.K) THEN
                  ISWAP = KACTIV(LINACT)
                  KACTIV(LINACT) = I
                  KACTIV(K) = ISWAP
               END IF
            ELSE
C
C              Record the old position of the 1st. nonlinear constraint.
C
               IF (K1.GT.NACTIV) K1 = K
            END IF
  100    CONTINUE
C
         IF (NACTIV.LE.1) CALL F06FLF(NCQP,ANORM,1,ASIZE,AMIN)
C
C        Compute the absolute values of the nonlinear constraints in
C        the working set.  Use DX as workspace.
C
         DO 120 K = LINACT + 1, NACTIV
            J = N + KACTIV(K)
            IF (ISTATE(J).EQ.1) DX(K) = ABS(QPBL(J))
            IF (ISTATE(J).GE.2) DX(K) = ABS(QPBU(J))
  120    CONTINUE
C
C        Sort the elements of KACTIV corresponding to nonlinear
C        constraints in descending order of violation (i.e.,
C        the first element of KACTIV for a nonlinear constraint
C        is associated with the most violated constraint.)
C        In this way, the rows of the Jacobian corresponding
C        to the more violated constraints tend to be included
C        in the  TQ  factorization.
C
C        The sorting procedure is taken from the simple insertion
C        sort in D. Knuth, ACP Volume 3, Sorting and Searching,
C        Page 81.  It should be replaced by a faster sort if the
C        number of active nonlinear constraints becomes large.
C
         DO 160 K = LINACT + 2, NACTIV
            L = K
            VIOL = DX(L)
            KVIOL = KACTIV(L)
C           WHILE (L .GT. LINACT+1  .AND.  DX(L-1) .LT. VIOL) DO
  140       IF (L.GT.LINACT+1) THEN
               IF (DX(L-1).LT.VIOL) THEN
                  DX(L) = DX(L-1)
                  KACTIV(L) = KACTIV(L-1)
                  L = L - 1
                  GO TO 140
               END IF
C              END WHILE
            END IF
            DX(L) = VIOL
            KACTIV(L) = KVIOL
  160    CONTINUE
C
         K2 = NACTIV
         NACTIV = K1 - 1
         NZ = NFREE - NACTIV
C
C        Update the factors  R,  T  and  Q  to include constraints
C        K1  through  K2.
C
         IF (K1.LE.K2) CALL E04NCY(UNITQ,VERTEX,INFORM,K1,K2,NACTIV,
     *                             NARTIF,NZ,NFREE,NRANK,NREJTD,NRPQ,
     *                             NGQ,N,LDZY,LDAQP,LDR,LDT,ISTATE,
     *                             KACTIV,KX,CONDMX,AQP,R,W(LT),W(LRPQ),
     *                             W(LGQ),W(LZY),W(LWRK1),DX,W(LRLAM),
     *                             MSGQP)
      END IF
C
C     ==================================================================
C     Solve for DX, the vector of minimum two-norm that satisfies the
C     constraints in the working set.
C     ==================================================================
      CALL E04UCM(UNITQ,NCQP,NACTIV,NFREE,NZ,N,NLNX,NCTOTL,LDZY,LDAQP,
     *            LDR,LDT,ISTATE,KACTIV,KX,DXNORM,GDX,AQP,ADX,QPBL,QPBU,
     *            W(LRPQ),W(LRPQ0),DX,W(LGQ),R,W(LT),W(LZY),W(LWRK1))
C
C     ==================================================================
C     Solve a quadratic program for the search direction  DX  and
C     multiplier estimates  CLAMDA.
C     ==================================================================
C     If there is no feasible point for the subproblem,  the sum of
C     infeasibilities is minimized subject to the linear constraints
C     (1  thru  JINF)  being satisfied.
C
      JINF = N + NCLIN
C
      NTRY = 1
C     +    REPEAT
  180 CALL E04NCZ('QP subproblem',QPNAMD,NAMES,LINOBJ,UNITQ,NQPERR,
     *            MINITS,JINF,NCQP,NCTOTL,NACTIV,NFREE,NRANK,NZ,NZ1,N,
     *            LDAQP,LDR,ISTATE,KACTIV,KX,GDX,SSQ,SSQ1,SUMINF,NUMINF,
     *            DXNORM,QPBL,QPBU,AQP,CLAMDA,ADX,QPTOL,R,DX,W)
C
      NVIOL = 0
      IF (NUMINF.GT.0) THEN
C
C           Count the violated linear constraints.
C
         DO 200 J = 1, NPLIN
            IF (ISTATE(J).LT.0) NVIOL = NVIOL + 1
  200    CONTINUE
C
         IF (NVIOL.GT.0) THEN
            NTRY = NTRY + 1
            UNITQ = .TRUE.
            NACTIV = 0
            NFREE = N
            NZ = N
            CALL F06DBF(NCTOTL,(0),ISTATE,1)
C
            CALL E04UCM(UNITQ,NCQP,NACTIV,NFREE,NZ,N,NLNX,NCTOTL,LDZY,
     *                  LDAQP,LDR,LDT,ISTATE,KACTIV,KX,DXNORM,GDX,AQP,
     *                  ADX,QPBL,QPBU,W(LRPQ),W(LRPQ0),DX,W(LGQ),R,W(LT)
     *                  ,W(LZY),W(LWRK1))
         END IF
      END IF
      IF ( .NOT. (NVIOL.EQ.0 .OR. NTRY.GT.2)) GO TO 180
C     +    UNTIL (    NVIOL .EQ. 0  .OR.  NTRY .GT. 2)
C
C     ==================================================================
C     Count the number of nonlinear constraint gradients in the  QP
C     working set.  Make sure that all small  QP  multipliers associated
C     with nonlinear inequality constraints have the correct sign.
C     ==================================================================
      NLNACT = 0
      IF (NACTIV.GT.0 .AND. NCNLN.GT.0) THEN
         DO 220 K = 1, NACTIV
            L = KACTIV(K)
            IF (L.GT.NCLIN) THEN
               NLNACT = NLNACT + 1
               J = N + L
               IF (ISTATE(J).EQ.1) CLAMDA(J) = MAX(ZERO,CLAMDA(J))
               IF (ISTATE(J).EQ.2) CLAMDA(J) = MIN(ZERO,CLAMDA(J))
            END IF
  220    CONTINUE
      END IF
C
      LINACT = NACTIV - NLNACT
C
C     ------------------------------------------------------------------
C     Extract various useful quantities from the QP solution.
C     ------------------------------------------------------------------
C     Compute  HPQ = R'R(pq)  from the transformed gradient of the QP
C     objective function and  R(pq)  from the transformed residual.
C
      CALL DSCAL(N,(-ONE),W(LRPQ),1)
      CALL DAXPY(N,(-ONE),W(LGQ),1,W(LHPQ),1)
      QPCURV = TWO*SSQ
C
      IF (NCNLN.GT.0) THEN
         IF (NUMINF.GT.0) THEN
            FEASQP = .FALSE.
            CALL F06FBF(NCTOTL,(ZERO),CLAMDA,1)
C
            IF (NZ.GT.0) THEN
C              ---------------------------------------------------------
C              Compute a null space element for the search direction
C              as the solution of  Z'HZ(pz) = -Z'g - Z'HY(py).
C              ---------------------------------------------------------
C              Overwrite DX with the transformed search direction
C              Q'(dx).  The first NZ elements of DX are zero.
C
               CALL E04NBW(6,N,NZ,NFREE,LDZY,UNITQ,KX,DX,W(LZY),W(LWRK1)
     *                     )
C
C              Overwrite the first NZ elements of DX with the solution
C              of  (Rz)u = -(v + w),  where  (Rz)'w = Z'g  and  v  is
C              vector of first NZ elements of  R(pq).
C
               CALL DCOPY(NZ,W(LGQ),1,DX,1)
               CALL DTRSV('U','T','N',NZ,R,LDR,DX,1)
C
               CALL DAXPY(NZ,(ONE),W(LRPQ),1,DX,1)
C
               CALL DTRSV('U','N','N',NZ,R,LDR,DX,1)
               CALL DSCAL(NZ,(-ONE),DX,1)
C
C              Recompute RPQ, HPQ, GDX and QPCURV.
C
               CALL DCOPY(NLNX,DX,1,W(LRPQ),1)
               CALL DTRMV('U','N','N',NLNX,R,LDR,W(LRPQ),1)
               IF (NLNX.LT.N) CALL DGEMV('N',NLNX,N-NLNX,ONE,R(1,NLNX+1)
     *                                   ,LDR,DX(NLNX+1),1,ONE,W(LRPQ),
     *                                   1)
C
               GDX = DDOT(N,W(LGQ),1,DX,1)
               QPCURV = DDOT(N,W(LRPQ),1,W(LRPQ),1)
C
               CALL E04NBW(3,N,NZ,NFREE,LDZY,UNITQ,KX,DX,W(LZY),W(LWRK1)
     *                     )
C
C              ---------------------------------------------------------
C              Recompute ADX and the 2-norm of DX.
C              ---------------------------------------------------------
               DXNORM = DNRM2(N,DX,1)
               IF (NCQP.GT.0) CALL DGEMV('N',NCQP,N,ONE,AQP,LDAQP,DX,1,
     *                                   ZERO,ADX,1)
C
            END IF
C
            CALL DCOPY(NLNX,W(LRPQ),1,W(LHPQ),1)
            CALL DTRMV('U','T','N',NLNX,R,LDR,W(LHPQ),1)
            IF (NLNX.LT.N) CALL DGEMV('T',NLNX,N-NLNX,ONE,R(1,NLNX+1),
     *                                LDR,W(LRPQ),1,ZERO,W(LHPQ+NLNX),1)
         END IF
C
C        ===============================================================
C        For given values of the objective function and constraints,
C        attempt to minimize the merit function with respect to each
C        slack variable.
C        ===============================================================
         DO 260 I = 1, NCNLN
            J = NPLIN + I
            CON = C(I)
C
            IF ( .NOT. FEASQP .AND. VIOLN(I).NE.ZERO .AND. RHO(I)
     *          .LE.ZERO) RHO(I) = ONE
C
            QUOTNT = F06BLF(CMUL(I),SCALE*RHO(I),OVERFL)
C
C           Define the slack variable to be  CON - MULT / RHO.
C           Force each slack to lie within its upper and lower bounds.
C
            IF (BL(J).GT.BIGLOW) THEN
               IF (QPBL(J).GE.-QUOTNT) THEN
                  SLK(I) = BL(J)
                  GO TO 240
               END IF
            END IF
C
            IF (BU(J).LT.BIGUPP) THEN
               IF (QPBU(J).LE.-QUOTNT) THEN
                  SLK(I) = BU(J)
                  GO TO 240
               END IF
            END IF
C
            SLK(I) = CON - QUOTNT
C
C           The slack has been set within its bounds.
C
  240       CS(I) = CON - SLK(I)
C
C           ------------------------------------------------------------
C           Compute the search direction for the slacks and multipliers.
C           ------------------------------------------------------------
            DSLK(I) = ADX(NCLIN+I) + CS(I)
C
            IF (FEASQP) THEN
C
C              If any constraint is such that  (DLAM)*(C - S)  is
C              positive,  the merit function may be reduced immediately
C              by substituting the QP multiplier.
C
               DLAM(I) = CLAMDA(J) - CMUL(I)
               IF (DLAM(I)*CS(I).GE.ZERO) THEN
                  CMUL(I) = CLAMDA(J)
                  DLAM(I) = ZERO
               END IF
            ELSE
C
C              The  QP  subproblem was infeasible.
C
               DLAM(I) = ZERO
C
               IF (ISTATE(J).LT.0 .OR. VIOLN(I).NE.ZERO) DSLK(I) = ZERO
C
            END IF
  260    CONTINUE
C
         IF ( .NOT. FEASQP) RHONRM = DNRM2(NCNLN,RHO,1)
C
      END IF
C
      RETURN
C
C
C     End of  E04UCU. (NPIQP)
C
99999 FORMAT (/1X,79('-'),/' Start of major itn',I6)
      END
