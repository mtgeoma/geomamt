      SUBROUTINE E04UCR(NEEDFD,INFORM,N,NCNLN,LDCJ,LDCJU,NFUN,NGRAD,
     *                  NEEDC,CONFUN,OBJFUN,ALFA,ALFBND,ALFMAX,ALFSML,
     *                  DXNORM,EPSRF,ETA,GDX,GRDALF,GLF1,GLF,OBJF,
     *                  OBJALF,QPCURV,XNORM,C,C2,CJAC,CJACU,CJDX,CJDX2,
     *                  CMUL1,CMUL,CS1,CS,DX,DLAM,DSLK,GRAD,GRADU,QPMUL,
     *                  RHO,SLK1,SLK,X1,X,WORK,IUSER,USER)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1086 (JUL 1993).
C     MARK 17 REVISED. IER-1606 (JUN 1995).
C
C     ==================================================================
C     E04UCR finds the steplength ALFA that gives sufficient decrease in
C     the augmented Lagrangian merit function.
C
C     On exit, if INFORM = 1, 2 or 3,  ALFA will be a nonzero steplength
C     with an associated merit function value  OBJALF  which is lower
C     than that at the base point. If  INFORM = 4, 5, 6, 7 or 8,  ALFA
C     is zero and  OBJALF will be the merit value at the base point.
C
C     Original version written  27-May-1985.
C     Level 2 BLAS added 12-June-1986.
C     This version of E04UCR dated  14-Sep-92.
C     ==================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF, ONE, TWO
      PARAMETER         (ZERO=0.0D+0,HALF=0.5D+0,ONE=1.0D+0,TWO=2.0D+0)
      DOUBLE PRECISION  TOLG
      PARAMETER         (TOLG=1.0D-1)
      DOUBLE PRECISION  RMU
      PARAMETER         (RMU=1.0D-4)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, ALFBND, ALFMAX, ALFSML, DXNORM, EPSRF,
     *                  ETA, GDX, GLF, GLF1, GRDALF, OBJALF, OBJF,
     *                  QPCURV, XNORM
      INTEGER           INFORM, LDCJ, LDCJU, N, NCNLN, NFUN, NGRAD
      LOGICAL           NEEDFD
C     .. Array Arguments ..
      DOUBLE PRECISION  C(*), C2(*), CJAC(LDCJ,*), CJACU(LDCJU,*),
     *                  CJDX(*), CJDX2(*), CMUL(*), CMUL1(*), CS(*),
     *                  CS1(*), DLAM(*), DSLK(*), DX(N), GRAD(N),
     *                  GRADU(N), QPMUL(*), RHO(*), SLK(*), SLK1(*),
     *                  USER(*), WORK(*), X(N), X1(N)
      INTEGER           IUSER(*), NEEDC(*)
C     .. Subroutine Arguments ..
      EXTERNAL          CONFUN, OBJFUN
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPSPT3, EPSPT5, EPSPT8, EPSPT9, RHODMP, RHOMAX,
     *                  RHONRM, SCALE
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
      LOGICAL           INCRUN
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALFBST, CS1JDX, CSJDX, CURVC, CURVLF, EPSAF,
     *                  EPSMCH, FBEST, FTERM, FTRY, G0, GBEST, GTRY,
     *                  OLDF, OLDG, Q, RHOBFS, S, T, TARGTG, TGDX, TGLF,
     *                  TOBJ, TOBJM, TOLABS, TOLAX, TOLREL, TOLRX,
     *                  TOLTNY
      INTEGER           J, MAXF, MODE, NSTATE, NUMF
      LOGICAL           DEBUG, DONE, FIRST, IMPRVD
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, E04UCJ, E04UCK, F06DBF,
     *                  F06FCF, F06QFF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /DE04UC/RHOMAX, RHONRM, RHODMP, SCALE, INCRUN
C     .. Save statement ..
      SAVE              /AX02ZA/
C     .. Executable Statements ..
C
      EPSMCH = WMACH(3)
C
      IF ( .NOT. NEEDFD .AND. NCNLN.GT.0) CS1JDX = DDOT(NCNLN,CS1,1,
     *    CJDX,1)
C
C     ------------------------------------------------------------------
C     Set the input parameters and tolerances for E04UCK and E04UCJ.
C
C     TOLRX   is the tolerance on relative changes in DX resulting from
C             changes in ALFA.
C
C     TOLAX   is the tolerance on absolute changes in DX resulting from
C             changes in ALFA.
C
C     TOLABS  is the tolerance on absolute changes in ALFA.
C
C     TOLREL  is the tolerance on relative changes in ALFA.
C
C     TOLTNY  is the magnitude of the smallest allowable value of ALFA.
C             if  M(TOLABS) - M(0) .gt. EPSAF,  the linesearch tries
C             steps in the range  TOLTNY .le. ALFA .le. TOLABS.
C     ------------------------------------------------------------------
      NSTATE = 0
      DEBUG = .FALSE.
C
      IF (NEEDFD) THEN
         MAXF = 15
      ELSE
         MAXF = 10
      END IF
C
      EPSAF = EPSRF*(ONE+ABS(OBJALF))
      TOLAX = EPSPT8
      TOLRX = EPSPT8
C
      IF (TOLRX*XNORM+TOLAX.LT.DXNORM*ALFMAX) THEN
         TOLABS = (TOLRX*XNORM+TOLAX)/DXNORM
      ELSE
         TOLABS = ALFMAX
      END IF
      TOLREL = MAX(TOLRX,EPSMCH)
C
      T = ZERO
      DO 20 J = 1, N
         S = ABS(DX(J))
         Q = ABS(X(J))*TOLRX + TOLAX
         IF (S.GT.T*Q) T = S/Q
   20 CONTINUE
C
      IF (T*TOLABS.GT.ONE) THEN
         TOLTNY = ONE/T
      ELSE
         TOLTNY = TOLABS
      END IF
C
      OLDF = OBJALF
      OLDG = GRDALF
      ALFBST = ZERO
      FBEST = ZERO
      GBEST = (ONE-RMU)*OLDG
      TARGTG = (RMU-ETA)*OLDG
      G0 = GBEST
C
      IF (NCNLN.GT.0) CALL F06DBF(NCNLN,(1),NEEDC,1)
C
      IF (NEEDFD) THEN
         MODE = 0
      ELSE
         MODE = 2
      END IF
C
      FIRST = .TRUE.
C
C     ------------------------------------------------------------------
C     Commence main loop, entering E04UCK or E04UCJ two or more times.
C     FIRST = .TRUE. for the first entry, .FALSE. for subsequent entries
C     DONE  = .TRUE. indicates termination, in which case the value of
C     INFORM gives the result of the search.
C     INFORM = 1 if the search is successful and ALFA < ALFMAX.
C            = 2 if the search is successful and ALFA = ALFMAX.
C            = 3 if a better point was found but too many functions
C                were needed (not sufficient decrease).
C            = 4 if ALFMAX < TOLABS (too small to do a search).
C            = 5 if ALFA < ALFSML (E04UCJ only -- maybe want to switch
C                to central differences to get a better direction).
C            = 6 if the search found that there is no useful step.
C                The interval of uncertainty is less than 2*TOLABS.
C                The minimizer is very close to ALFA = zero
C                or the gradients are not sufficiently accurate.
C            = 7 if there were too many function calls.
C            = 8 if the input parameters were bad
C                (ALFMAX le TOLTNY  or  OLDG ge 0).
C     ------------------------------------------------------------------
C     +    repeat
   40 IF (NEEDFD) THEN
         CALL E04UCJ(FIRST,DEBUG,DONE,IMPRVD,INFORM,MAXF,NUMF,IPRINT,
     *               ALFMAX,ALFSML,EPSAF,G0,TARGTG,FTRY,TOLABS,TOLREL,
     *               TOLTNY,ALFA,ALFBST,FBEST)
      ELSE
         CALL E04UCK(FIRST,DEBUG,DONE,IMPRVD,INFORM,MAXF,NUMF,IPRINT,
     *               ALFMAX,EPSAF,G0,TARGTG,FTRY,GTRY,TOLABS,TOLREL,
     *               TOLTNY,ALFA,ALFBST,FBEST,GBEST)
      END IF
C
      IF (IMPRVD) THEN
         OBJF = TOBJ
         OBJALF = TOBJM
C
         IF (NCNLN.GT.0) CALL DCOPY(NCNLN,C2,1,C,1)
C
         IF ( .NOT. NEEDFD) THEN
            CALL DCOPY(N,GRADU,1,GRAD,1)
            GDX = TGDX
            GLF = TGLF
C
            IF (NCNLN.GT.0) THEN
               CALL DCOPY(NCNLN,CJDX2,1,CJDX,1)
               CALL F06QFF('General',NCNLN,N,CJACU,LDCJU,CJAC,LDCJ)
            END IF
         END IF
      END IF
C
C     ---------------------------------------------------------------
C     If DONE = .FALSE.,  the problem functions must be computed for
C     the next entry to E04UCK or E04UCJ.
C     If DONE = .TRUE.,   this is the last time through.
C     ---------------------------------------------------------------
      IF ( .NOT. DONE) THEN
C
         CALL DCOPY(N,X1,1,X,1)
         CALL DAXPY(N,ALFA,DX,1,X,1)
C
         IF (NCNLN.GT.0) THEN
C
C           Compute the new estimates of the multipliers and slacks.
C           If the step length is greater than one,  the multipliers
C           are fixed as the QP-multipliers.
C
            IF (ALFA.LE.ONE) THEN
               CALL DCOPY(NCNLN,CMUL1,1,CMUL,1)
               CALL DAXPY(NCNLN,ALFA,DLAM,1,CMUL,1)
            END IF
            CALL DCOPY(NCNLN,SLK1,1,SLK,1)
            CALL DAXPY(NCNLN,ALFA,DSLK,1,SLK,1)
C
C           ---------------------------------------------------------
C           Compute the new constraint vector and Jacobian.
C           ---------------------------------------------------------
            CALL CONFUN(MODE,NCNLN,N,LDCJU,NEEDC,X,C2,CJACU,NSTATE,
     *                  IUSER,USER)
            IF (MODE.LT.0) GO TO 60
C
            CALL DCOPY(NCNLN,C2,1,CS,1)
            CALL DAXPY(NCNLN,(-ONE),SLK,1,CS,1)
C
            CALL DCOPY(NCNLN,CS,1,WORK,1)
            CALL F06FCF(NCNLN,RHO,1,WORK,1)
C
            FTERM = DDOT(NCNLN,CMUL,1,CS,1) - HALF*SCALE*DDOT(NCNLN,
     *              WORK,1,CS,1)
         END IF
C
C        ------------------------------------------------------------
C        Compute the value and gradient of the objective function.
C        ------------------------------------------------------------
         CALL OBJFUN(MODE,N,X,TOBJ,GRADU,NSTATE,IUSER,USER)
         IF (MODE.LT.0) GO TO 60
C
         IF (NCNLN.GT.0) THEN
            TOBJM = TOBJ - FTERM
         ELSE
            TOBJM = TOBJ
         END IF
C
         FTRY = TOBJM - OLDF - RMU*OLDG*ALFA
C
C        ------------------------------------------------------------
C        Compute auxiliary gradient information.
C        ------------------------------------------------------------
         IF ( .NOT. NEEDFD) THEN
            GTRY = DDOT(N,GRADU,1,DX,1)
            TGDX = GTRY
            TGLF = GTRY
            IF (NCNLN.GT.0) THEN
C
C              Compute the Jacobian times the search direction.
C
               CALL DGEMV('N',NCNLN,N,ONE,CJACU,LDCJU,DX,1,ZERO,CJDX2,1)
C
               CALL DCOPY(NCNLN,CJDX2,1,WORK,1)
               CALL DAXPY(NCNLN,(-ONE),DSLK,1,WORK,1)
C
               GTRY = GTRY - DDOT(NCNLN,CMUL,1,WORK,1)
               IF (ALFA.LE.ONE) GTRY = GTRY - DDOT(NCNLN,DLAM,1,CS,1)
C
               CALL F06FCF(NCNLN,RHO,1,WORK,1)
               GTRY = GTRY + SCALE*DDOT(NCNLN,WORK,1,CS,1)
C
               TGLF = TGDX - DDOT(NCNLN,CJDX2,1,QPMUL,1)
C
C              ------------------------------------------------------
C              If ALFBND .le. ALFA .lt. ALFMAX and the norm of the
C              quasi-Newton update is bounded, set ALFMAX to be ALFA.
C              This will cause the linesearch to stop if the merit
C              function is decreasing at the boundary.
C              ------------------------------------------------------
               IF (ALFBND.LE.ALFA .AND. ALFA.LT.ALFMAX) THEN
C
                  CSJDX = DDOT(NCNLN,CS,1,CJDX2,1)
C
                  CURVLF = TGLF - GLF1
                  CURVC = ABS(CSJDX-CS1JDX)
                  RHOBFS = MAX(QPCURV*TOLG-CURVLF,ZERO)
                  IF (RHOBFS.LE.CURVC*RHOMAX) THEN
                     ALFMAX = ALFA
                  ELSE
                     ALFBND = MIN(TWO*ALFA,ALFMAX)
                  END IF
               END IF
            END IF
C
            GTRY = GTRY - RMU*OLDG
C
         END IF
      END IF
C     +    until (      done)
      IF ( .NOT. DONE) GO TO 40
C
      NFUN = NFUN + NUMF
      IF ( .NOT. NEEDFD) NGRAD = NGRAD + NUMF
      ALFA = ALFBST
C
      IF ( .NOT. IMPRVD) THEN
         CALL DCOPY(N,X1,1,X,1)
         CALL DAXPY(N,ALFA,DX,1,X,1)
         IF (NCNLN.GT.0) THEN
            IF (ALFA.LE.ONE) THEN
               CALL DCOPY(NCNLN,CMUL1,1,CMUL,1)
               CALL DAXPY(NCNLN,ALFA,DLAM,1,CMUL,1)
            END IF
            CALL DCOPY(NCNLN,SLK1,1,SLK,1)
            CALL DAXPY(NCNLN,ALFA,DSLK,1,SLK,1)
            CALL DCOPY(NCNLN,C,1,CS,1)
            CALL DAXPY(NCNLN,(-ONE),SLK,1,CS,1)
         END IF
      END IF
C
      RETURN
C
C     The user wants to stop.  Who am I to object?
C
   60 INFORM = MODE
      RETURN
C
C
C     End of  E04UCR. (NPSRCH)
C
      END
