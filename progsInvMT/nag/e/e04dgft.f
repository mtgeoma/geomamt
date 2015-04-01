      SUBROUTINE E04DGF(N,FUNGRD,ITER,OBJF,OBJGRD,X,IWORK,WORK,IUSER,
     *                  USER,IFAIL)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     MARK 16 REVISED. IER-1053 (JUL 1993).
C     MARK 17 REVISED. IER-1555 (JUN 1995).
C
C     ==================================================================
C     E04DGF is a pre-conditioned, limited memory quasi-Newton method,
C     based on method PLMA of "Conjugate gradient methods for
C     large-scale nonlinear optimization", Technical Report SOL 79-15.
C     ==================================================================
C
C     -----------
C     Parameters
C     -----------
C
C     N      - INTEGER.
C           On entry, N must contain the number of variables in the
C           problem.
C
C     FUNGRD - SUBROUTINE, supplied by the user.
C           FUNGRD must be declared as EXTERNAL in the routine that
C           calls E04DGF.
C           If MODE = 0 then FUNGRD must calculate the objective
C           function, if MODE = 2 it must also calculate the gradient
C           vector.
C
C          It's specification is:
C          SUBROUTINE FUNGRD( MODE, N, X, OBJF, G, NSTATE, IUSER, USER )
C          INTEGER    MODE, N, NSTATE, IUSER(*)
C          REAL       X(N), OBJF, G(N), USER(*)
C
C          MODE   - INTEGER.
C          MODE indicates which parameter values within FUNGRD need to
C          be set.
C          If MODE = 0, FUNGRD must only return the value of the
C          objective function (OBJF) at X.
C          If MODE = 2, FUNGRD must also return the gradient ( G(1)...
C          ...G(N) ).
C          If MODE is negative on exit from FUNGRD, the execution of
C          E04DGF is terminated with IFAIL set to MODE.
C
C          N      - INTEGER.
C          On entry, N specifies the number of variables as input to
C          E04DGF.
C          N must not be altered by FUNGRD.
C
C          X      - REAL array of dimension (N).
C          On entry, X contains the point at which the objective
C          function is to be evaluated.
C          X must not be altered by FUNGRD.
C
C          OBJF   - REAL.
C          On exit, OBJF must contain the value of the objective
C          function.
C
C          G      - REAL array of dimension (N).
C          On exit, if MODE = 2 G(1) ... G(N) must contain the gradient
C          vector of the objective function. The j-th component of G
C          must contain the partial derivative of the function with
C          respect to the j-th variable.
C
C          NSTATE - INTEGER.
C          On entry, NSTATE will be set to 1 on the first call of FUNGRD
C          by E04DGF, and is set to 0 for all subsequent calls. Thus, if
C          the user wishes, NSTATE may be tested within FUNGRD in order
C          to perform cerain calculations once only. For example, the
C          user may read data or initialise COMMON blocks when NSTATE
C          = 1.
C
C          IUSER  - INTEGER array of dimension as least (1).
C          USER   - REAL    array of dimension at least (1).
C          FUNGRD is called from E04DGF with the parameters IUSER and
C          USER as supplied to E04DGF. The user is free to use arrays
C          IUSER and USER to supply information to FUNGRD as an
C          alternative to using COMMON.
C
C     ITER   - INTEGER.
C           On exit, ITER is set to the number of iterations performed
C           by E04DGF.
C
C     OBJF   - REAL.
C           On exit, contains the value of the objective function at the
C           final iterate.
C
C     OBJGRD - REAL array of length at least (N).
C           On exit, OBJGRD contains the gradient vector at the final
C           iterate.
C
C     X      - REAL array of length at least (N).
C           On entry, X must contain the initial estimate of the
C           solution.
C           On exit, X contains the final estimate of the solution.
C
C     IWORK  - INTEGER array of dimension at least (N+1).
C           Work array.
C
C     WORK   - REAL array of dimension at least (13*N).
C           Work array.
C
C     IFAIL  - INTEGER.
C           On entry, IFAIL must be set to 0,-1 or 1. The recommended
C           value for IFAIL is -1.
C           On exit, IFAIL represents a diagnostic indicator. The values
C           of IFAIL on exit are as follows:-
C
C           IFAIL < 0  the user requested termination by setting MODE
C                      negative within routine FUNGRD.
C
C           IFAIL = 0  indicates successful termination.
C
C           IFAIL = 3  maximum number of function evaluations have
C                      been performed.
C
C           IFAIL = 4  the computed upper bound on the step length
C                      taken during the linesearch was too small.
C                      A re-run with a larger value assigned to the
C                      maximum step length ( i.e. BIGDX ) may be
C                      successful. If BIGDX is already large ( .ge.
C                      the default value ) then current point cannot
C                      be improved upon.
C
C           IFAIL = 6  the current point cannot be improved upon.
C                      A sufficient decrease in the function value
C                      could not be attained during the final
C                      linesearch. If the subroutine OBJFUN computes
C                      the function and gradients correctly, then
C                      this may occur because an overly stringent
C                      accuracy has been requested i.e. Optimality
C                      tolerance (FTOL) is too small or if the
C                      minimum lies close to a step length of zero.
C
C           IFAIL = 7  large errors were found in the derivatives of
C                      the objective function. This value of IFAIL
C                      will occur if the verification process
C                      indicated that at least one gradient component
C                      had no correct figures. The user should refer
C                      to the printed output to determine which
C                      elements are suspect to be in error.
C
C           IFAIL = 8  indicates that the gradient (g) at the starting
C                      point is too small. The value g(trans)*g is
C                      less than machine precision. The problem
C                      should be rerun at a different starting point.
C
C           IFAIL = 9  an input parameter is invalid.
C
C
C     -------------------
C     Optional Parameters
C     -------------------
C
C     ITMAX   an INTEGER, the maximum number of iterations to be
C          performed. If ITMAX .lt. 0, then a default value of
C          max(50,5*N) is taken.
C
C     BIGDX   a REAL variable, the maximum allowable step length. The
C          default value is 1.0E+20.
C
C     FGUESS  a REAL variable, user supplied guess of optimum objective
C          function value. FGUESS is F(est) of CG paper (p.10 ).
C
C     EPSRF   a REAL variable which specifies the relative precision
C          of the objective function at the starting point.
C
C     ETA     is a REAL variable used to specify the accuracy of the
C          linesearch.
C
C     FTOL    is a REAL variable which indicates the correct number of
C          figures the user desires in the optimum function value.
C          For this purpose, leading zeros after the decimal are
C          considered as candidates for correct figures.
C
C     MSGCG   is an INTEGER variable which indicates the amount of
C          output required by the user.
C          The values of MSGCG are as follows:-
C          MSGCG = 0   No printout.
C          MSGCG = 1   Final solution only.
C          MSGCG = 5   One line of output for each iteration.
C          MSGCG = 10  Final solution and brief line of output for each
C                      iteration.
C
C     LDBGCG  is an INTEGER variable controlling debug output. If LDBGCG
C          > 0, then debug printout will be produced.
C
C     IDBGCG  is an INTEGER variable which indicates the number of the
C          iteration from which debug print should start.
C
C     ******************************************************************
C
C     -- Written on 4th June 1986.
C     Sven Hammarling and Janet Welding, NAG Central Office.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E04DGF')
      INTEGER           LDBG
      PARAMETER         (LDBG=5)
      INTEGER           MXPARM, NIPARM, NRPARM
      PARAMETER         (MXPARM=30,NIPARM=8,NRPARM=5)
      DOUBLE PRECISION  BIGBND
      PARAMETER         (BIGBND=1.0D+20)
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
      DOUBLE PRECISION  POINT3, POINT8, POINT9
      PARAMETER         (POINT3=3.3D-1,POINT8=0.8D+0,POINT9=0.9D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  OBJF
      INTEGER           IFAIL, ITER, N
C     .. Array Arguments ..
      DOUBLE PRECISION  OBJGRD(N), USER(*), WORK(13*N), X(N)
      INTEGER           IUSER(*), IWORK(N+1)
C     .. Subroutine Arguments ..
      EXTERNAL          FUNGRD
C     .. Scalars in Common ..
      DOUBLE PRECISION  BIGDX, EPSPT3, EPSPT5, EPSPT8, EPSPT9, EPSRF,
     *                  ETA, FGUESS, FTOL
      INTEGER           IDBGCG, IPRINT, ISUMM, ITMAX, JVRFY1, JVRFY2,
     *                  LDBGCG, LFDSET, LINES1, LINES2, LVERFY, LVLDIF,
     *                  LVRFYC, MSGCG, NCDIFF, NFDIFF, NN, NOUT
      LOGICAL           DEBUG
C     .. Arrays in Common ..
      DOUBLE PRECISION  RPADCG(MXPARM-NRPARM), RPSVCG(MXPARM), WMACH(15)
      INTEGER           INPDBG(LDBG), IPADCG(MXPARM-NIPARM),
     *                  IPSVCG(MXPARM), JVERFY(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  CDINT, EPSMCH, FDCHK, FDINT, FDNORM, RTEPS,
     *                  XNORM
      INTEGER           I, INFORM, K, LDIAGB, LHG, LHYK, LHYR, LOLDDB,
     *                  LOLDG, LPK, LSK, LSR, LVLDER, LWGRAD, LWX, LYK,
     *                  LYR, MSG1, NCNLN, NERR, NFEVAL, NFUN, NGRAD,
     *                  NROWJ, NROWUJ, NSTATE
      LOGICAL           ITPRT, SLPRT
      CHARACTER*12      TITLE
C     .. Local Arrays ..
      DOUBLE PRECISION  CDUM(1,1), RPRMCG(MXPARM), UDUM(1,1), W(1)
      INTEGER           IPRMCG(MXPARM)
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2
      INTEGER           P01ABF
      EXTERNAL          DNRM2, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, E04DGR, E04DGZ, E04UCY, E04UDM, F06DFF,
     *                  F06FBF, X02ZAZ, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         MOD, SQRT
C     .. Common blocks ..
      COMMON            /AE04DG/IPSVCG, IDBGCG, ITMAX, JVRFY1, JVRFY2,
     *                  LDBGCG, LVERFY, MSGCG, NN, IPADCG
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
      COMMON            /BE04DG/RPSVCG, BIGDX, EPSRF, ETA, FGUESS, FTOL,
     *                  RPADCG
      COMMON            /BE04UC/LVLDIF, NCDIFF, NFDIFF, LFDSET
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
      COMMON            /CE04UC/LVRFYC, JVERFY
      COMMON            /FE04UC/INPDBG, DEBUG
C     .. Equivalences ..
      EQUIVALENCE       (IPRMCG(1),IDBGCG), (RPRMCG(1),BIGDX)
C     .. Save statement ..
      SAVE              /AE04DG/, /BE04DG/, /AX02ZA/
C     .. Data statements ..
      DATA              TITLE/' *** E04DGF '/
C     .. Executable Statements ..
C
C     Initialise variables.
C
      CALL X02ZAZ
      NOUT = WMACH(11)
      NERR = WMACH(12)
C
      EPSMCH = WMACH(3)
      RTEPS = WMACH(4)
      EPSPT3 = EPSMCH**POINT3
      EPSPT5 = RTEPS
      EPSPT8 = EPSMCH**POINT8
      EPSPT9 = EPSMCH**POINT9
C
C     Test input parameters.
C
      IF (N.LE.0) THEN
         IF (IFAIL.EQ.0 .OR. IFAIL.EQ.1) THEN
            WRITE (REC,FMT=99990) N
            CALL X04BAY(NERR,2,REC)
         END IF
         INFORM = 9
         GO TO 60
      END IF
C
      DEBUG = LDBGCG .GT. 0
      INFORM = 0
      NSTATE = 1
      NFEVAL = 0
      SLPRT = .FALSE.
      ITPRT = .FALSE.
C
C     Set the default values for the parameters.
C
      CALL E04DGR(N,TITLE)
C
      XNORM = DNRM2(N,X,1)
      IF (LVERFY.GE.0) THEN
C
C        Set variables for gradient check by E04UCY.
C
         NFUN = 0
         NGRAD = 0
         NROWJ = 1
         NROWUJ = 1
         NCNLN = 0
C        ***    The following variables should be set in the
C        defaults routine if finite-difference approximations to
C        the gradient are incorporated.
         JVERFY(1) = JVRFY1
         JVERFY(2) = JVRFY2
         JVERFY(3) = 1
         JVERFY(4) = N
         LVLDER = 3
         LFDSET = 0
         LVLDIF = 0
         K = 1
         MSG1 = IDBGCG
         DO 20 I = 1, LDBG
            INPDBG(I) = MOD(MSG1/K,10)
            K = K*10
   20    CONTINUE
C
C        ***     End of default variables.
C
         IF (LFDSET.EQ.0) FDCHK = SQRT(EPSRF)
C
C        Set up dummy arrays bu, bl in work array.
C
         DO 40 I = 1, N
            WORK(I) = -BIGBND
            WORK(N+I) = BIGBND
   40    CONTINUE
         LVRFYC = LVERFY
C
C        Check gradients
C
         CALL E04UCY(INFORM,MSGCG,NSTATE,LVLDER,NFUN,NGRAD,NROWJ,NROWUJ,
     *               N,NCNLN,E04UDM,FUNGRD,IWORK,BIGBND,EPSRF,CDINT,
     *               FDINT,FDCHK,FDNORM,OBJF,XNORM,WORK,WORK(N+1),
     *               WORK(6*N+1),WORK(6*N+2),CDUM,UDUM,WORK(6*N+3),
     *               WORK(2*N+1),WORK(3*N+1),OBJGRD,W,W,X,WORK(4*N+1),
     *               WORK(5*N+1),IUSER,USER)
         NSTATE = 0
         IF (INFORM.NE.0) THEN
            IF (INFORM.GT.0) INFORM = 7
            GO TO 60
         END IF
      END IF
C
C     Set the work arrays for E04DGZ.
C
      LYK = 1
      LDIAGB = LYK + N
      LSR = LDIAGB + N
      LYR = LSR + N
      LOLDG = LYR + N
      LHG = LOLDG + N
      LHYK = LHG + N
      LPK = LHYK + N
      LHYR = LPK + N
      LSK = LHYR + N
      LWX = LSK + N
      LWGRAD = LWX + N
      LOLDDB = LWGRAD + N
C
      CALL F06FBF(N,ONE,WORK(LDIAGB),1)
C
      CALL E04DGZ(N,X,XNORM,FUNGRD,NSTATE,ITER,ITMAX,NFEVAL,BIGDX,
     *            INFORM,FTOL,EPSRF,OBJF,IDBGCG,MSGCG,SLPRT,ITPRT,ETA,
     *            DEBUG,FGUESS,WORK(LSK),WORK(LYK),WORK(LDIAGB),
     *            WORK(LOLDDB),WORK(LSR),WORK(LYR),WORK(LOLDG),WORK(LHG)
     *            ,WORK(LHYK),WORK(LPK),WORK(LHYR),OBJGRD,WORK(LWX),
     *            WORK(LWGRAD),IUSER,USER)
C
C     Print messages if required.
C
   60 IF (MSGCG.GT.0) THEN
         IF (INFORM.LT.0) WRITE (REC,FMT=99999)
         IF (INFORM.EQ.0) WRITE (REC,FMT=99998)
         IF (INFORM.EQ.3) WRITE (REC,FMT=99997)
         IF (INFORM.EQ.4) WRITE (REC,FMT=99996)
         IF (INFORM.EQ.6) WRITE (REC,FMT=99995)
         IF (INFORM.EQ.7) WRITE (REC,FMT=99991)
         IF (INFORM.EQ.8) WRITE (REC,FMT=99994)
         IF (INFORM.EQ.9) WRITE (REC,FMT=99993)
         CALL X04BAY(NOUT,2,REC)
C
         IF (INFORM.GE.0 .AND. INFORM.LE.6) THEN
            WRITE (REC,FMT=99992) OBJF
            CALL X04BAY(NOUT,2,REC)
         END IF
      END IF
C
C     Recover the optional parameters set by the user.
C
      CALL F06DFF(MXPARM,IPSVCG,1,IPRMCG,1)
      CALL DCOPY(MXPARM,RPSVCG,1,RPRMCG,1)
C
      IF (INFORM.NE.0 .AND. (IFAIL.EQ.0 .OR. IFAIL.EQ.-1)) THEN
         IF (INFORM.LT.0) WRITE (REC,FMT=99989)
         IF (INFORM.EQ.3) WRITE (REC,FMT=99988)
         IF (INFORM.EQ.4) WRITE (REC,FMT=99987)
         IF (INFORM.EQ.6) WRITE (REC,FMT=99986)
         IF (INFORM.EQ.7) WRITE (REC,FMT=99985)
         IF (INFORM.EQ.8) WRITE (REC,FMT=99984)
         IF (INFORM.EQ.9) WRITE (REC,FMT=99983)
         CALL X04BAY(NERR,2,REC)
      END IF
      IFAIL = P01ABF(IFAIL,INFORM,SRNAME,0,REC)
C
      RETURN
C
C     End of E04DGF (CGSOL).
C
99999 FORMAT (/' Exit E04DGF - User requested termination.')
99998 FORMAT (/' Exit E04DGF - Optimal solution found.')
99997 FORMAT (/' Exit E04DGF - Too many iterations.')
99996 FORMAT (/' Exit E04DGF - Computed upper bound on step length is ',
     *       'too small.')
99995 FORMAT (/' Exit E04DGF - Current point cannot be improved upon.')
99994 FORMAT (/' Exit E04DGF - Gradient at the starting point is too s',
     *       'mall.')
99993 FORMAT (/' Exit E04DGF - 1 error found in the input parameters. ',
     *       ' Problem abandoned.')
99992 FORMAT (/' Final objective value =',G16.7)
99991 FORMAT (/' Exit E04DGF - Large errors found in the derivatives.')
99990 FORMAT (/' ** On entry, N.le.0:',/'    N = ',I16)
99989 FORMAT (/' ** User requested termination.')
99988 FORMAT (/' ** Too many iterations.')
99987 FORMAT (/' ** Computed upper bound on step length is too small.')
99986 FORMAT (/' ** Current point cannot be improved upon.')
99985 FORMAT (/' ** Large errors found in the derivatives.')
99984 FORMAT (/' ** Gradient at the starting point is too small.')
99983 FORMAT (/' ** 1 error found in the input parameters.  Problem ab',
     *       'andoned.')
      END
