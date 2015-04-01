      SUBROUTINE E04DGU(DEBUG,IDBGCG,ITER,NEEDFD,OBJFUN,N,ALFMAX,EPSRF,
     *                  ETA,DXNORM,XNORM,DX,ALFA,OBJF,GDX,GRAD,X,INFORM,
     *                  NFUN,UGRAD,X1,IUSER,USER)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     MARK 16 REVISED. IER-1056 (JUL 1993).
C
C     E04DGU is the line search routine for E04DGF.
C
C     -- Original version (Mark 12) written on 4-June-1986.
C     Sven Hammarling and Janet Welding, NAG Central Office.
C
C     -- This version (Mark 16) written on 15-January-1993.
C     Alan Brown, NAG Ltd.
C
C     .. Parameters ..
      DOUBLE PRECISION  RMU
      PARAMETER         (RMU=1.0D-4)
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALFA, ALFMAX, DXNORM, EPSRF, ETA, GDX, OBJF,
     *                  XNORM
      INTEGER           IDBGCG, INFORM, ITER, N, NFUN
      LOGICAL           DEBUG, NEEDFD
C     .. Array Arguments ..
      DOUBLE PRECISION  DX(N), GRAD(N), UGRAD(N), USER(*), X(N), X1(N)
      INTEGER           IUSER(*)
C     .. Subroutine Arguments ..
      EXTERNAL          OBJFUN
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPSPT3, EPSPT5, EPSPT8, EPSPT9
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALFBST, ALFSML, EPSAF, EPSMCH, FBEST, FTRY, G0,
     *                  GBEST, GTRY, OLDF, OLDG, Q, S, T, TARGTG, TOBJ,
     *                  TOBJM, TOLABS, TOLAX, TOLREL, TOLRX, TOLTNY
      INTEGER           J, MAXF, MODE, NSTATE, NUMF
      LOGICAL           DBGSR, DONE, FIRST, IMPRVD
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, E04UCJ, E04UCK
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
C     .. Save statement ..
      SAVE              /AX02ZA/
C     .. Executable Statements ..
C
      EPSMCH = WMACH(3)
C
      NSTATE = 0
C
      IF (NEEDFD) THEN
         MAXF = 15
      ELSE
         MAXF = 10
      END IF
C
      EPSAF = EPSRF*(ONE+ABS(OBJF))
      TOLRX = EPSPT8
      TOLAX = EPSPT8
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
      OLDF = OBJF
      OLDG = GDX
      ALFBST = ZERO
      FBEST = ZERO
      GBEST = (ONE-RMU)*OLDG
      TARGTG = (RMU-ETA)*OLDG
      G0 = GBEST
C
      IF (DEBUG .AND. ITER.GE.IDBGCG) THEN
         DBGSR = .TRUE.
      ELSE
         DBGSR = .FALSE.
      END IF
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
C     FIRST = true for the first entry,  false for subsequent entries.
C     DONE  = true indicates termination, in which case the value of
C     INFORM gives the result of the search.
C     ------------------------------------------------------------------
C     +    REPEAT
   40 IF (NEEDFD) THEN
         CALL E04UCJ(FIRST,DBGSR,DONE,IMPRVD,INFORM,MAXF,NUMF,IPRINT,
     *               ALFMAX,ALFSML,EPSAF,G0,TARGTG,FTRY,TOLABS,TOLREL,
     *               TOLTNY,ALFA,ALFBST,FBEST)
      ELSE
         CALL E04UCK(FIRST,DBGSR,DONE,IMPRVD,INFORM,MAXF,NUMF,IPRINT,
     *               ALFMAX,EPSAF,G0,TARGTG,FTRY,GTRY,TOLABS,TOLREL,
     *               TOLTNY,ALFA,ALFBST,FBEST,GBEST)
      END IF
C
      IF (IMPRVD) THEN
         OBJF = TOBJ
C
         IF ( .NOT. NEEDFD) CALL DCOPY(N,UGRAD,1,GRAD,1)
C
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
C        ------------------------------------------------------------
C        Compute the value and gradient of the objective function.
C        ------------------------------------------------------------
         CALL OBJFUN(MODE,N,X,TOBJ,UGRAD,NSTATE,IUSER,USER)
         IF (MODE.LT.0) GO TO 60
C
         TOBJM = TOBJ
C
         FTRY = TOBJM - OLDF - RMU*OLDG*ALFA
C
C        ------------------------------------------------------------
C        Compute auxiliary gradient information.
C        ------------------------------------------------------------
         IF ( .NOT. NEEDFD) GTRY = DDOT(N,UGRAD,1,DX,1)
      END IF
C     +    UNTIL (      DONE)
      IF ( .NOT. DONE) GO TO 40
C
      NFUN = NFUN + NUMF
      ALFA = ALFBST
C
      IF ( .NOT. IMPRVD) THEN
         CALL DCOPY(N,X1,1,X,1)
         CALL DAXPY(N,ALFA,DX,1,X,1)
      END IF
C
      RETURN
C
C     The user wants to stop.
C
   60 INFORM = MODE
      RETURN
C
C     End of E04DGU (SEARCH).
C
      END
