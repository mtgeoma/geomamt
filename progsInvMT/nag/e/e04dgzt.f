      SUBROUTINE E04DGZ(N,X,XNORM,FUNGRD,NSTATE,ITER,ITMAX,NFEVAL,BIGDX,
     *                  INFORM,FTOL,EPSRF,OBJF,IDBGCG,MSGCG,SLPRT,ITPRT,
     *                  ETA,DEBUG,FGUESS,SK,YK,DIAGB,OLDDB,SR,YR,OLDG,
     *                  HG,HYK,PK,HYR,G,WX,WGRAD,IUSER,USER)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     MARK 16A REVISED. IER-987 (SEP 1993).
C
C     Variables here such as SK and YK correspond to the same
C     variables in the CG paper.
C
C     EPSREF is a reference value which is used to update the value of
C     EPSAF later.
C
C     THETAQ is the theta( q ) on p.25 of the CG paper.
C
C     -- Written on 4-June-1986.
C     Sven Hammarling and Janet Welding, NAG Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGDX, EPSRF, ETA, FGUESS, FTOL, OBJF, XNORM
      INTEGER           IDBGCG, INFORM, ITER, ITMAX, MSGCG, N, NFEVAL,
     *                  NSTATE
      LOGICAL           DEBUG, ITPRT, SLPRT
C     .. Array Arguments ..
      DOUBLE PRECISION  DIAGB(N), G(N), HG(N), HYK(N), HYR(N), OLDDB(N),
     *                  OLDG(N), PK(N), SK(N), SR(N), USER(*), WGRAD(N),
     *                  WX(N), X(N), YK(N), YR(N)
      INTEGER           IUSER(*)
C     .. Subroutine Arguments ..
      EXTERNAL          FUNGRD
C     .. Scalars in Common ..
      INTEGER           IPRINT, ISUMM, LINES1, LINES2, NOUT
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALFA, ALFMAX, CUBERT, DENOM, DIFNEW, DIFOLD,
     *                  EPSAF, EPSMCH, FKEEP, GHYK, GHYR, GNORM, GSK,
     *                  GSR, GTG, GTPNEW, OLDF, OLDGTP, PE, PNORM,
     *                  RTFTOL, SPE, THETAQ, XKXKNM, YKHYK, YKHYR, YKSK,
     *                  YKSR, YRHYR, YRSR
      INTEGER           I, ICYCLE, MODE, NM1
      LOGICAL           CN, ERROR, FAIL, FONLY
C     .. Local Arrays ..
      CHARACTER*100     REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2, E04DGW, F06BLF
      LOGICAL           E04DGV
      EXTERNAL          DDOT, DNRM2, E04DGW, F06BLF, E04DGV
C     .. External Subroutines ..
      EXTERNAL          DCOPY, E04DGT, E04DGU, E04DGX, E04DGY, F06FBF, 
     *                  F06FGF, X04BAY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Common blocks ..
      COMMON            /AE04NB/NOUT, IPRINT, ISUMM, LINES1, LINES2
      COMMON            /AX02ZA/WMACH
C     .. Save statement ..
      SAVE              /AX02ZA/
C     .. Executable Statements ..
C
C     Initialization.
C
      FONLY = .FALSE.
      MODE = 2
      CALL FUNGRD(MODE,N,X,OBJF,G,NSTATE,IUSER,USER)
      NSTATE = 0
      NFEVAL = NFEVAL + 1
      IF (MODE.LT.0) THEN
         INFORM = MODE
         GO TO 140
      END IF
      GTG = DDOT(N,G,1,G,1)
      EPSMCH = WMACH(3)
      EPSAF = EPSRF*(ONE+ABS(OBJF))
      IF (FTOL.LT.EPSRF) THEN
         IF (MSGCG.GT.0) THEN
            WRITE (REC,FMT=99999) EPSRF
            CALL X04BAY(NOUT,2,REC)
         END IF
         FTOL = EPSRF
      END IF
      RTFTOL = SQRT(FTOL)
      CUBERT = FTOL**(1.0D0/3.0D0)
      ITER = 0
      NM1 = N - 1
      ICYCLE = NM1
      THETAQ = 5.0D-2
      IF (GTG.LT.EPSAF) THEN
         INFORM = 8
         GO TO 140
      END IF
      INFORM = 0
      CALL DCOPY(N,G,1,PK,1)
      GTPNEW = -GTG
      PNORM = SQRT(GTG)
      DIFNEW = ZERO
      FKEEP = OBJF
      GNORM = PNORM
      CALL F06FBF(N,ZERO,YR,1)
      CALL F06FBF(N,ZERO,SR,1)
C
C     ***** Start of main iteration loop *****
C
C     Exit if the maximum number of iterations is reached.
C
   20 CONTINUE
      IF (ITER.GE.ITMAX) THEN
         INFORM = 3
         GO TO 140
      END IF
      CALL F06FGF(N,PK,1)
      OLDF = OBJF
      OLDGTP = GTPNEW
      CALL DCOPY(N,G,1,OLDG,1)
      PE = PNORM + EPSMCH
      SPE = BIGDX/PE
      ALFA = E04DGW(OBJF,FGUESS,OLDGTP,SPE)
      IF (DEBUG .AND. ITER.GE.IDBGCG) THEN
         WRITE (REC,FMT=99998) ITER, ICYCLE, EPSAF, THETAQ, GTG, PNORM,
     *     GNORM
         CALL X04BAY(NOUT,3,REC)
      END IF
      IF (MSGCG.GE.5 .AND. ITER.EQ.0) THEN
C
C        Print the details of initial values.
C
         ITPRT = .TRUE.
         CALL E04DGT(MSGCG,ITPRT,SLPRT,ITER,NFEVAL,ALFA,OBJF,GNORM,G,X,
     *               XNORM,XKXKNM,N)
      END IF
C
C     Perform the line search.
C
      XNORM = DNRM2(N,X,1)
      PNORM = DNRM2(N,PK,1)
C
      ALFMAX = F06BLF(BIGDX,PNORM,FAIL)
      CALL DCOPY(N,X,1,WX,1)
      CALL E04DGU(DEBUG,IDBGCG,ITER,FONLY,FUNGRD,N,ALFMAX,EPSRF,ETA,
     *            PNORM,XNORM,PK,ALFA,OBJF,GTPNEW,G,X,INFORM,NFEVAL,
     *            WGRAD,WX,IUSER,USER)
      ITER = ITER + 1
C
C     Exit if the linesearch routine failed to find a better point.
C
      ERROR = (INFORM.EQ.7) .OR. (INFORM.EQ.8)
      IF (ERROR) INFORM = 6
      IF (INFORM.LT.0 .OR. INFORM.GT.3) GO TO 140
C
      EPSAF = EPSRF*(ONE+ABS(OBJF))
      GTG = DDOT(N,G,1,G,1)
      GNORM = SQRT(GTG)
      XNORM = DNRM2(N,X,1)
      IF (MSGCG.GE.5) THEN
C
C        Print details of this iteration.
C
         XKXKNM = ALFA*PNORM
         ITPRT = .TRUE.
         CALL E04DGT(MSGCG,ITPRT,SLPRT,ITER,NFEVAL,ALFA,OBJF,GNORM,G,X,
     *               XNORM,XKXKNM,N)
      END IF
C
C     Test for convergence.
C
      INFORM = 0
      CN = E04DGV(OBJF,OLDF,ALFA,PNORM,RTFTOL,CUBERT,FTOL,XNORM,GNORM,
     *     EPSAF)
C
C     Exit from the main iteration if algorithm has converged.
C
      IF (CN) GO TO 140
C
      DIFOLD = DIFNEW
      DIFNEW = OLDF - OBJF
C
C     If this is the first iteration if a new cycle compute the
C     percentage reduction factor for the resetting test.
C     The formula for the new value of theta(q) is given on page
C     25 of the CG paper.
C
      IF (ICYCLE.EQ.1) THEN
         IF (DIFNEW.GT.(2*DIFOLD)) THETAQ = 2*THETAQ
         IF (DIFNEW.LT.(DIFOLD/2)) THETAQ = THETAQ/2
      END IF
C
C     Compute the change in the iterates and the corresponding change
C     in the gradients. Update diagonals of approximate Hessian.
C
      DO 40 I = 1, N
         YK(I) = G(I) - OLDG(I)
         SK(I) = PK(I)*ALFA
   40 CONTINUE
      CALL DCOPY(N,DIAGB,1,OLDDB,1)
      CALL E04DGY(N,DIAGB,OLDDB,ALFA,PK,YK,OLDG,ITER,DEBUG,IDBGCG)
C
C     Test for start of new cycle, then compute new search
C     direction.  The test is that of ( 6.4 ), p.25 of CG paper.
C
C     The new search direction is calculated from ( 3.8 ), p.16 of
C     the CG paper. The u( 1 ) there corresponds to the inverse of
C     the diagonals here, y(t) and s(t) there are accumulated
C     gradient difference and step, that is, YR and SR here.
C
C     The approximate inverse Hessian is never actually computed.
C     Instead, routine E04DGX will return the product of the
C     updated Hessian with an arbitrary vector. It therefore
C     requires as input the product of the vector to be multiplied
C     with other vectors or matrices. The formula being used
C     is that for h( k + 1 ) on page 16 of the CG paper, with
C     gamma equal to one. By multiplying this expression on
C     the right by the vector v, say one can tell what input
C     the E04DGX routine will need.
C
C     In this particular case, we want to start with the inverse
C     of the diagonals and perform two updates ( as in formula
C     ( 3.8 ) ), and compute the new search direction as the
C     product of the new matrix with the gradient -- this direction
C     is negated in the call to F06FGF at the start of the main
C     iteration loop. Basically, this
C     multiplication is carried out in the third call to E04DGX
C     below. The first two calls to E04DGX compute vectors
C     needed as input to the third call. Thus, for example
C     in terms of the notation on page 16 of the paper we wish
C     to compute h( k + 1 )*g. h( k ) in that formula will
C     be u( 2 ), i.e. the result of updating the inverse of
C     the diagonals. Therefore we must compute u( 2 )*g
C     for input to this third call to E04DGX, and this is
C     done in the first call to E04DGX -- the code below calls
C     this product HG. Similarly, the second call to E04DGX
C     computes u( 2 )*yk, called HYK below.
C
      YKSK = DDOT(N,YK,1,SK,1)
      GSK = DDOT(N,G,1,SK,1)
      YRSR = DDOT(N,YR,1,SR,1)
      IF ((ICYCLE.NE.NM1) .AND. (DIFNEW.GE.(THETAQ*(FKEEP-OBJF)))
     *    .AND. (YRSR.GT.ZERO)) THEN
C
         DO 60 I = 1, N
            DENOM = DIAGB(I)
            HG(I) = G(I)/DENOM
            HYK(I) = YK(I)/DENOM
            HYR(I) = YR(I)/DENOM
   60    CONTINUE
         YKSR = DDOT(N,YK,1,SR,1)
         YKHYR = DDOT(N,YK,1,HYR,1)
         GSR = DDOT(N,G,1,SR,1)
         GHYR = DDOT(N,G,1,HYR,1)
         YRHYR = DDOT(N,YR,1,HYR,1)
C
         CALL E04DGX(N,ONE,SR,YR,HG,HYR,YRSR,YRHYR,GSR,GHYR,HG,ITER,
     *               DEBUG,IDBGCG)
         CALL E04DGX(N,ONE,SR,YR,HYK,HYR,YRSR,YRHYR,YKSR,YKHYR,HYK,ITER,
     *               DEBUG,IDBGCG)
C
         YKHYK = DDOT(N,HYK,1,YK,1)
         GHYK = DDOT(N,HYK,1,G,1)
         CALL E04DGX(N,ONE,SK,YK,HG,HYK,YKSK,YKHYK,GSK,GHYK,PK,ITER,
     *               DEBUG,IDBGCG)
         GTPNEW = -DDOT(N,G,1,PK,1)
         PNORM = DNRM2(N,PK,1)
C
C        Compute the accumulated step and its corresponding
C        gradient difference.
C
         DO 80 I = 1, N
            SR(I) = SR(I) + SK(I)
            YR(I) = YR(I) + YK(I)
   80    CONTINUE
         ICYCLE = ICYCLE + 1
      ELSE
C
C        A new cycle is starting, per test on page 25 of CG paper.
C        Now compute the new search direction, pk. This requires
C        computing hg and hy, which are needed as input to the
C        E04DGX routine. There is no accumulated step at this
C        point, so only one call to E04DGX is needed. That is,
C        the inverse of diagonals is used to produce a BFGS updated
C        matrix, and new pk = this matrix * g( which will be negated
C        to give descent direction at start of main iteration loop ).
C
         IF (DEBUG .AND. ITER.GE.IDBGCG) THEN
            WRITE (REC,FMT=99997)
            CALL X04BAY(NOUT,2,REC)
         END IF
         DO 100 I = 1, N
            DENOM = DIAGB(I)
            HG(I) = G(I)/DENOM
            HYK(I) = YK(I)/DENOM
  100    CONTINUE
         YKHYK = DDOT(N,YK,1,HYK,1)
         GHYK = DDOT(N,G,1,HYK,1)
         CALL E04DGX(N,ONE,SK,YK,HG,HYK,YKSK,YKHYK,GSK,GHYK,PK,ITER,
     *               DEBUG,IDBGCG)
         GTPNEW = -DDOT(N,G,1,PK,1)
         PNORM = DNRM2(N,PK,1)
         GNORM = DNRM2(N,G,1)
C
C        Initialize the sum of all the changes in X.
C
         DO 120 I = 1, N
            SR(I) = SK(I)
            YR(I) = YK(I)
  120    CONTINUE
         FKEEP = OBJF
         ICYCLE = 1
      END IF
      GO TO 20
C
C     ***** End of main iteration *****
C
  140 CONTINUE
C
      IF (INFORM.GE.0) THEN
         IF (MSGCG.EQ.1 .OR. MSGCG.EQ.10) THEN
C
C           Print the solution.
C
            SLPRT = .TRUE.
            CALL E04DGT(MSGCG,ITPRT,SLPRT,ITER,NFEVAL,ALFA,OBJF,GNORM,G,
     *                  X,XNORM,XKXKNM,N)
         END IF
      END IF
      RETURN
C
C
C     End of E04DGZ (CGCORE).
C
99999 FORMAT (/' XXX  Warning: FTOL increased to ',G16.8)
99998 FORMAT (/' //E04DGZ// - ITER =',I5,' ICYCLE =',I3,' EPSAF =',
     *       G16.8,' THETAQ =',G16.8,/'              GTG  =',G16.8,' P',
     *       'NORM =',G16.8,' GNORM =',G16.8)
99997 FORMAT (/' //E04DGZ// - A new cycle is starting ..')
      END
