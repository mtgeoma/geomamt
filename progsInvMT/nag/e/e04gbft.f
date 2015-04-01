      SUBROUTINE E04GBF(M,N,LSQLIN,LSFJC,LSMON,IPRINT,MAXCAL,ETA,XTOL,
     *                  STEPMX,X,FSUMSQ,FVEC,FJAC,LJ,S,VT,LVT,NITER,
     *                  NFTOTL,IW,LIW,W,LW,IFAIL)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     MARK 16 REVISED. IER-989 (JUN 1993).
C
C     **************************************************************
C
C     E04GBF IS A COMPREHENSIVE QUASI-NEWTON ALGORITHM FOR FINDING
C     AN UNCONSTRAINED MINIMUM OF A SUM OF SQUARES OF M NONLINEAR
C     FUNCTIONS IN N VARIABLES (M .GE. N).  FIRST DERIVATIVES ARE
C     REQUIRED.
C
C     THE ROUTINE IS ESSENTIALLY IDENTICAL TO THE SUBROUTINE LSFDQ
C     IN THE NPL ALGORITHMS LIBRARY (REF. NO. E4/21/F). E04GBF
C     WILL NORMALLY BE CALLED WITH E04FCV (WHICH CALLS E04ABZ)
C     OR E04HEV (WHICH CALLS E04BBZ) AS THE PARAMETER LSQLIN.
C
C     PHILIP E. GILL, SUSAN M. PICKEN, WALTER MURRAY
C     BRIAN T. HINDE AND NICHOLAS I. M. GOULD
C     D.N.A.C., NATIONAL PHYSICAL LABORATORY, ENGLAND.
C
C     **************************************************************
C
C     Modified to output explanatory messages.
C     Peter Mayes, NAG Central Office, December 1987.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E04GBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ETA, FSUMSQ, STEPMX, XTOL
      INTEGER           IFAIL, IPRINT, LIW, LJ, LVT, LW, M, MAXCAL, N,
     *                  NFTOTL, NITER
C     .. Array Arguments ..
      DOUBLE PRECISION  FJAC(LJ,N), FVEC(M), S(N), VT(LVT,N), W(LW),
     *                  X(N)
      INTEGER           IW(LIW)
C     .. Subroutine Arguments ..
      EXTERNAL          LSFJC, LSMON, LSQLIN
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA, EPSEPS, EPSMCH, GTG, PEPS, PNORM, RTEPS,
     *                  RTOL, SIGMA, SSQNEW, SSQOLD, TAU, TLEPSQ, TOL,
     *                  TOLEPS, U, XNORM
      INTEGER           IFLAG, IGRADE, LEND, LGSSQ, LH, LHK, LP, LP2,
     *                  LPH, LPHESD, LPHESL, LRHS, LUTF, LW1, LW2,
     *                  MAXRNK, NPHI, NREC, NS, NWHY, NWY
      LOGICAL           GAUSS, NOMOVE
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2
      INTEGER           P01ABF
      EXTERNAL          DDOT, DNRM2, P01ABF
C     .. External Subroutines ..
      EXTERNAL          E04GBX, E04GBZ, E04GDQ, E04GDR, E04GDS, E04GDV,
     *                  E04GDX, E04GDY, E04GDZ, E04HEY, F04AQZ, F04JAZ,
     *                  F06FBF, DCOPY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
C
C     SET ADDRESSES FOR REAL WORKSPACE.
C
C     **************************************************************
C
C     THE FOLLOWING ADDRESSES AND ARRAY LENGTHS DEFINE BLOCKS OF THE
C     WORKSPACE ARRAY W(LW).
C
C     W(1), LENGTH 2*N + M + M*N   LOCAL WORKSPACE FOR AUXILIARY
C                                  SUBROUTINES.
C
C     W(LRHS), LENGTH N            THE VECTOR NEEDED TO CALCULATE
C                                  P2, THE CORRECTION TO THE
C                                  GRADED GAUSS-NEWTON DIRECTION.
C
C     W(LPHESL), LENGTH N*(N - 1)/2 THE LOWER TRIANGULAR MATRIX
C                                  ASSOCIATED WITH THE CHOLESKY
C                                  FACTORIZATION USED TO CALCULATE
C                                  P2.
C
C     W(LPHESD), LENGTH N          THE DIAGONAL MATRIX ASSOCIATED
C                                  WITH THE CHOLESKY FACTORIZATION
C                                  USED TO CALCULATE P2.
C
C     W(LP), LENGTH N              THE SEARCH DIRECTION VECTOR P.
C
C     W(LP2), LENGTH N             THE PROJECTION OF P2.
C
C     W(LGSSQ), LENGTH N           THE GRADIENT VECTOR OF THE SUM
C                                  OF SQUARES.
C
C     W(LUTF), LENGTH M            THE VECTOR OBTAINED FROM
C     .                            U TRANSPOSE * VECTOR OF FUNCTION
C     .                            VALUES.  U IS ASSOCIATED WITH
C     .                            THE SINGULAR VALUE DECOMPOSITION
C     .                                         T
C     .                            J = U * S * V .
C
C     W(LHK), LENGTH N*(N + 1)/2   THE SECOND TERM OF THE HESSIAN
C                                  MATRIX STORED BY ROWS AS A LOWER
C                                  TRIANGLE.
C
C     THE TOTAL LENGTH OF W SHOULD NOT BE LESS THAN
C     N*N + 7*N + N*M + 2*M OR 9 + 3*M IF N IS 1.
C     IW SHOULD HAVE AT LEAST ONE ELEMENT.
C
C     **************************************************************
C
      LH = N*(N+1)/2
      LRHS = 2*N + M + M*N + 1
      LPHESL = LRHS + N
      LPHESD = LPHESL + N*(N-1)/2
      IF (N.EQ.1) LPHESD = LPHESD + 1
      LP = LPHESD + N
      LP2 = LP + N
      LGSSQ = LP2 + N
      LUTF = LGSSQ + N
      LHK = LUTF + M
      LEND = LHK + LH
      NWHY = 1
      IF (LW.LT.LEND-1) GO TO 140
C
C     ALL THE WORKSPACE ADDRESSES ARE NOW ALLOCATED.
C
      CALL E04GDZ(M,N,E04GBX,LSFJC,MAXCAL,ETA,XTOL,STEPMX,X,SSQNEW,FVEC,
     *            FJAC,LJ,LVT,NITER,NFTOTL,NWHY,RTEPS,RTOL,LEND,TOL,
     *            IGRADE,PEPS,W(LGSSQ),SSQOLD,GAUSS,TAU,ALPHA,PNORM,
     *            EPSMCH,NOMOVE,IW,LIW,W,LW)
      IF (NWHY.NE.0) GO TO 140
      TOLEPS = RTOL + EPSMCH
      TLEPSQ = TOLEPS*TOLEPS
      EPSEPS = EPSMCH*EPSMCH
C
C     INITIALIZE THE APPROXIMATION TO THE SECOND DERIVATIVE TERM OF
C     THE HESSIAN MATRIX TO ZERO.
C
      CALL F06FBF(LH,0.0D0,W(LHK),1)
C
C     ..................START OF THE ITERATION LOOP.................
C
   20 NWHY = 2
      IF (NFTOTL.GT.MAXCAL) GO TO 100
      NWHY = 0
C
C     OVERALL CONVERGENCE CRITERION.
C
      GTG = DDOT(N,W(LGSSQ),1,W(LGSSQ),1)
      XNORM = DNRM2(N,X,1)
      U = 1.0D+0 + SSQNEW
      IF (NITER.GT.0 .AND. ALPHA*PNORM.LT.TOLEPS*(1.0D+0+XNORM)
     *    .AND. ABS(SSQOLD-SSQNEW).LT.TLEPSQ*U .AND. GTG.LT.PEPS*U*
     *    U .OR. SSQNEW.LT.EPSEPS .OR. GTG.LT.SQRT(SSQNEW)*EPSMCH)
     *    GO TO 100
C
C     COMPUTE THE SINGULAR VALUE DECOMPOSITION OF THE JACOBIAN
C     MATRIX.  DETERMINE THE GRADE OF THE JACOBIAN AND WHETHER THE
C     GAUSS-NEWTON DIRECTION SHOULD BE CORRECTED.
C
   40 IF (NOMOVE) CALL E04GDV(N,VT,LVT)
      CALL E04GDY(M,N,MAXRNK,IGRADE,NWHY,TAU,EPSMCH,FVEC,FJAC,LJ,S,VT,
     *            LVT,W(LUTF),GAUSS,NOMOVE,W,LW)
      IF (NWHY.NE.0) GO TO 140
C
C     COMPUTE THE GRADED GAUSS-NEWTON DIRECTION
C
      IF (IGRADE.EQ.0) CALL F06FBF(N,0.0D0,W(LP),1)
      IF (IGRADE.NE.0) CALL F04JAZ(M,N,IGRADE,S,N,W(LUTF),VT,LVT,W(LP),
     *                             SIGMA,W)
C
C     TRANSPOSE VT
C
      CALL E04GDV(N,VT,LVT)
C
C     IF NECESSARY, FIND A CORRECTION TO THE GAUSS-NEWTON DIRECTION
C     OF SEARCH.
C
      IF (GAUSS) GO TO 60
      NS = N - IGRADE
      IF (NS.EQ.0) GO TO 60
      LPH = NS*(NS-1)/2
      IF (NS.EQ.1) LPH = 1
C
C     FORM A PROJECTION OF THE QUASI-NEWTON APPROXIMATION TO THE SUM
C     OF SQUARES.
C
      CALL E04HEY(N,LH,LPH,NS,IGRADE,VT,LVT,W(LP),W(LPHESL),W(LPHESD),
     *            W(LRHS),W(LHK),W)
C
C     ADD THE SQUARE OF THE SINGULAR VALUE TO THE APPROPRIATE
C     DIAGONAL ELEMENT OF VTHV. MODIFY THE RIGHT-HAND-SIDE VECTOR
C     AND FORM THE MODIFIED LDLT FACTORIZATION.
C
      CALL E04GDS(M,N,NS,LPH,IGRADE,EPSMCH,W(LUTF),S,W(LPHESL),W(LPHESD)
     *            ,W(LRHS),NPHI)
C
C     USE THE CHOLESKY FACTORIZATION TO FIND THE CORRECTION TO THE
C     GAUSS-NEWTON DIRECTION OF SEARCH. FIRST FORM A PROJECTION
C     OF THE REQUIRED DIRECTION.
C
      CALL F04AQZ(NS,LPH,W(LPHESL),W(LPHESD),W(LRHS),W(LP2))
C
C     ADD V2*P2 TO THE VECTOR P.
C
      CALL E04GDR(N,NS,VT(1,IGRADE+1),LVT,W(LP2),W(LP))
C
C     THE GRADED GAUSS-NEWTON DIRECTION HAS BEEN CORRECTED.
C     STORE THE OLD GRADIENT OF THE SUM OF SQUARES SINCE IT IS
C     REQUIRED BY E04GBZ.
C
   60 CALL DCOPY(N,W(LGSSQ),1,W(LPHESD),1)
C
C     COMPUTE THE STEPLENGTH ALPHA. IF NECESSARY, PRINT DETAILS OF
C     THE CURRENT ITERATION.
C
      CALL E04GDQ(M,N,LSQLIN,E04GBX,LSFJC,RTEPS,ETA,IGRADE,IPRINT,LSMON,
     *            STEPMX,EPSMCH,TAU,XNORM,NOMOVE,X,FVEC,FJAC,LJ,W(LGSSQ)
     *            ,W(LP),S,ALPHA,PNORM,SSQNEW,SSQOLD,NITER,NFTOTL,NWHY,
     *            IFLAG,IW,LIW,W,LW)
      IF (NWHY.LT.0) GO TO 140
      NWHY = 0
C
C     IF IFLAG = 0 (DERIVATIVE LINEAR SEARCH) OR IFLAG = 2
C     (NON-DERIVATIVE LINEAR SEARCH), A LOWER POINT HAS
C     BEEN FOUND.
C
      IF (IFLAG.NE.0 .AND. IFLAG.NE.2) GO TO 80
C
C     UPDATE THE QUASI-NEWTON APPROXIMATION TO THE SECOND TERM OF
C     THE HESSIAN MATRIX.
C
      CALL E04GBZ(M,N,LH,W(LHK),ALPHA,W(LP),W(LGSSQ),W(LPHESD),FJAC,LJ,
     *            W,LW)
      GO TO 20
C
C     A LOWER POINT COULD NOT BE FOUND.
C
   80 NOMOVE = .TRUE.
      IF (IGRADE.GT.0) GO TO 40
      NWHY = 3
C
C     .....................END OF ITERATION LOOP....................
C
  100 FSUMSQ = SSQNEW
C
C     CALCULATE SINGULAR VALUES AT SOLUTION IF IT IS A NEW POINT,
C     AND TRANSPOSE ELEMENTS OF VT TO GIVE MATRIX V ON EXIT.
C     PRINT OUT DETAILS OF SOLUTION.
C
      IF (NOMOVE) GO TO 120
      IFLAG = NWHY
      NWHY = 4
      NWY = 1
      CALL E04GDX(NWY,M,N,FJAC,LJ,VT,LVT,FVEC,S,.TRUE.,W(LUTF),W,LW)
      IF (NWY.NE.0) GO TO 140
      NWHY = IFLAG
      CALL E04GDV(N,VT,LVT)
  120 IF (IPRINT.GE.0 .AND. NWHY.NE.1) CALL LSMON(M,N,X,FVEC,FJAC,LJ,S,
     *    IGRADE,NITER,NFTOTL,IW,LIW,W,LW)
  140 IF (NWHY.LT.0) THEN
         P01REC(1) = ' ** Negative value of IFLAG set in LSQFUN by user'
         NREC = 1
      ELSE IF (NWHY.EQ.1) THEN
         LW1 = 7*N + M*N + 2*M + N*N
         LW2 = 9 + 3*M
         IF (N.LT.1) THEN
            WRITE (P01REC,FMT=99999) N
            NREC = 1
         ELSE IF (M.LT.N) THEN
            WRITE (P01REC,FMT=99998) M, N
            NREC = 1
         ELSE IF (MAXCAL.LT.1) THEN
            WRITE (P01REC,FMT=99997) MAXCAL
            NREC = 1
         ELSE IF (ETA.LT.0.0D0 .OR. ETA.GE.1.0D0) THEN
            WRITE (P01REC,FMT=99996) ETA
            NREC = 1
         ELSE IF (XTOL.LT.0.0D0) THEN
            WRITE (P01REC,FMT=99995) XTOL
            NREC = 1
         ELSE IF (STEPMX.LT.XTOL) THEN
            WRITE (P01REC,FMT=99994) STEPMX, XTOL
            NREC = 2
         ELSE IF (LJ.LT.M) THEN
            WRITE (P01REC,FMT=99993) LJ, M
            NREC = 1
         ELSE IF (LVT.LT.N) THEN
            WRITE (P01REC,FMT=99992) LVT, N
            NREC = 1
         ELSE IF (LIW.LT.1) THEN
            WRITE (P01REC,FMT=99991) LIW
            NREC = 1
         ELSE IF (LW.LT.LW1 .AND. N.GT.1) THEN
            WRITE (P01REC,FMT=99990) LW, LW1
            NREC = 2
         ELSE IF (LW.LT.LW2 .AND. N.EQ.1) THEN
            WRITE (P01REC,FMT=99989) LW, LW2
            NREC = 2
         END IF
      ELSE IF (NWHY.EQ.2) THEN
         WRITE (P01REC,FMT=99988) MAXCAL
         NREC = 1
      ELSE IF (NWHY.EQ.3) THEN
         P01REC(1) =
     *   ' ** The conditions for a minimum have not all been satisfied,'
         P01REC(2) = ' ** but a lower point could not be found'
         NREC = 2
      ELSE IF (NWHY.EQ.4) THEN
         P01REC(1) = ' ** Failure in computing SVD of Jacobian matrix'
         NREC = 1
      END IF
      IFAIL = P01ABF(IFAIL,NWHY,SRNAME,NREC,P01REC)
      RETURN
C
C     END OF E04GBF   (LSFDQ)
C
99999 FORMAT (' ** On entry, N must be at least 1: N =',I16)
99998 FORMAT (' ** On entry, M must be at least N: M =',I16,', N =',I16)
99997 FORMAT (' ** On entry, MAXCAL must be at least 1: MAXCAL =',I16)
99996 FORMAT (' ** On entry, ETA must satisfy 0.le.ETA.lt.1: ETA =',1P,
     *       D13.5)
99995 FORMAT (' ** On entry, XTOL must be at least 0.0: XTOL =',1P,
     *       D13.5)
99994 FORMAT (' ** On entry, STEPMX must be at least XTOL:',/' ** STEP',
     *       'MX =',1P,D13.5,', XTOL =',1P,D13.5)
99993 FORMAT (' ** On entry, LJ must be at least M: LJ =',I16,', M =',
     *       I16)
99992 FORMAT (' ** On entry, LV must be at least N: LV =',I16,', N =',
     *       I16)
99991 FORMAT (' ** On entry, LIW must be at least 1: LIW =',I16)
99990 FORMAT (' ** On entry, LW must be at least 7*N + M*N + 2*M + N*N',
     *       ' if N.gt.1:',/' ** LW =',I16,', LW must be at least',I16)
99989 FORMAT (' ** On entry, LW must be at least 9 + 3*M if N.eq.1:',
     *       /' ** LW =',I16,', LW must be at least',I16)
99988 FORMAT (' ** There have been MAXCAL calls to LSQFUN: MAXCAL =',
     *       I16)
      END
