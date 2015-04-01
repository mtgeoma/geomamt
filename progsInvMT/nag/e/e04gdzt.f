      SUBROUTINE E04GDZ(M,N,LSFJSA,LSF,MAXFUN,ETA,XTOL,STEPMX,X,SSQNEW,
     *                  FVEC,FJAC,LJ,LVT,NITER,NFTOTL,NWHY,RTEPS,RTOL,
     *                  LEND,TOL,IGRADE,PEPS,G,SSQOLD,GAUSS,TAU,ALPHA,
     *                  PNORM,EPSMCH,NOMOVE,IW,LIW,W,LW)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C
C     **************************************************************
C
C     CHECK INPUT PARAMETERS, SET INITIAL VALUES OF LOCAL VARIABLES
C     AND, IF NECESSARY, CALL PRINT ROUTINE.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN
C     AND BRIAN T. HINDE
C     D.N.A.C., NATIONAL PHYSICAL LABORATORY, ENGLAND.
C
C     **************************************************************
C
C     Modified to call BLAS.
C     Peter Mayes, NAG Central Office, October 1987.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA, EPSMCH, ETA, PEPS, PNORM, RTEPS, RTOL,
     *                  SSQNEW, SSQOLD, STEPMX, TAU, TOL, XTOL
      INTEGER           IGRADE, LEND, LIW, LJ, LVT, LW, M, MAXFUN, N,
     *                  NFTOTL, NITER, NWHY
      LOGICAL           GAUSS, NOMOVE
C     .. Array Arguments ..
      DOUBLE PRECISION  FJAC(LJ,N), FVEC(M), G(N), W(LW), X(N)
      INTEGER           IW(LIW)
C     .. Subroutine Arguments ..
      EXTERNAL          LSF, LSFJSA
C     .. Local Scalars ..
      DOUBLE PRECISION  RMIN
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, X02AJF, X02AMF
      EXTERNAL          DDOT, X02AJF, X02AMF
C     .. External Subroutines ..
      EXTERNAL          DGEMV
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      EPSMCH = X02AJF()
      RMIN = X02AMF()
      RTEPS = SQRT(EPSMCH)
      RTOL = XTOL
      NITER = 0
      NFTOTL = 0
      NWHY = 1
C
C     CHECK FOR AN ERROR IN THE INPUT PARAMETERS.
C
      IF (LW.LT.LEND-1 .OR. LIW.LT.1 .OR. N.LT.1 .OR. M.LT.N .OR.
     *    RTOL.LT.0.0D+0 .OR. ETA.GE.1.0D+0 .OR. ETA.LT.0.0D+0 .OR.
     *    STEPMX.LT.RTOL .OR. MAXFUN.LT.1 .OR. LJ.LT.M .OR. LVT.LT.N)
     *    RETURN
      TOL = 1.0D+1*EPSMCH
      IF (RTOL.LT.TOL) RTOL = TOL
      TOL = RMIN/EPSMCH
      IGRADE = N
      PEPS = EPSMCH**0.6666D+0
      NWHY = 2
C
C     EVALUATE THE FUNCTION VALUES AND JACOBIAN MATRIX AT THE
C     STARTING POINT.
C
      CALL LSFJSA(NWHY,M,N,LSF,X,FVEC,FJAC,LJ,IW,LIW,W,LW)
      NFTOTL = 1
      IF (NWHY.LT.0) RETURN
      NWHY = 0
C
C     COMPUTE THE FUNCTION VALUE.
C
      SSQNEW = DDOT(M,FVEC,1,FVEC,1)
C
C     COMPUTE AND STORE THE GRADIENT VECTOR IN G.
C
      CALL DGEMV('Transpose',M,N,2.0D0,FJAC,LJ,FVEC,1,0.0D0,G,1)
      SSQOLD = SSQNEW
C
C     GAUSS = .TRUE. INDICATES THAT A GAUSS-NEWTON STEP WAS TAKEN
C     LAST ITERATION.
C
      NOMOVE = .FALSE.
      GAUSS = .TRUE.
      TAU = 1.0D+1
      ALPHA = 0.0D+0
      PNORM = 0.0D+0
      RETURN
C
C     END OF E04GDZ   (SETUP2)
C
      END
