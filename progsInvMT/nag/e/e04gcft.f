      SUBROUTINE E04GCF(M,N,X,FSUMSQ,IW,LIW,W,LW,IFAIL)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C
C     **************************************************************
C
C     E04GCF IS AN EASY-TO-USE QUASI-NEWTON ALGORITHM FOR FINDING
C     AN UNCONSTRAINED MINIMUM OF A SUM OF SQUARES OF M NONLINEAR
C     FUNCTIONS IN N VARIABLES (M .GE. N).  FIRST DERIVATIVES ARE
C     REQUIRED.
C
C     THE ROUTINE IS ESSENTIALLY IDENTICAL TO THE SUBROUTINE LSFDQ2
C     IN THE NPL ALGORITHMS LIBRARY (REF. NO. E4/21/F) AND CALLS
C     E04GBF WITH SUITABLE DEFAULT SETTINGS FOR PARAMETERS.  IT
C     CALLS THE USER-SUPPLIED ROUTINE LSFUN2.
C
C     N.B. LSFUN2 IS A DESIGNATED NAME.
C     --------------------------------
C
C     GIVEN AN INITIAL APPROXIMATION TO THE MINIMUM, E04GCF COMPUTES
C     THE POSITION OF THE MINIMUM AND THE CORRESPONDING FUNCTION
C     VALUE.  THE REAL ARRAY W AND THE INTEGER ARRAY IW ARE USED
C     AS WORKSPACE.  W MUST BE DIMENSIONED AT LEAST
C     8*N + N*N + 2*M*N + 3*M + N*(N - 1)/2 + N*(N + 1)/2 OR
C     11 + 5*M IF N = 1 AND IW MUST BE HAVE AT LEAST ONE ELEMENT.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN
C     AND NICHOLAS I. M. GOULD.
C     D.N.A.C., NATIONAL PHYSICAL LABORATORY, ENGLAND.
C
C     **************************************************************
C
C     Modified to call BLAS.
C     Peter Mayes, NAG Central Office, October 1987.
C
C     Modified to output explanatory messages.
C     Peter Mayes, NAG Central Office, December 1987.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E04GCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  FSUMSQ
      INTEGER           IFAIL, LIW, LW, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  W(LW), X(N)
      INTEGER           IW(LIW)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPSMCH, ETA, GTG, PEPS, STEPMX, U, XTOL
      INTEGER           IPRINT, LFJAC, LFVEC, LS, LV, LW1, LW2, MAXFUN,
     *                  NFTOTL, NITER, NREC, NWHY
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, X02AJF
      INTEGER           P01ABF
      EXTERNAL          DDOT, X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          E04FDZ, E04GBF, E04GCZ, E04HEV, E04YAF, DGEMV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      EPSMCH = X02AJF()
      PEPS = EPSMCH**0.66666D+0
C
C     CHECK THAT THE WORKSPACE HAS BEEN ALLOCATED CORRECTLY.
C
      NWHY = 1
      IF (LW.LT.8*N+N*N+2*M*N+3*M+N*(N-1)/2+N*(N+1)
     *    /2 .AND. N.GT.1 .OR. LW.LT.11+5*M .AND. N.EQ.1 .OR. LIW.LT.
     *    1 .OR. M.LT.N) GO TO 20
      NWHY = 0
C
C     SET UP INPUT PARAMETERS FOR E04GBF.
C
C     SUPPRESS THE PRINT FREQUENCY.
C
      IPRINT = 0
C
C     SET THE MAXIMUM NUMBER OF FUNCTION EVALUATIONS ALLOWED
C     AS 50*N.
C
      MAXFUN = 50*N
C
C     SET THE LINEAR SEARCH PARAMETER.
C
      ETA = 9.0D-1
      IF (N.EQ.1) ETA = 0.0D+0
C
C     SPECIFY THE OVERALL CONVERGENCE CRITERION.
C
      XTOL = 1.0D+1*SQRT(EPSMCH)
C
C     SPECIFY THE BOUND ON THE STEPLENGTH.
C
      STEPMX = 1.0D+1
C
C     COMPUTE THE ADDRESSES FOR THE ARRAYS USED IN E04GBF.
C
      LS = 7*N + M*N + 2*M + N*(N-1)/2 + N*(N+1)/2 + 1
      IF (N.EQ.1) LS = LS + 1
      LV = LS + N
      LFVEC = LV + N*N
      LFJAC = LFVEC + M
C
C     CHECK THAT THE JACOBIAN MATRIX OF THE FUNCTION VECTOR HAS
C     BEEN PROGRAMMED CORRECTLY BY CALLING E04YAF WITH SOFT
C     FAILURE OPTION.
C
      NWHY = 1
      CALL E04YAF(M,N,E04GCZ,X,W(LFVEC),W(LFJAC),M,IW,LIW,W,LW,NWHY)
      IF (NWHY.EQ.2) NWHY = 9
      IF (NWHY.NE.0) GO TO 20
C
C     CALL E04GBF (WITH SOFT FAILURE OPTION) TO MINIMIZE
C     THE FUNCTION.
C
      NWHY = 1
      CALL E04GBF(M,N,E04HEV,E04GCZ,E04FDZ,IPRINT,MAXFUN,ETA,XTOL,
     *            STEPMX,X,FSUMSQ,W(LFVEC),W(LFJAC),M,W(LS),W(LV),N,
     *            NITER,NFTOTL,IW,LIW,W,LW,NWHY)
      IF (NWHY.LE.2) GO TO 20
      U = 1.0D+0 + ABS(FSUMSQ)
      U = PEPS*U*U
      CALL DGEMV('Transpose',M,N,2.0D0,W(LFJAC),M,W(LFVEC),1,0.0D0,W,1)
      GTG = DDOT(N,W,1,W,1)
C
C     ATTEMPT TO DETERMINE WHETHER FAILURE WAS DUE TO XTOL BEING SET
C     TOO SMALL. THE FOLLOWING RULE IS USED -
C
C     NWHY = 3  -  THE MINIMIZATION HAS FAILED.
C     NWHY = 5  -  THE MINIMIZATION HAS PROBABLY WORKED.
C     NWHY = 6  -  THE MINIMIZATION HAS POSSIBLY WORKED.
C     NWHY = 7  -  THE MINIMIZATION IS UNLIKELY TO HAVE WORKED.
C     NWHY = 8  -  THE MINIMIZATION IS VERY UNLIKELY TO HAVE WORKED.
C
      IF (GTG.LE.1.0D+3*U) NWHY = 8
      IF (GTG.LE.1.0D+2*U) NWHY = 7
      IF (GTG.LE.1.0D+1*U) NWHY = 6
      IF (GTG.LE.U) NWHY = 5
   20 IF (NWHY.EQ.1) THEN
         LW1 = 8*N + 2*N*N + 2*M*N + 3*M
         LW2 = 11 + 5*M
         IF (N.LT.1) THEN
            WRITE (P01REC,FMT=99999) N
            NREC = 1
         ELSE IF (M.LT.N) THEN
            WRITE (P01REC,FMT=99998) M, N
            NREC = 1
         ELSE IF (LIW.LT.1) THEN
            WRITE (P01REC,FMT=99997) LIW
            NREC = 1
         ELSE IF (LW.LT.LW1 .AND. N.GT.1) THEN
            WRITE (P01REC,FMT=99996) LW, LW1
            NREC = 2
         ELSE IF (LW.LT.LW2 .AND. N.EQ.1) THEN
            WRITE (P01REC,FMT=99995) LW, LW2
            NREC = 2
         END IF
      ELSE IF (NWHY.EQ.2) THEN
         P01REC(1) = ' ** There have been 50*N calls to LSFUN2'
         NREC = 1
      ELSE IF (NWHY.EQ.3) THEN
         P01REC(1) =
     *   ' ** The conditions for a minimum have not all been satisfied,'
         P01REC(2) = ' ** but a lower point could not be found'
         NREC = 2
      ELSE IF (NWHY.EQ.4) THEN
         P01REC(1) = ' ** Failure in computing SVD of Jacobian matrix'
         NREC = 1
      ELSE IF (NWHY.EQ.5) THEN
         P01REC(1) =
     *         ' ** It is probable that a local minimum has been found,'
         P01REC(2) = ' ** but it cannot be guaranteed'
         NREC = 2
      ELSE IF (NWHY.EQ.6) THEN
         P01REC(1) =
     *         ' ** It is possible that a local minimum has been found,'
         P01REC(2) = ' ** but it cannot be guaranteed'
         NREC = 2
      ELSE IF (NWHY.EQ.7) THEN
         P01REC(1) =
     *          ' ** It is unlikely that a local minimum has been found'
         NREC = 1
      ELSE IF (NWHY.EQ.8) THEN
         P01REC(1) =
     *     ' ** It is very unlikely that a local minimum has been found'
         NREC = 1
      ELSE IF (NWHY.EQ.9) THEN
         P01REC(1) =
     *        ' ** It is very likely that the user has made an error in'
         P01REC(2) = ' ** forming the derivatives in LSFUN2'
         NREC = 2
      END IF
      IFAIL = P01ABF(IFAIL,NWHY,SRNAME,NREC,P01REC)
      RETURN
C
C     END OF E04GCF   (LSFDQ2)
C
99999 FORMAT (' ** On entry, N must be at least 1: N =',I16)
99998 FORMAT (' ** On entry, M must be at least N: M =',I16,', N =',I16)
99997 FORMAT (' ** On entry, LIW must be at least 1: LIW =',I16)
99996 FORMAT (' ** On entry, LW must be at least 8*N + 2*N*N + 2*M*N +',
     *       ' 3*M if N.gt.1:',/' ** LW =',I16,', LW must be at least',
     *       I16)
99995 FORMAT (' ** On entry, LW must be at least 11 + 5*M if N.eq.1:',
     *       /' ** LW =',I16,', LW must be at least',I16)
      END
