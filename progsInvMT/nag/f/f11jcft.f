      SUBROUTINE F11JCF(METHOD,N,NNZ,A,LA,IROW,ICOL,IPIV,ISTR,B,TOL,
     *                  MAXITN,X,RNORM,ITN,WORK,LWORK,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Purpose
C     =======
C     This black box routine uses F11GAF, F11GBF, F11JBF and F11XEF to
C     solve a sparse symmetric linear system
C
C               A x = b
C
C     F11JCF uses the incomplete Cholesky factorization determined by
C     F11JAF as the preconditioning matrix. A call to F11JCF must always
C     be preceeded by a call to F11JAF. Alternative preconditioners for
C     the same storage scheme are available by calling F11JEF.
C
C     The matrix A, and the preconditioning matrix M are represented in
C     symmetric coordinate storage (SCS) format in the arrays A, IROW
C     and ICOL, as returned from F11JAF. The array A holds the non-zero
C     entries in the lower triangular parts of these matrices, while
C     IROW and ICOL hold the corresponding row and column indices.
C
C     Arguments
C     =========
C
C     METHOD (input) CHARACTER*(*)
C            On entry, the iterative method to be used:
C               'CG'      ==> Conjugate gradient method
C               'SYMMLQ'  ==> Lanczos method (SYMMLQ).
C            METHOD = 'CG' or 'SYMMLQ'.
C
C     N      (input) INTEGER
C            On entry, the order of the matrix A.
C            N >= 1.
C
C     NNZ    (input) INTEGER
C            On entry, the number of non-zero elements in the lower
C            triangular part of A.
C            1 <= NNZ <= N*(N+1)/2.
C
C     A      (input) DOUBLE PRECISION array, dimension (LA)
C            On entry, A must be unchanged from a previous call to
C            F11JAF.
C
C     LA     (input) INTEGER
C            On entry, the dimension of the arrays A, IROW and ICOL
C            as declared in the (sub)program from which F11JCF is
C            called.
C            LA >= NNZ + NNZC, where NNZC is the value returned
C            from F11JAF.
C
C     IROW   (input) INTEGER array, dimension (LA)
C     ICOL   (input) INTEGER array, dimension (LA)
C     IPIV   (input) INTEGER array, dimension (N)
C     ISTR   (input) INTEGER array, dimension (N+1)
C            On entry, IROW, ICOL, IPIV and ISTR must be unchanged from
C            a previous call to F11JAF.
C
C     B      (input) DOUBLE PRECISION array, dimension (N)
C            On entry, the right-hand side vector b.
C
C     TOL    (input) DOUBLE PRECISION
C            On entry, the tolerance required. The iteration is judged
C            to have converged at step k if:
C
C               || r || <= TOL*(|| b || + || A ||*|| x ||).
C                   k                                 k
C            If TOL = 0.0 the default value of SQRT(machine precision)
C            is used.
C            TOL < 1.0.
C
C     MAXITN (input) INTEGER
C            On entry, the maximum number of iterations allowed.
C            MAXITN >= 1
C
C     X      (input/output) DOUBLE PRECISION array, dimension (N)
C            On entry, an initial approximation of the vector x.
C            On exit, an improved approximation to the vector x.
C
C     RNORM  (output) DOUBLE PRECISION
C            On exit, the final residual norm.
C
C     ITN    (output) INTEGER
C            On exit, the actual number of iterations taken.
C
C     WORK   (workspace) DOUBLE PRECISION array, dimension (LWORK)
C
C     LWORK  (input) INTEGER
C            On entry, the dimension of the array WORK as declared in
C            the (sub)program from which F11JCF is called.
C               METHOD = 'CG'      ==> LWORK >= 6*N
C               METHOD = 'SYMMLQ'  ==> LWORK >= 7*N.
C
C     IFAIL  (input/output) INTEGER
C            On entry, IFAIL must be -1, 0, or 1.
C            On exit, the following values may occur:
C               IFAIL = 0 => no error detected.
C               IFAIL = 1 => METHOD invalid, or
C                            N < 1, or
C                            NNZ < 1, or
C                            NNZ > N*(N+1)/2, or
C                            LA too small, or
C                            TOL >= 1.0, or
C                            MAXITN < 1, or
C                            LWORK too small.
C               IFAIL = 2 => SCS representation of A invalid.
C               IFAIL = 3 => SCS representation of M invalid.
C               IFAIL = 4 => Required accuracy not obtained. However,
C                            a reasonable accuracy has been obtained
C                            and further iterations could not improve
C                            the result.
C               IFAIL = 5 => Required accuracy not obtained in MAXITN
C                            iterations.
C               IFAIL = 6 => The preconditioner appears not to be
C                            positive definite.
C               IFAIL = 7 => The matrix A appears not to be positive
C                            definite (CG only).
C               IFAIL = 8 => A serious error has occurred in an internal
C                            call to F11GAF, F11GBF, or F11GCF.
C
C     ==================================================================
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F11JCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RNORM, TOL
      INTEGER           IFAIL, ITN, LA, LWORK, MAXITN, N, NNZ
      CHARACTER*(*)     METHOD
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LA), B(N), WORK(LWORK), X(N)
      INTEGER           ICOL(LA), IPIV(N), IROW(LA), ISTR(N+1)
C     .. Local Scalars ..
      DOUBLE PRECISION  ABA, ANORM, SIGERR, SIGMAX, SIGTOL, STPRHS
      INTEGER           I, IBAD, ICI, IERGAF, IERGBF, IERGCF, IERJAZ,
     *                  IERJBF, IERR, IERXEF, IREVCM, IRI, ISTW, ITERM,
     *                  ITS, LW, LWREQ, MAXITS, MONIT, NREC
      LOGICAL           MOK, SYM
      CHARACTER         CHECKA, CHECKM, NORM, PR, SIGCMP, WEIGHT
C     .. Local Arrays ..
      CHARACTER*80      REC(4)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, F06FBF, F11GAF, F11GBF, F11GCF, F11JAZ,
     *                  F11JBF, F11XEF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Executable Statements ..
C
      IERR = 0
C
C     Check the validity of the SCS representation of A.
C
      SYM = .TRUE.
      CALL F11JAZ(N,NNZ,IROW,ICOL,SYM,IBAD,IERJAZ)
C
      IF (IERJAZ.EQ.1) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99999) N
         GO TO 80
      END IF
C
      IF (IERJAZ.EQ.2) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99998) NNZ
         GO TO 80
      END IF
C
      IF (IERJAZ.EQ.3) THEN
         IERR = 1
         NREC = 2
         WRITE (REC,FMT=99997) NNZ, N
         GO TO 80
      END IF
C
      IF (IERJAZ.EQ.4) THEN
         IERR = 2
         NREC = 4
         WRITE (REC,FMT=99996) IBAD, IROW(IBAD), N
         GO TO 80
      END IF
C
      IF (IERJAZ.EQ.5) THEN
         IERR = 2
         NREC = 4
         WRITE (REC,FMT=99995) IBAD, ICOL(IBAD), IROW(IBAD)
         GO TO 80
      END IF
C
      IF (IERJAZ.EQ.6) THEN
         IERR = 2
         NREC = 4
         WRITE (REC,FMT=99994) IBAD
         GO TO 80
      END IF
C
      IF (IERJAZ.EQ.7) THEN
         IERR = 2
         NREC = 4
         WRITE (REC,FMT=99993) IBAD
         GO TO 80
      END IF
C
C     Check LA.
C
      IF (LA.LT.2*NNZ) THEN
         IERR = 1
         NREC = 2
         WRITE (REC,FMT=99992) LA, NNZ
         GO TO 80
      END IF
C
C     Find the infinity norm of A.
C
      CALL F06FBF(N,0.D0,WORK,1)
      DO 20 I = 1, NNZ
         ABA = ABS(A(I))
         IRI = IROW(I)
         ICI = ICOL(I)
         WORK(IRI) = WORK(IRI) + ABA
         IF (IRI.NE.ICI) WORK(ICI) = WORK(ICI) + ABA
   20 CONTINUE
      ANORM = 0.D0
      DO 40 I = 1, N
         ANORM = MAX(ANORM,WORK(I))
   40 CONTINUE
C
C     Copy right-hand side vector b into WORK.
C
      CALL DCOPY(N,B,1,WORK,1)
      ISTW = N + 1
      LW = LWORK - N
C
C     Call F11GAF to initialize solver.
C
      PR = 'P'
      SIGCMP = 'N'
      NORM = 'I'
      WEIGHT = 'N'
      ITERM = 1
      SIGMAX = 0.D0
      SIGTOL = 0.0D0
      MAXITS = 1
      MONIT = 0
      MOK = .TRUE.
      IERGAF = 1
C
      CALL F11GAF(METHOD,PR,SIGCMP,NORM,WEIGHT,ITERM,N,TOL,MAXITN,ANORM,
     *            SIGMAX,SIGTOL,MAXITS,MONIT,LWREQ,IERGAF)
C
      IF (IERGAF.EQ.-1) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99991) METHOD
         GO TO 80
      ELSE IF (IERGAF.EQ.-8) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99990) TOL
         GO TO 80
      ELSE IF (IERGAF.EQ.-9) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99989) MAXITN
         GO TO 80
      ELSE IF (IERGAF.NE.0) THEN
         IERR = 8
         NREC = 2
         WRITE (REC,FMT=99982) 'F11GAF', IERGAF
         GO TO 80
      END IF
C
C     Call F11GBF to solve the linear system.
C
      IREVCM = 0
      CHECKM = 'C'
C
   60 CONTINUE
C
      IERGBF = 1
      CALL F11GBF(IREVCM,X,WORK(1),WORK(ISTW),LW,IERGBF)
C
      IF (IREVCM.EQ.1) THEN
C
C         Compute matrix vector product.
C
         IERXEF = 1
         CHECKA = 'N'
         CALL F11XEF(N,NNZ,A,IROW,ICOL,CHECKA,X,WORK(1),IERXEF)
         GO TO 60
C
      ELSE IF (IREVCM.EQ.2) THEN
C
C         Solve linear system involving preconditioning matrix.
C
         IERJBF = 1
         CALL F11JBF(N,A,LA,IROW,ICOL,IPIV,ISTR,CHECKM,X,WORK(1),IERJBF)
         CHECKM = 'N'
         IF (IERJBF.NE.0) THEN
            MOK = .FALSE.
            IREVCM = 5
         END IF
C
         GO TO 60
C
      ELSE IF (IREVCM.EQ.4) THEN
C
C         Check error exits from F11JBF.
C
         IF ( .NOT. MOK) THEN
            IERR = 3
            NREC = 3
            WRITE (REC,FMT=99988)
            GO TO 80
         END IF
C
C         Check error exits from F11GBF.
C
         IF (IERGBF.EQ.-5) THEN
            IERR = 1
            NREC = 2
            WRITE (REC,FMT=99987) LWORK, ISTW + LWREQ - 1
            GO TO 80
         ELSE IF (IERGBF.EQ.2) THEN
            IERR = 4
            NREC = 2
            WRITE (REC,FMT=99986)
            GO TO 80
         ELSE IF (IERGBF.EQ.5) THEN
            IERR = 5
            NREC = 1
            WRITE (REC,FMT=99985) MAXITN
            GO TO 80
         ELSE IF (IERGBF.EQ.6) THEN
            IERR = 6
            NREC = 2
            WRITE (REC,FMT=99984)
            GO TO 80
         ELSE IF (IERGBF.EQ.7) THEN
            IERR = 7
            NREC = 2
            WRITE (REC,FMT=99983)
            GO TO 80
         ELSE IF (IERGBF.NE.0) THEN
            IERR = 8
            NREC = 2
            WRITE (REC,FMT=99982) 'F11GBF', IERGBF
            GO TO 80
         END IF
C
         IERGCF = 0
         CALL F11GCF(ITN,RNORM,STPRHS,ANORM,SIGMAX,ITS,SIGERR,IERGCF)
C
         IF (IERGCF.NE.0) THEN
            IERR = 8
            NREC = 2
            WRITE (REC,FMT=99982) 'F11GCF', IERGCF
            GO TO 80
         END IF
C
      ELSE
C
         IERR = 8
         NREC = 2
         WRITE (REC,FMT=99981) IREVCM
         GO TO 80
C
      END IF
C
   80 CONTINUE
C
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, N .lt. 1: N =',I16,'.')
99998 FORMAT (1X,'** On entry, NNZ .lt. 1: NNZ =',I16,'.')
99997 FORMAT (1X,'** On entry, NNZ .gt. N*(N+1)/2:',/4X,'NNZ =',I16,
     *       ' N =',I16,'.')
99996 FORMAT (1X,'** On entry, IROW(I) .lt. 1 or IROW(I) .gt. N:',/4X,
     *       'I = ',I16,', IROW(I) = ',I16,', N = ',I16,'.',/4X,'Check',
     *       ' that A, IROW, ICOL, IPIV and ISTR have not been',/4X,
     *       'corrupted between calls to F11JAF and F11JCF.')
99995 FORMAT (1X,'** On entry, ICOL(I) .lt. 1 or ICOL(I) .gt. IROW(I):',
     *       /4X,'I =',I16,', ICOL(I) =',I16,', IROW(I) =',I16,'.',/4X,
     *       'Check that A, IROW, ICOL, IPIV and ISTR have not been',
     *       /4X,'corrupted between calls to F11JAF and F11JCF.')
99994 FORMAT (1X,'** On entry, A(I) is out of order:',/4X,'I =',I16,'.',
     *       /4X,'Check that A, IROW, ICOL, IPIV and ISTR have not been'
     *       ,/4X,'corrupted between calls to F11JAF and F11JCF.')
99993 FORMAT (1X,'** On entry, the location (IROW(I), ICOL(I)) is a ',
     *       'duplicate:',/4X,'I =',I16,'.',/4X,'Check that A, IROW, I',
     *       'COL, IPIV and ISTR have not been',/4X,'corrupted between',
     *       ' calls to F11JAF and F11JCF.')
99992 FORMAT (1X,'** On entry, LA .lt. 2*NNZ:',/4X,'LA =',I16,' NNZ =',
     *       I16,'.')
99991 FORMAT (1X,'** On entry, METHOD .ne. ''CG'' or ''SYMMLQ'': ','ME',
     *       'THOD = ''',A,'''.')
99990 FORMAT (1X,'** On entry, TOL .ge. 1.0: TOL =',E12.4,'.')
99989 FORMAT (1X,'** On entry, MAXITN .lt. 1: MAXITN =',I16,'.')
99988 FORMAT (1X,'** The SCS representation of the preconditioner is ',
     *       'invalid.',/4X,'Check that A, IROW, ICOL, IPIV and ISTR h',
     *       'ave not been',/4X,'corrupted between calls to F11JAF and',
     *       ' F11JCF.')
99987 FORMAT (1X,'** On entry, LWORK is too small: LWORK =',I16,'.',/4X,
     *       'Minimum required value of LWORK =',I16,'.')
99986 FORMAT (1X,'** The required accuracy could not be obtained.',/4X,
     *       'However a reasonable accuracy has been achieved.')
99985 FORMAT (1X,'** The solution has not converged after',I7,' iterat',
     *       'ions.')
99984 FORMAT (1X,'** The preconditioner appears not to be positive-',
     *       'definite.',/4X,'The computation cannot continue.')
99983 FORMAT (1X,'** The matrix of the coefficients A appears not to ',
     *       'be positive-definite.',/4X,'The computation cannot conti',
     *       'nue.')
99982 FORMAT (1X,'** A serious error has occurred in an internal call ',
     *       'to ',A6,': IFAIL =',I6,/4X,'Check all subroutine ','call',
     *       's and array sizes. Seek expert help.')
99981 FORMAT (1X,'** A serious error has occurred in an internal call ',
     *       'to F11GBF: IREVCM =',I6,/4X,'Check all subroutine ','cal',
     *       'ls and array sizes. Seek expert help.')
      END
