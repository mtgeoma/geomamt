      SUBROUTINE F11JEF(METHOD,PRECON,N,NNZ,A,IROW,ICOL,OMEGA,B,TOL,
     *                  MAXITN,X,RNORM,ITN,WORK,LWORK,IWORK,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Purpose
C     =======
C     This black box routine uses F11GAF, F11GBF and F11XEF to solve a
C     sparse symmetric linear system
C
C               A x = b
C
C     by (preconditioned) conjugate gradients or SYMMLQ. The routine
C     allows the following choices for the preconditioner:
C
C               no preconditioning
C               Jacobi preconditioning
C               SSOR preconditioning
C
C     The matrix A is represented in symmetric coordinate storage (SCS)
C     format in the arrays A, IROW and ICOL. A holds the non-zero
C     entries in the lower triangular part of the matrix, while IROW
C     and ICOL hold the corresponding row and column indices.
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
C     PRECON (input) CHARACTER*1
C            On entry, determines the type of preconditioning used:
C               'N'       ==> no preconditioning
C               'J'       ==> Jacobi preconditioning.
C               'S'       ==> symmetric SOR preconditioning.
C            PRECON = 'N', 'J', 'S'.
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
C     A      (input) DOUBLE PRECISION array, dimension (NNZ)
C            On entry, the non-zero elements of the lower triangular
C            part of the matrix A, ordered by increasing row index
C            and by increasing column index within each row. Multiple
C            elements with the same row and column indices are not
C            allowed.
C
C     IROW   (input) INTEGER array, dimension (NNZ)
C     ICOL   (input) INTEGER array, dimension (NNZ)
C            On entry, the row and column indices corresponding to the
C            non-zero elements given in the array A.
C            IROW and ICOL must satisfy the following constraints:
C            1 <= IROW(i) <= N, and 1 <= ICOL(i) <= IROW(i), for
C            i = 1,2,...,NNZ.
C            IROW(i-1) < IROW(i), or
C            IROW(i-1) = IROW(i) and ICOL(i-1) < ICOL(i), for
C            i = 2,3,...,NNZ.
C
C     OMEGA  (input) DOUBLE PRECISION
C            On entry, if PRECON = 'S' then OMEGA is the relaxation
C            parameter to be used. Otherwise it is not referenced.
C            0.0 <= OMEGA <= 2.0.
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
C            the (sub)program from which F11JEF is called.
C               METHOD = 'CG'      ==> LWORK >= 6*N + M
C               METHOD = 'SYMMLQ'  ==> LWORK >= 7*N + M
C            where M = N for PRECON = 'J', or 'S' and 0 otherwise.
C
C     IWORK  (workspace) INTEGER array, dimension (N+1)
C
C     IFAIL  (input/output) INTEGER
C            On entry, IFAIL must be -1, 0, or 1.
C            On exit, the following values may occur:
C               IFAIL = 0 => no error detected.
C               IFAIL = 1 => METHOD invalid, or
C                            PRECON invalid, or
C                            N < 1, or
C                            NNZ < 1, or
C                            NNZ > N*(N+1)/2, or
C                            OMEGA not in [0,2], or
C                            TOL >= 1.0, or
C                            MAXITN < 1, or
C                            LWORK too small.
C               IFAIL = 2 => the arrays IROW and ICOL fail to satisfy
C                            the following constraints:
C                            1 <= IROW(i) <= N,
C                            1 <= ICOL(i) <= IROW(i), i = 1,2,...,NNZ.
C                            IROW(i-1) < IROW(i), or
C                            IROW(i-1) = IROW(i) and ICOL(i-1) < ICOL(i)
C                            for i = 2,3,...,NNZ.
C               IFAIL = 3 => A has a zero diagonal element.
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
      PARAMETER         (SRNAME='F11JEF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  OMEGA, RNORM, TOL
      INTEGER           IFAIL, ITN, LWORK, MAXITN, N, NNZ
      CHARACTER         PRECON
      CHARACTER*(*)     METHOD
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NNZ), B(N), WORK(LWORK), X(N)
      INTEGER           ICOL(NNZ), IROW(NNZ), IWORK(N+1)
C     .. Local Scalars ..
      DOUBLE PRECISION  ABA, ANORM, SIGERR, SIGMAX, SIGTOL, STPRHS
      INTEGER           I, IBAD, ICI, IERGAF, IERGBF, IERGCF, IERJAZ,
     *                  IERR, IERXEF, IPREC, IREVCM, IRI, IRIM1, ISTW,
     *                  ITERM, ITS, J, LW, LWNEED, MAXITS, MONIT, NREC
      LOGICAL           SYM
      CHARACTER         CHECKA, NORM, PR, SIGCMP, WEIGHT
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, F06DBF, F06FBF, F06FCF, F11GAF, F11GBF,
     *                  F11GCF, F11JAZ, F11JDZ, F11XEF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Executable Statements ..
C
      IERR = 0
C
C     Check N, NNZ and the validity of the SCS representation of A.
C
      SYM = .TRUE.
      CALL F11JAZ(N,NNZ,IROW,ICOL,SYM,IBAD,IERJAZ)
C
      IF (IERJAZ.EQ.1) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99999) N
         GO TO 120
      END IF
C
      IF (IERJAZ.EQ.2) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99998) NNZ
         GO TO 120
      END IF
C
      IF (IERJAZ.EQ.3) THEN
         IERR = 1
         NREC = 2
         WRITE (REC,FMT=99997) NNZ, N
         GO TO 120
      END IF
C
      IF (IERJAZ.EQ.4) THEN
         IERR = 2
         NREC = 2
         WRITE (REC,FMT=99996) IBAD, IROW(IBAD), N
         GO TO 120
      END IF
C
      IF (IERJAZ.EQ.5) THEN
         IERR = 2
         NREC = 2
         WRITE (REC,FMT=99995) IBAD, ICOL(IBAD), IROW(IBAD)
         GO TO 120
      END IF
C
      IF (IERJAZ.EQ.6) THEN
         IERR = 2
         NREC = 2
         WRITE (REC,FMT=99994) IBAD
         GO TO 120
      END IF
C
      IF (IERJAZ.EQ.7) THEN
         IERR = 2
         NREC = 2
         WRITE (REC,FMT=99993) IBAD
         GO TO 120
      END IF
C
C     Check PRECON.
C
      IPREC = -1
      IF (PRECON.EQ.'N' .OR. PRECON.EQ.'n') IPREC = 0
      IF (PRECON.EQ.'J' .OR. PRECON.EQ.'j') IPREC = 1
      IF (PRECON.EQ.'S' .OR. PRECON.EQ.'s') IPREC = 2
C
      IF (IPREC.EQ.-1) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99992) PRECON
         GO TO 120
      END IF
C
C     Check OMEGA.
C
      IF (IPREC.EQ.2) THEN
         IF (OMEGA.LT.0.D0 .OR. OMEGA.GT.2.D0) THEN
            IERR = 1
            NREC = 1
            WRITE (REC,FMT=99991) OMEGA
            GO TO 120
         END IF
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
C     Initialize data required for JACOBI and SSOR preconditioning.
C
      ISTW = N + 1
      LW = LWORK - N
C
      IF (IPREC.EQ.1 .OR. IPREC.EQ.2) THEN
C
         ISTW = 2*N + 1
         LW = LWORK - 2*N
C
C         Store start addresses for each row in IWORK.
C
         CALL F06DBF(N+1,0,IWORK,1)
         IWORK(1) = 1
         DO 60 I = 2, NNZ
            IRI = IROW(I)
            IRIM1 = IROW(I-1)
            IF (IRI.GT.IRIM1) IWORK(IRI) = I
   60    CONTINUE
         IWORK(N+1) = NNZ + 1
C
C         Store reciprocal diagonal matrix elements in WORK.
C
         DO 80 I = 1, N
            J = IWORK(I+1) - 1
            IF (J.EQ.-1) THEN
               IERR = 3
               NREC = 1
               WRITE (REC,FMT=99990) I + 1
               GO TO 120
            END IF
            IF (IROW(J).NE.I .OR. ICOL(J).NE.I) THEN
               IERR = 3
               NREC = 1
               WRITE (REC,FMT=99990) I
               GO TO 120
            END IF
            IF (A(J).EQ.0.D0) THEN
               IERR = 3
               NREC = 1
               WRITE (REC,FMT=99989) I
               GO TO 120
            ELSE
               WORK(N+I) = 1.D0/A(J)
            END IF
   80    CONTINUE
C
      END IF
C
C     Copy right-hand side vector b into WORK.
C
      CALL DCOPY(N,B,1,WORK,1)
C
C     Call F11GAF to initialize solver.
C
      PR = 'P'
      IF (IPREC.EQ.0) PR = 'N'
      SIGCMP = 'N'
      NORM = 'I'
      WEIGHT = 'N'
      ITERM = 1
      SIGMAX = 0.D0
      SIGTOL = 0.D0
      MAXITS = 1
      MONIT = 0
      IERGAF = 1
C
      CALL F11GAF(METHOD,PR,SIGCMP,NORM,WEIGHT,ITERM,N,TOL,MAXITN,ANORM,
     *            SIGMAX,SIGTOL,MAXITS,MONIT,LWNEED,IERGAF)
C
      IF (IERGAF.EQ.-1) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99988) METHOD
         GO TO 120
      ELSE IF (IERGAF.EQ.-8) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99987) TOL
         GO TO 120
      ELSE IF (IERGAF.EQ.-9) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99986) MAXITN
         GO TO 120
      ELSE IF (IERGAF.NE.0) THEN
         IERR = 8
         NREC = 2
         WRITE (REC,FMT=99980) 'F11GAF', IERGAF
         GO TO 120
      END IF
C
C     Call F11GBF to solve the linear system.
C
      IREVCM = 0
C
  100 CONTINUE
C
      IERGBF = 1
      CALL F11GBF(IREVCM,X,WORK(1),WORK(ISTW),LW,IERGBF)
C
      IF (IREVCM.EQ.1) THEN
C
C         Compute matrix vector product.
C
         IERXEF = 0
         CHECKA = 'N'
         CALL F11XEF(N,NNZ,A,IROW,ICOL,CHECKA,X,WORK(1),IERXEF)
         GO TO 100
C
      ELSE IF (IREVCM.EQ.2) THEN
C
C         Solve linear system involving preconditioning matrix.
C
         IF (IPREC.EQ.1) THEN
C
C             Jacobi preconditioning.
C
            CALL DCOPY(N,X,1,WORK(1),1)
            CALL F06FCF(N,WORK(N+1),1,WORK(1),1)
C
         ELSE IF (IPREC.EQ.2) THEN
C
C             SSOR preconditioning.
C
            CALL F11JDZ(N,NNZ,A,ICOL,IWORK,WORK(N+1),OMEGA,X,WORK(1))
C
         END IF
C
         GO TO 100
C
      ELSE IF (IREVCM.EQ.4) THEN
C
C         Check error exits from F11GBF.
C
         IF (IERGBF.EQ.-5) THEN
            IERR = 1
            NREC = 2
            WRITE (REC,FMT=99985) LWORK, ISTW + LWNEED - 1
            GO TO 120
         ELSE IF (IERGBF.EQ.2) THEN
            IERR = 4
            NREC = 2
            WRITE (REC,FMT=99984)
            GO TO 120
         ELSE IF (IERGBF.EQ.5) THEN
            IERR = 5
            NREC = 1
            WRITE (REC,FMT=99983) MAXITN
            GO TO 120
         ELSE IF (IERGBF.EQ.6) THEN
            IERR = 6
            NREC = 2
            WRITE (REC,FMT=99982)
            GO TO 120
         ELSE IF (IERGBF.EQ.7) THEN
            IERR = 7
            NREC = 2
            WRITE (REC,FMT=99981)
            GO TO 120
         ELSE IF (IERGBF.NE.0) THEN
            IERR = 8
            NREC = 2
            WRITE (REC,FMT=99980) 'F11GBF', IERGBF
            GO TO 120
         END IF
C
         IERGCF = 0
         CALL F11GCF(ITN,RNORM,STPRHS,ANORM,SIGMAX,ITS,SIGERR,IERGCF)
C
         IF (IERGCF.NE.0) THEN
            IERR = 8
            NREC = 2
            WRITE (REC,FMT=99980) 'F11GCF', IERGCF
            GO TO 120
         END IF
C
      ELSE
C
         IERR = 8
         NREC = 2
         WRITE (REC,FMT=99979) IREVCM
         GO TO 120
C
      END IF
C
  120 CONTINUE
C
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, N .lt. 1: N =',I16,'.')
99998 FORMAT (1X,'** On entry, NNZ .lt. 1: NNZ =',I16,'.')
99997 FORMAT (1X,'** On entry, NNZ .gt. N*(N+1)/2:',/4X,'NNZ =',I16,
     *       ' N =',I16,'.')
99996 FORMAT (1X,'** On entry, IROW(I) .lt. 1 or IROW(I) .gt. N:',/4X,
     *       'I = ',I16,', IROW(I) = ',I16,', N = ',I16,'.')
99995 FORMAT (1X,'** On entry, ICOL(I) .lt. 1 or ICOL(I) .gt. IROW(I):',
     *       /4X,'I =',I16,', ICOL(I) =',I16,', IROW(I) =',I16,'.')
99994 FORMAT (1X,'** On entry, A(I) is out of order:',/4X,'I =',I16,'.')
99993 FORMAT (1X,'** On entry, the location (IROW(I), ICOL(I)) is a ',
     *       'duplicate:',/4X,'I =',I16,'.')
99992 FORMAT (1X,'** On entry, PRECON .ne. ''N'', ''J'' or ''S'':',' P',
     *       'RECON = ''',A,'''.')
99991 FORMAT (1X,'** On entry, OMEGA .lt. 0.0 or OMEGA .gt. 2.0:',' OM',
     *       'EGA = ',E12.4,'.')
99990 FORMAT (1X,'** The matrix A has no diagonal entry in row',I16,'.')
99989 FORMAT (1X,'** The matrix A has a zero diagonal entry in row',I16,
     *       '.')
99988 FORMAT (1X,'** On entry, METHOD .ne. ''CG'' or ''SYMMLQ'':',' ME',
     *       'THOD = ''',A,'''.')
99987 FORMAT (1X,'** On entry, TOL .ge. 1.0: TOL = ',E12.4,'.')
99986 FORMAT (1X,'** On entry, MAXITN .lt. 1: MAXITN =',I16,'.')
99985 FORMAT (1X,'** On entry, LWORK is too small: LWORK =',I16,'.',/4X,
     *       'Minimum required value of LWORK =',I16,'.')
99984 FORMAT (1X,'** The required accuracy could not be obtained.',/4X,
     *       'However a reasonable accuracy has been achieved.')
99983 FORMAT (1X,'** The solution has not converged after',I7,' iterat',
     *       'ions.')
99982 FORMAT (1X,'** The preconditioner appears not to be positive-',
     *       'definite.',/4X,'The computation cannot continue.')
99981 FORMAT (1X,'** The matrix of the coefficients A appears not ',
     *       'to be positive-definite.',/4X,'The computation cannot co',
     *       'ntinue.')
99980 FORMAT (1X,'** A serious error has occurred in an internal call ',
     *       'to ',A6,': IFAIL =',I6,/4X,'Check all subroutine calls',
     *       ' and array sizes. Seek expert help.')
99979 FORMAT (1X,'** A serious error has occurred in an internal call ',
     *       'to F11GBF: IREVCM =',I6,/4X,'Check all subroutine calls ',
     *       'and array sizes. Seek expert help.')
      END
