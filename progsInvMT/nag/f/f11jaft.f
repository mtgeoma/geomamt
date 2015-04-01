      SUBROUTINE F11JAF(N,NNZ,A,LA,IROW,ICOL,LFILL,DTOL,MIC,DSCALE,
     *                  PSTRAT,IPIV,ISTR,NNZC,NPIVM,IWORK,LIWORK,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Purpose
C     =======
C     F11JAF computes an incomplete Cholesky factorization:
C
C               A  =  M + R
C
C                            T T
C               M  =  P L D L P
C
C     of a real NxN symmetric sparse matrix A, of arbitrary sparsity
C     pattern, stored in symmetric coordinate storage (SCS) format. It
C     is designed for positive-definite matrices, but may also work for
C     some mildly indefinite cases. In the above decomposition L is
C     lower triangular with unit diagonal elements, D is diagonal, P is
C     a permutation matrix and R is a remainder matrix.
C
C     The matrix M is returned in terms of the SCS representation of
C     the lower triangular matrix
C
C                          -1
C               C  =  L + D  - I .
C
C     The amount of fill in the decomposition may be controlled either
C     by specifying the maximum level of fill LFILL, or by supplying a
C     drop tolerance DTOL below which fill-in elements are omitted.
C
C     By setting MIC = 'M' the factorization may be modified in order
C     to preserve the row sums of the original matrix, and diagonal
C     elements may be scaled prior to factorization to ensure that the
C     matrix M is positive-definite. For certain mesh-based problems
C     setting MIC = 'M' and choosing the scaling parameter DSCALE
C     appropriately can reduce the order of magnitude of the condition
C     number of the preconditioned matrix as a function of the mesh
C     steplength.
C
C     The argument PSTRAT defines the pivoting strategy to be used. The
C     available options are no pivoting, user-defined diagonal pivoting,
C     and diagonal pivoting based on the Markowitz startegy, aimed at
C     minimizing fill-in.
C
C     Arguments
C     =========
C     N      (input) INTEGER
C            On entry, the order of the matrix A.
C            N >= 1.
C
C     NNZ    (input) INTEGER
C            On entry, the number of non-zero elements in the lower
C            triangular part of A.
C            1 <= NNZ <= N*(N+1)/2.
C
C     A      (input/output) DOUBLE PRECISION array, dimension (LA)
C            On entry, the non-zero elements in the lower triangular
C            part of the matrix A, ordered by increasing row index
C            and by increasing column index within each row. Multiple
C            entries for the same row and column indices are not
C            allowed.
C            On exit, the first NNZ elements of A are unchanged and the
C            next NNZC contain the non-zero elements of C.
C
C     LA     (input) INTEGER
C            On entry, the dimension of the arrays A, IROW and ICOL
C            as declared in the (sub)program from which F11JAF is
C            called. LA must be an estimate of NNZ + NNZC.
C            LA >= 2*NNZ.
C
C     IROW   (input/output) INTEGER array, dimension (LA)
C     ICOL   (input/output) INTEGER array, dimension (LA)
C            On entry, the row and column indices corresponding to the
C            non-zero elements given in the array A.
C            IROW and ICOL must satisfy the following constraints:
C            1 <= IROW(i) <= N, and 1 <= ICOL(i) <= IROW(i), for
C            i = 1,2,...,NNZ.
C            IROW(i-1) < IROW(i), or
C            IROW(i-1) = IROW(i) and ICOL(i-1) < ICOL(i), for
C            i = 2,3,...,NNZ.
C            On exit, the row and column indices corresponding to the
C            non-zero elements returned in the array A.
C
C     LFILL  (input) INTEGER
C            On entry, if LFILL >= 0 it gives the required level of
C            fill of the decomposition. If LFILL < 0 it indicates that
C            DTOL will be used to control the fill instead.
C
C     DTOL   (input) DOUBLE PRECISION
C            On entry, if LFILL >= 0 DTOL is not used. Otherwise DTOL
C            gives the drop tolerance of the decomposition. A fill
C            element a   in row i and column j will be dropped if
C                     ij
C
C                |a  | < DTOL*SQRT(| a   * a   |).
C                  ij                 ii    jj
C
C            DTOL >= 0.0
C
C     MIC    (input) CHARACTER*1
C            On entry, indicates whether or not the factorization should
C            be modified to preserve row sums:
C               MIC = 'M'  => the factorization is modified (MIC);
C               MIC = 'N'  => the factorization is not modified.
C            MIC = 'M', or 'N'.
C
C     DSCALE (input) DOUBLE PRECISION
C            On entry, the diagonal scaling parameter. Diagonal elements
C            are multiplied by the factor (1 + DSCALE) at the start of
C            the factorization. This can be used to ensure that the
C            matrix M generated is positive-definite. For certain
C            mesh-based problems setting MIC = 'M' and choosing DSCALE
C            appropriately can reduce the order of magnitude of the
C            condition number of the preconditioned matrix as a function
C            of the mesh steplength.
C
C     PSTRAT (input) CHARACTER*1
C            On entry, specifies the pivoting strategy to be adopted.
C               PSTRAT = 'N'  => no pivoting.
C               PSTRAT = 'U'  => user defined pivoting specified by the
C                                input value of IPIV.
C               PSTRAT = 'M'  => diagonal pivoting aimed at minimizing
C                                fill, based on the Markowitz strategy.
C            PSTRAT = 'N', 'U', or 'M'.
C
C     IPIV   (input/output) INTEGER array, dimension (N)
C            On entry, if PSTRAT = 'U' or 'u' then IPIV(i) must specify
C            the row index of the diagonal element used as a pivot at
C            elimination stage i. Otherwise IPIV need not be defined.
C            On exit, the pivot indices. If IPIV(i) = j then the
C            diagonal element in row j was used as the pivot at
C            elimination stage i.
C
C     ISTR   (output) INTEGER array, dimension (N+1)
C            On exit, ISTR(i) gives the starting address in the arrays
C            A, IROW and ICOL, of row i of the matrix C. ISTR(N+1)
C            holds NNZ+NNZC+1.
C
C     NNZC   (output) INTEGER
C            On exit, the number of non-zero elements in the lower
C            triangular matrix C.
C
C     NPIVM  (output) INTEGER
C            On exit, the number of pivots which were modified during
C            the factorization to ensure that M is positive-definite.
C            If A is an M-matrix then no pivot modifications will be
C            required, and NPIVM = 0. For other cases the quality of
C            the preconditioner will generally depend on the returned
C            value of NPIVM. If NPIVM is large the preconditioner may
C            not be satisfactory. In this case it may be advantageous
C            to call F11JAF again with an increased value of DSCALE.
C
C     IWORK  (workspace) INTEGER array, dimension (LIWORK)
C
C     LIWORK (input) INTEGER
C            On entry, the dimension of the array IWORK as declared in
C            the (sub)program from which F11JAF is called.
C            LIWORK >= 2*LA - 3*NNZ + 7*N + 1, for LFILL >= 0, or
C            LIWORK >=   LA -   NNZ + 7*N + 1, for LFILL <  0.
C
C     IFAIL  (input/output) INTEGER
C            On entry, IFAIL must be -1, 0, or 1.
C            On exit, the following values may occur:
C               IFAIL = 0 => no error detected.
C               IFAIL = 1 => N < 1, or
C                            NNZ < 1, or
C                            NNZ > N*(N+1)/2, or
C                            LA < 2*NNZ, or
C                            DTOL < 0.0, or
C                            MIC not equal to 'M', or 'N', or
C                            PSTRAT not equal to 'N', 'U', or 'M', or
C                            LIWORK < 2*LA-3*NNZ+7*N+1, LFILL >= 0, or
C                            LIWORK < LA-NNZ+7*N+1, LFILL <  0.
C               IFAIL = 2 => the arrays IROW and ICOL fail to satisfy
C                            the following constraints:
C                            1 <= IROW(i) <= N,
C                            1 <= ICOL(i) <= IROW(i), i = 1,2,...,NNZ.
C                            IROW(i-1) < IROW(i), or
C                            IROW(i-1) = IROW(i) and ICOL(i-1) < ICOL(i)
C                            for i = 2,3,...,NNZ.
C               IFAIL = 3 => PSTRAT = 'U' or 'u' and IPIV does not
C                            represent a valid permutation of the
C                            integers on [1,N]. An element of IPIV
C                            is out of range or repeated.
C               IFAIL = 4 => LA is too small, resulting in insufficient
C                            space for fill-in elements.
C                            The decomposition has been terminated
C                            before completion. Either increase LA,
C                            or reduce the fill by setting PSTRAT = 'M',
C                            reducing LFILL, or increasing DTOL.
C               IFAIL = 5 => a serious error has occurred in an internal
C                            call to F11ZBF.
C
C     Further Details
C     ===============
C     The time taken for a call to F11JAF is roughly proportional to
C     NNZC*NNZC/N.
C
C     This routine is a driver for F11JAY.
C
C     ==================================================================
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F11JAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DSCALE, DTOL
      INTEGER           IFAIL, LA, LFILL, LIWORK, N, NNZ, NNZC, NPIVM
      CHARACTER         MIC, PSTRAT
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LA)
      INTEGER           ICOL(LA), IPIV(N), IROW(LA), ISTR(N+1),
     *                  IWORK(LIWORK)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA
      INTEGER           I, IBAD, IC, IDLEVF, IERR, IERROR, IW1, IW2,
     *                  IW3, IW4, IW5, IW6, IW7, IW8, IW9, IWE, MAXF,
     *                  NREC
      LOGICAL           SYM
C     .. Local Arrays ..
      CHARACTER*80      REC(4)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, F06DFF, F11JAY, F11JAZ, M01ZBF
C     .. Executable Statements ..
C
C     Check input arguments LA, DTOL, MIC and PSTRAT.
C
      IERR = 0
      NNZC = 0
      IF (LA.LT.2*NNZ) THEN
         IERR = 1
         NREC = 2
         WRITE (REC,FMT=99999) LA, NNZ
         GO TO 40
      END IF
      IF (LFILL.LT.0) THEN
         IF (DTOL.LT.0.D0) THEN
            IERR = 1
            NREC = 1
            WRITE (REC,FMT=99998) DTOL
            GO TO 40
         END IF
      END IF
      IF (MIC.EQ.'N' .OR. MIC.EQ.'n') THEN
         ALPHA = 0.D0
      ELSE IF (MIC.EQ.'M' .OR. MIC.EQ.'m') THEN
         ALPHA = 1.D0
      ELSE
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99997) MIC
         GO TO 40
      END IF
      IF (PSTRAT.NE.'N' .AND. PSTRAT.NE.'n' .AND. PSTRAT.NE.'U' .AND.
     *    PSTRAT.NE.'u' .AND. PSTRAT.NE.'M' .AND. PSTRAT.NE.'m') THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99996) PSTRAT
         GO TO 40
      END IF
C
C     Check matrix.
C
      SYM = .TRUE.
      CALL F11JAZ(N,NNZ,IROW,ICOL,SYM,IBAD,IERROR)
C
      IF (IERROR.EQ.1) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99995) N
         GO TO 40
      END IF
C
      IF (IERROR.EQ.2) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99994) NNZ
         GO TO 40
      END IF
C
      IF (IERROR.EQ.3) THEN
         IERR = 1
         NREC = 2
         WRITE (REC,FMT=99993) NNZ, N
         GO TO 40
      END IF
C
      IF (IERROR.EQ.4) THEN
         IERR = 2
         NREC = 2
         WRITE (REC,FMT=99992) IBAD, IROW(IBAD), N
         GO TO 40
      END IF
C
      IF (IERROR.EQ.5) THEN
         IERR = 2
         NREC = 2
         WRITE (REC,FMT=99991) IBAD, ICOL(IBAD), IROW(IBAD)
         GO TO 40
      END IF
C
      IF (IERROR.EQ.6) THEN
         IERR = 2
         NREC = 2
         WRITE (REC,FMT=99990) IBAD
         GO TO 40
      END IF
C
      IF (IERROR.EQ.7) THEN
         IERR = 2
         NREC = 2
         WRITE (REC,FMT=99989) IBAD
         GO TO 40
      END IF
C
C     Check IPIV if it is user-supplied.
C
      IF (PSTRAT.EQ.'U' .OR. PSTRAT.EQ.'u') THEN
         IERROR = 1
         CALL M01ZBF(IPIV,1,N,IERROR)
         IF (IERROR.EQ.2) THEN
            IERR = 3
            NREC = 1
            WRITE (REC,FMT=99988)
            GO TO 40
         END IF
         IF (IERROR.EQ.3) THEN
            IERR = 3
            NREC = 1
            WRITE (REC,FMT=99987)
            GO TO 40
         END IF
      END IF
C
C     Determine starting index in the arrays A, IROW and ICOL for
C     the matrix C. Copy the matrix A to this location.
C
      IC = NNZ + 1
      CALL DCOPY(NNZ,A,1,A(IC),1)
      CALL F06DFF(NNZ,IROW,1,IROW(IC),1)
      CALL F06DFF(NNZ,ICOL,1,ICOL(IC),1)
C
C     Split up workspace.
C
      MAXF = LA - 2*NNZ
      IDLEVF = 1
      IF (LFILL.GE.0) IDLEVF = MAXF + 1
C
      IW1 = 1
      IW2 = IW1 + N
      IW3 = IW2 + IDLEVF
      IW4 = IW3 + N
      IW5 = IW4 + NNZ + MAXF
      IW6 = IW5 + N
      IW7 = IW6 + N
      IW8 = IW7 + N
      IW9 = IW8 + N
      IWE = IW9 + N
C
      IF (LIWORK.LT.IWE-1) THEN
         IERR = 1
         NREC = 2
         WRITE (REC,FMT=99986) LIWORK, IWE - 1
         GO TO 40
      END IF
C
C     Factorization.
C
      CALL F11JAY(N,NNZ,NNZC,MAXF,A(IC),IROW(IC),ICOL(IC),LFILL,DTOL,
     *            DSCALE,ALPHA,PSTRAT,IPIV,NPIVM,IWORK(IW1),IWORK(IW2),
     *            IDLEVF,IWORK(IW3),IWORK(IW4),ISTR,IWORK(IW5),
     *            IWORK(IW6),IWORK(IW7),IWORK(IW8),IWORK(IW9),IERROR)
C
C     Shift elements of ISTR.
C
      DO 20 I = 1, N + 1
         ISTR(I) = ISTR(I) + NNZ
   20 CONTINUE
C
      IF (IERROR.EQ.1) THEN
         IERR = 4
         NREC = 4
         WRITE (REC,FMT=99985)
         GO TO 40
      END IF
C
      IF (IERROR.EQ.2) THEN
         IERR = 5
         NREC = 2
         WRITE (REC,FMT=99984)
         GO TO 40
      END IF
C
   40 CONTINUE
C
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, LA .lt. 2*NNZ:',/4X,'LA =',I16,' NNZ =',
     *       I16,'.')
99998 FORMAT (1X,'** On entry, DTOL .lt. 0.0: DTOL =',E12.4,'.')
99997 FORMAT (1X,'** On entry, MIC .ne. ''M'' or ''N'': MIC = ''',A,
     *       '''.')
99996 FORMAT (1X,'** On entry, PSTRAT .ne. ''N'', ''U'' or ''M'':',' P',
     *       'STRAT = ''',A,'''.')
99995 FORMAT (1X,'** On entry, N .lt. 1: N =',I16,'.')
99994 FORMAT (1X,'** On entry, NNZ .lt. 1: NNZ =',I16,'.')
99993 FORMAT (1X,'** On entry, NNZ .gt. N*(N+1)/2:',/4X,'NNZ =',I16,
     *       ' N =',I16,'.')
99992 FORMAT (1X,'** On entry, IROW(I) .lt. 1 or IROW(I) .gt. N:',/4X,
     *       'I = ',I16,', IROW(I) = ',I16,', N = ',I16,'.')
99991 FORMAT (1X,'** On entry, ICOL(I) .lt. 1 or ICOL(I) .gt. IROW(I):',
     *       /4X,'I =',I16,', ICOL(I) =',I16,', IROW(I) =',I16,'.')
99990 FORMAT (1X,'** On entry, A(I) is out of order:',/4X,'I =',I16,'.')
99989 FORMAT (1X,'** On entry, the location (IROW(I), ICOL(I)) is a ',
     *       'duplicate:',/4X,'I =',I16,'.')
99988 FORMAT (1X,'** On entry, a user-supplied value of IPIV lies',' o',
     *       'utside the range [1,N].')
99987 FORMAT (1X,'** On entry, a user-supplied value of IPIV is',' rep',
     *       'eated.')
99986 FORMAT (1X,'** On entry, LIWORK is too small: LIWORK =',I16,'.',
     *       /4X,'Minimum required value of LIWORK =',I16,'.')
99985 FORMAT (1X,'** The number of non-zero entries in the',' decompos',
     *       'ition is too large.',/4X,'The decomposition',' has been ',
     *       'terminated before completion.',/4X,'Either',' increase L',
     *       'A, or reduce the fill by setting PSTRAT',' = ''M'',',/4X,
     *       'reducing LFILL, or increasing DTOL.')
99984 FORMAT (1X,'** A serious error has occurred in an internal call',
     *       ' to F11ZBF.',/4X,'Check all subroutine calls and array',
     *       ' sizes. Seek expert help.')
      END
