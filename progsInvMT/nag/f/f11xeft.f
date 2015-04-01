      SUBROUTINE F11XEF(N,NNZ,A,IROW,ICOL,CHECK,X,Y,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Purpose
C     =======
C     Computes the matrix-vector product
C
C              y  =  Ax
C
C     where A is an NxN symmetric sparse matrix, of arbitrary
C     sparsity pattern, stored in symmetric coordinate storage (SCS)
C     format.
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
C     A      (input) DOUBLE PRECISION array, dimension (NNZ)
C            On entry, the non-zero elements of the lower triangular
C            part of the matrix A, ordered by increasing row index,
C            and by increasing column index within each row. Multiple
C            entries for the same row and column indices are not
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
C     CHECK  (input) CHARACTER*1
C            On entry, specifies whether or not the input data should
C            be checked.
C               CHECK = 'C'  => checks are carried on the values of
C                               N, NNZ, IROW and ICOL.
C               CHECK = 'N'  => no checks are carried out.
C            CHECK = 'C' or 'N'.
C
C     X      (input) DOUBLE PRECISION array, dimension (N)
C            On entry, the vector x.
C
C     Y      (output) DOUBLE PRECISION array, dimension (N)
C            On exit, the vector y.
C
C     IFAIL  (input/output) INTEGER
C            On entry, IFAIL must be -1, 0, or 1.
C            On exit, the following values may occur:
C               IFAIL = 0 => no error detected.
C               IFAIL = 1 => CHECK .ne. 'C' or 'N'.
C               IFAIL = 2 => N < 1, or
C                            NNZ < 1, or
C                            NNZ > N*(N+1)/2.
C               IFAIL = 3 => the arrays IROW and ICOL fail to satisfy
C                            the following constraints:
C                            1 <= IROW(i) <= N,
C                            1 <= ICOL(i) <= IROW(i), i = 1,2,...,NNZ.
C                            IROW(i-1) < IROW(i), or
C                            IROW(i-1) = IROW(i) and ICOL(i-1) < ICOL(i)
C                            for i = 2,3,...,NNZ.
C
C     ==================================================================
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F11XEF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, N, NNZ
      CHARACTER         CHECK
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NNZ), X(N), Y(N)
      INTEGER           ICOL(NNZ), IROW(NNZ)
C     .. Local Scalars ..
      INTEGER           I, IBAD, ICI, IERR, IERROR, IRI, NREC
      LOGICAL           CK, SYM
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06FBF, F11JAZ
C     .. Executable Statements ..
C
      IERR = 0
C
C     Check CHECK.
C
      IF (CHECK.EQ.'C' .OR. CHECK.EQ.'c') THEN
         CK = .TRUE.
      ELSE IF (CHECK.EQ.'N' .OR. CHECK.EQ.'n') THEN
         CK = .FALSE.
      ELSE
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99999) CHECK
         GO TO 40
      END IF
C
      IF (CK) THEN
C
C         Check matrix.
C
         SYM = .TRUE.
         CALL F11JAZ(N,NNZ,IROW,ICOL,SYM,IBAD,IERROR)
C
         IF (IERROR.EQ.1) THEN
            IERR = 2
            NREC = 1
            WRITE (REC,FMT=99998) N
            GO TO 40
         END IF
C
         IF (IERROR.EQ.2) THEN
            IERR = 2
            NREC = 1
            WRITE (REC,FMT=99997) NNZ
            GO TO 40
         END IF
C
         IF (IERROR.EQ.3) THEN
            IERR = 2
            NREC = 2
            WRITE (REC,FMT=99996) NNZ, N
            GO TO 40
         END IF
C
         IF (IERROR.EQ.4) THEN
            IERR = 3
            NREC = 2
            WRITE (REC,FMT=99995) IBAD, IROW(IBAD), N
            GO TO 40
         END IF
C
         IF (IERROR.EQ.5) THEN
            IERR = 3
            NREC = 2
            WRITE (REC,FMT=99994) IBAD, ICOL(IBAD), IROW(IBAD)
            GO TO 40
         END IF
C
         IF (IERROR.EQ.6) THEN
            IERR = 3
            NREC = 2
            WRITE (REC,FMT=99993) IBAD
            GO TO 40
         END IF
C
         IF (IERROR.EQ.7) THEN
            IERR = 3
            NREC = 2
            WRITE (REC,FMT=99992) IBAD
            GO TO 40
         END IF
C
      END IF
C
C     Compute matrix vector product.
C
      CALL F06FBF(N,0.D0,Y,1)
C
      DO 20 I = 1, NNZ
         IRI = IROW(I)
         ICI = ICOL(I)
         Y(IRI) = Y(IRI) + A(I)*X(ICI)
         IF (IRI.NE.ICI) THEN
            Y(ICI) = Y(ICI) + A(I)*X(IRI)
         END IF
   20 CONTINUE
C
   40 CONTINUE
C
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, CHECK .ne. ''C'' or ''N'': CHECK = ''',A,
     *       '''.')
99998 FORMAT (1X,'** On entry, N .lt. 1: N =',I16,'.')
99997 FORMAT (1X,'** On entry, NNZ .lt. 1: NNZ =',I16,'.')
99996 FORMAT (1X,'** On entry, NNZ .gt. N*(N+1)/2:',/4X,'NNZ =',I16,
     *       ' N =',I16,'.')
99995 FORMAT (1X,'** On entry, IROW(I) .lt. 1 or IROW(I) .gt. N:',/4X,
     *       'I = ',I16,', IROW(I) = ',I16,', N = ',I16,'.')
99994 FORMAT (1X,'** On entry, ICOL(I) .lt. 1 or ICOL(I) .gt. IROW(I):',
     *       /4X,'I =',I16,', ICOL(I) =',I16,', IROW(I) =',I16,'.')
99993 FORMAT (1X,'** On entry, A(I) is out of order:',/4X,'I =',I16,'.')
99992 FORMAT (1X,'** On entry, the location (IROW(I), ICOL(I)) is a ',
     *       'duplicate:',/4X,'I =',I16,'.')
      END
