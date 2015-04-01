      SUBROUTINE F11ZBF(N,NNZ,A,IROW,ICOL,DUP,ZERO,ISTR,IWORK,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Purpose
C     =======
C     Sorts the non-zero elements of a real NxN sparse matrix A into
C     increasing row order and increasing column order within each row.
C     Removes or sums values of any non-zero elements with duplicate
C     row and column indices. Optionally removes any zero elements.
C     Returns a pointer ISTR to the starting index of each row in A.
C
C     Arguments
C     =========
C     N      (input) INTEGER
C            On entry, the order of the matrix A.
C            N >= 1.
C
C     NNZ    (input/output) INTEGER
C            On entry, the number of lower triangular non-zero elements
C            supplied.
C            On exit, the number of lower triangular non-zero elements
C            with unique row and column indices.
C            0 <= NNZ.
C
C     A      (input/output) DOUBLE PRECISION array, dimension (*)
C            Note: the dimension of A must be at least NNZ.
C            On entry, the non-zero elements of the lower triangular
C            part of the matrix A. These may be in any order and there
C            may be multiple non-zero elements with the same row and
C            column indices.
C            On exit, the lower triangular non-zero elements ordered by
C            increasing row index, and by increasing column index within
C            each row. Each non-zero element has a unique row and column
C            index.
C
C     IROW   (input/output) INTEGER array, dimension (*)
C            Note: the dimension of IROW must be at least NNZ.
C            On entry, the row indices corresponding to the non-zero
C            elements supplied in the array A.
C            1 <= IROW(i) <= N, for i = 1,2,...,N.
C            On exit, the row indices corresponding to the non-zero
C            elements returned in the array A.
C
C     ICOL   (input/output) INTEGER array, dimension (*)
C            Note: the dimension of ICOL must be at least NNZ.
C            On entry, the column indices corresponding to the non-zero
C            elements supplied in the array A.
C            1 <= IROW(i) <= ICOL(i), for i = 1,2,...,N.
C            On exit, the column indices corresponding to the non-zero
C            elements returned in the array A.
C
C     DUP    (input) CHARACTER*1
C            On entry, indicates how non-zero elements with duplicate
C            row and column indices are to be treated,
C               DUP = 'R', or 'r' ==> duplicates are removed;
C               DUP = 'S', or 's' ==> duplicates are summed.
C            DUP = 'R', 'r', 'S', or 's'.
C
C     ZERO   (input) CHARACTER*1
C            On entry, indicates how elements with zero coefficient
C            values are to be treated,
C               ZERO = 'R', or 'r' ==> zeros are removed;
C               ZERO = 'K', or 'k' ==> zeros are kept.
C            DUP = 'R', 'r', 'K', or 'k'.
C
C     ISTR   (output) INTEGER array, dimension (N+1)
C            On exit, ISTR(i) holds the start address in the arrays
C            A, IROW and ICOL of row i of the matrix, for i=1,2,...,N.
C            ISTR(N+1) holds NNZ+1.
C
C     IWORK  (workspace) INTEGER array, dimension (N)
C
C     IFAIL  (input/output) INTEGER
C            On entry, IFAIL must be -1, 0, or 1.
C            On exit, the following values may occur:
C               IFAIL = 0 => no error detected.
C               IFAIL = 1 => N < 1, or
C                            NNZ < 0, or
C                            DUP invalid.
C               IFAIL = 2 => An element of IROW, or ICOL is out of
C                            range.
C
C     Further Details
C     ===============
C     This routine uses a simple sorting algorithm which is more
C     efficient for the given task than a call to M01DBF, requiring
C     O(NNZ) operations as opposed to O(NNZ*log(NNZ)).
C
C     ==================================================================
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F11ZBF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, N, NNZ
      CHARACTER         DUP, ZERO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(*)
      INTEGER           ICOL(*), IROW(*), ISTR(N+1), IWORK(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  DIDUP
      INTEGER           I, IDUP, IERR, IRI, IZERO, J, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06DBF, F11ZBY, F11ZBZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
C     Check the SCS representation of A.
C
      IERR = 0
C
      IF (N.LT.1) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99999) N
         GO TO 120
      END IF
C
      IF (NNZ.LT.0) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99998) NNZ
         GO TO 120
      END IF
C
      DO 20 I = 1, NNZ
         IF (IROW(I).LT.1 .OR. IROW(I).GT.N) THEN
            IERR = 2
            NREC = 2
            WRITE (REC,FMT=99997) I, IROW(I), N
            GO TO 120
         END IF
         IF (ICOL(I).LT.1 .OR. ICOL(I).GT.IROW(I)) THEN
            IERR = 2
            NREC = 2
            WRITE (REC,FMT=99996) I, ICOL(I), IROW(I)
            GO TO 120
         END IF
   20 CONTINUE
C
C     Check DUP.
C
      IDUP = -1
      IF (DUP.EQ.'R' .OR. DUP.EQ.'r') IDUP = 0
      IF (DUP.EQ.'S' .OR. DUP.EQ.'s') IDUP = 1
      IF (IDUP.EQ.-1) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99995) DUP
         GO TO 120
      END IF
C
C     Check ZERO.
C
      IZERO = -1
      IF (ZERO.EQ.'R' .OR. ZERO.EQ.'r') IZERO = 0
      IF (ZERO.EQ.'K' .OR. ZERO.EQ.'k') IZERO = 1
      IF (IZERO.EQ.-1) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99994) ZERO
         GO TO 120
      END IF
C
      IF (NNZ.GT.1) THEN
C
C         Sort non-zero elements into row order, with column order
C         within each row.
C
         CALL F11ZBY(N,NNZ,A,IROW,ICOL,ISTR,IWORK)
         CALL F11ZBZ(N,NNZ,A,IROW,ICOL,ISTR,IWORK)
C
C         Remove or sum non-zero elements with duplicate row and column
C         indices.
C
         J = 1
         DIDUP = DBLE(IDUP)
         DO 40 I = 2, NNZ
            IF (IROW(I).EQ.IROW(J) .AND. ICOL(I).EQ.ICOL(J)) THEN
               A(J) = A(J) + DIDUP*A(I)
            ELSE
               J = J + 1
               A(J) = A(I)
               IROW(J) = IROW(I)
               ICOL(J) = ICOL(I)
            END IF
   40    CONTINUE
C
         NNZ = J
C
      END IF
C
      IF (IZERO.EQ.0) THEN
C
C         Remove any elements with zero values.
C
         J = 0
         DO 60 I = 1, NNZ
            IF (A(I).NE.0.D0) THEN
               J = J + 1
               A(J) = A(I)
               IROW(J) = IROW(I)
               ICOL(J) = ICOL(I)
            END IF
   60    CONTINUE
C
         NNZ = J
C
      END IF
C
C     Find number of elements in each row.
C
      CALL F06DBF(N,0,IWORK,1)
C
      DO 80 I = 1, NNZ
         IRI = IROW(I)
         IWORK(IRI) = IWORK(IRI) + 1
   80 CONTINUE
C
C     Set up pointer to start of each row.
C
      ISTR(1) = 1
      DO 100 I = 1, N
         ISTR(I+1) = ISTR(I) + IWORK(I)
  100 CONTINUE
C
  120 CONTINUE
C
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, N .lt. 1: N =',I16,'.')
99998 FORMAT (1X,'** On entry, NNZ .lt. 0: NNZ =',I16,'.')
99997 FORMAT (1X,'** On entry, IROW(I) .lt. 1 or IROW(I) .gt. N:',/4X,
     *       'I = ',I16,', IROW(I) = ',I16,', N = ',I16,'.')
99996 FORMAT (1X,'** On entry, ICOL(I) .lt. 1 or ICOL(I) .gt. IROW(I):',
     *       /4X,'I =',I16,', ICOL(I) =',I16,', IROW(I) =',I16,'.')
99995 FORMAT (1X,'** On entry, DUP .ne. ''R'' or ''S'': DUP = ''',A,
     *       '''.')
99994 FORMAT (1X,'** On entry, ZERO .ne. ''R'' or ''K'': ZERO = ''',A,
     *       '''.')
      END
