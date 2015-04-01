      SUBROUTINE F11JAZ(N,NNZ,IROW,ICOL,SYM,IBAD,IERROR)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Purpose
C     =======
C     Checks the coordinate storage (CS) representation of an NxN sparse
C     matrix A. If SYM = .TRUE. it is assumed that A is symmetric and
C     only the lower triangular non-zero elements are stored.
C
C     Arguments
C     =========
C     N      (input) INTEGER
C            On entry, the order of the matrix A.
C            N >= 1.
C
C     NNZ    (input) INTEGER
C            On entry, the number of non-zero elements of A, or its
C            lower triangular part if SYM = .TRUE..
C            1 <= NNZ <= N**2      for SYM = .FALSE. and
C            1 <= NNZ <= N*(N+1)/2 for SYM = .TRUE.
C
C     IROW   (input) INTEGER array, dimension (*)
C     ICOL   (input) INTEGER array, dimension (*)
C            Note: the dimensions of IROW and ICOL must be at least NNZ.
C            On entry, the row and column indices of the non-zero
C            elements of the matrix A, or its lower triangular part if
C            SYM = .TRUE.. These elements should be without duplicate
C            entries, ordered by increasing row index and by increasing
C            column index within each row.
C            1 <= IROW(i) <= N,
C            1 <= ICOL(i) <= N       for SYM = .FALSE., and
C            1 <= ICOL(i) <= IROW(i) for SYM = .TRUE., for i = 1,...,N.
C            IROW(i-1) < IROW(i), or
C            IROW(i-1) = IROW(i) and ICOL(i-1) < ICOL(i), for
C            i = 2,3,...,NNZ.
C
C     SYM    (input) LOGICAL
C            On entry, indicates whether or not the matrix is stored in
C            symmetric form.
C                SYM = .TRUE.  => only the lower triangular non-zero
C                                 elements are stored.
C                SYM = .FALSE. => all non-zero elements are stored.
C
C     IBAD   (output) INTEGER
C            On exit, if IERROR = 4, 5, 6, or 7 IBAD holds the index of
C            the element which is out of range.
C
C     IERROR (output) INTEGER
C            On exit, the following values may occur:
C               IERROR = 0 => no error detected.
C               IERROR = 1 => N < 1.
C               IERROR = 2 => NNZ < 1.
C               IERROR = 3 => NNZ is too large.
C               IERROR = 4 => IROW(IBAD) is out of range.
C               IERROR = 5 => ICOL(IBAD) is out of range.
C               IERROR = 6 => the element IBAD is out of order.
C               IERROR = 7 => the element IBAD has the same row and
C                             column indices as a previous element.
C
C     ==================================================================
C
C     .. Scalar Arguments ..
      INTEGER           IBAD, IERROR, N, NNZ
      LOGICAL           SYM
C     .. Array Arguments ..
      INTEGER           ICOL(*), IROW(*)
C     .. Local Scalars ..
      INTEGER           I, ICI, ICIM1, IRI, IRIM1, MAXIND, MAXNNZ
C     .. Executable Statements ..
C
      IERROR = 0
      IBAD = 0
C
C     Check N and NNZ.
C
      IF (N.LT.1) THEN
         IERROR = 1
         GO TO 40
      END IF
C
      IF (NNZ.LT.1) THEN
         IERROR = 2
         GO TO 40
      END IF
C
      MAXNNZ = N*N
      IF (SYM) MAXNNZ = N*(N+1)/2
C
      IF (NNZ.GT.MAXNNZ) THEN
         IERROR = 3
         GO TO 40
      END IF
C
C     Check the location, ordering and uniqueness
C     of all non-zero matrix elements.
C
      DO 20 I = 1, NNZ
C
         IRI = IROW(I)
         ICI = ICOL(I)
         IRIM1 = IROW(I-1)
         ICIM1 = ICOL(I-1)
C
C         Location.
C
         IF (IRI.LT.1 .OR. IRI.GT.N) THEN
            IERROR = 4
            IBAD = I
            GO TO 40
         END IF
C
         MAXIND = N
         IF (SYM) MAXIND = IRI
C
         IF (ICI.LT.1 .OR. ICI.GT.MAXIND) THEN
            IERROR = 5
            IBAD = I
            GO TO 40
         END IF
C
         IF (I.GT.1) THEN
C
C             Ordering.
C
            IF (IRI.LT.IRIM1) THEN
               IERROR = 6
               IBAD = I
               GO TO 40
            END IF
C
            IF (IRI.EQ.IRIM1) THEN
C
               IF (ICI.LT.ICIM1) THEN
                  IERROR = 6
                  IBAD = I
                  GO TO 40
               END IF
C
C                 Uniqueness.
C
               IF (ICI.EQ.ICIM1) THEN
                  IERROR = 7
                  IBAD = I
                  GO TO 40
               END IF
C
            END IF
C
         END IF
C
   20 CONTINUE
C
   40 CONTINUE
C
      RETURN
      END
