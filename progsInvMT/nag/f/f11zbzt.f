      SUBROUTINE F11ZBZ(N,NNZ,A,IROW,ICOL,ISTR,IWORK)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Purpose
C     =======
C     Sorts the non-zero elements of an NxN sparse matrix A into row
C     order and returns a pointer ISTR to the starting index of each
C     row in A.
C
C
C     Arguments
C     =========
C     N      (input) INTEGER
C            On entry, the order of the matrix A.
C
C     NNZ    (input) INTEGER
C            On entry, the number of non-zero elements stored.
C
C     A      (input/output) DOUBLE PRECISION array, dimension (NNZ)
C            On entry, the non-zero elements of the matrix A.
C            On exit, the non-zero elements in row order.
C
C     IROW   (input/output) INTEGER array, dimension (NNZ)
C            On entry, the row indices corresponding to the non-zero
C            elements stored in the array A.
C            On exit, the row indices in row order.
C
C     ICOL   (input/output) INTEGER array, dimension (NNZ)
C            On entry, the column indices corresponding to the non-zero
C            elements stored in the array A.
C            On exit, the column indices in row order.
C
C     ISTR   (output) INTEGER array, dimension (N+1)
C            On exit, ISTR(i) holds the start address in the arrays
C            A, IROW and ICOL of row i of the matrix, for i=1,2,...,N.
C            ISTR(N+1) holds NNZ+1.
C
C     IWORK  (workspace) INTEGER array, dimension (N)
C
C
C     Further Details
C     ===============
C     This routine sorts the non-zero matrix elements into row order,
C     i.e. IROW(i+1) >= IROW(i), but within each row the elements need
C     not be in column order. To achieve this the routine uses a simple
C     sorting algorithm which is more efficient for this purpose than a
C     call to M01DBF, requiring O(NNZ) operations as opposed to
C     O(NNZ*log(NNZ)).
C
C     ==================================================================
C
C     .. Scalar Arguments ..
      INTEGER           N, NNZ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NNZ)
      INTEGER           ICOL(NNZ), IROW(NNZ), ISTR(N+1), IWORK(N)
C     .. Local Scalars ..
      INTEGER           I, IFAIL, IRI, J
C     .. External Subroutines ..
      EXTERNAL          F06DBF, F06DFF, M01EAF, M01EBF
C     .. Executable Statements ..
C
C     Find number of elements in each row.
C
      CALL F06DBF(N,0,IWORK,1)
C
      DO 20 I = 1, NNZ
         IRI = IROW(I)
         IWORK(IRI) = IWORK(IRI) + 1
   20 CONTINUE
C
C     Set up pointer to start of each row.
C
      ISTR(1) = 1
      DO 40 I = 1, N
         ISTR(I+1) = ISTR(I) + IWORK(I)
   40 CONTINUE
C
C     Use IWORK as a temporary copy of ISTR.
C
      CALL F06DFF(N,ISTR,1,IWORK,1)
C
C     Find rank vector and store in IROW.
C
      DO 60 I = 1, NNZ
         J = IROW(I)
         IROW(I) = IWORK(J)
         IWORK(J) = IWORK(J) + 1
   60 CONTINUE
C
C     Reorder A and ICOL.
C
      IFAIL = 1
      CALL M01EAF(A,1,NNZ,IROW,IFAIL)
      IFAIL = 1
      CALL M01EBF(ICOL,1,NNZ,IROW,IFAIL)
C
C     Create ordered version of IROW.
C
      DO 100 I = 1, N
         DO 80 J = ISTR(I), ISTR(I+1) - 1
            IROW(J) = I
   80    CONTINUE
  100 CONTINUE
C
      RETURN
      END
