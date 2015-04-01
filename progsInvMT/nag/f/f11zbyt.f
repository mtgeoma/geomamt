      SUBROUTINE F11ZBY(N,NNZ,A,IROW,ICOL,ISTC,IWORK)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Purpose
C     =======
C     Sorts the non-zero elements of an NxN sparse matrix A into column
C     order and returns a pointer ISTC to the starting index of each
C     column in A.
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
C            On exit, the non-zero elements in column order.
C
C     IROW   (input/output) INTEGER array, dimension (NNZ)
C            On entry, the row indices corresponding to the non-zero
C            elements stored in the array A.
C            On exit, the row indices in column order.
C
C     ICOL   (input/output) INTEGER array, dimension (NNZ)
C            On entry, the column indices corresponding to the non-zero
C            elements stored in the array A.
C            On exit, the column indices in column order.
C
C     ISTC   (output) INTEGER array, dimension (N+1)
C            On exit, ISTC(i) holds the start address in the arrays
C            A, IROW and ICOL of column i of the matrix, for
C            i=1,2,...,N. ISTC(N+1) holds NNZ+1.
C
C     IWORK  (workspace) INTEGER array, dimension (N)
C
C
C     Further Details
C     ===============
C     This routine sorts the non-zero matrix elements into column order,
C     i.e. ICOL(i+1) >= ICOL(i), but within each column the elements
C     need not be in row order. To achieve this the routine uses a
C     simple sorting algorithm which is more efficient for this purpose
C     than a call to M01DBF, requiring O(NNZ) operations as opposed to
C     O(NNZ*log(NNZ)).
C
C     ==================================================================
C
C     .. Scalar Arguments ..
      INTEGER           N, NNZ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NNZ)
      INTEGER           ICOL(NNZ), IROW(NNZ), ISTC(N+1), IWORK(N)
C     .. Local Scalars ..
      INTEGER           I, ICI, IFAIL, IWJ, J
C     .. External Subroutines ..
      EXTERNAL          F06DBF, F06DFF, M01EAF, M01EBF
C     .. Executable Statements ..
C
C     Find number of elements in each column.
C
      CALL F06DBF(N,0,IWORK,1)
C
      DO 20 I = 1, NNZ
         ICI = ICOL(I)
         IWORK(ICI) = IWORK(ICI) + 1
   20 CONTINUE
C
C     Set up pointer to start of each column.
C
      ISTC(1) = 1
      DO 40 I = 1, N
         ISTC(I+1) = ISTC(I) + IWORK(I)
   40 CONTINUE
C
C     Use IWORK as a temporary copy of ISTC.
C
      CALL F06DFF(N,ISTC,1,IWORK,1)
C
C     Find rank vector and store in ICOL.
C
      DO 60 I = 1, NNZ
         J = ICOL(I)
         IWJ = IWORK(J)
         ICOL(I) = IWJ
         IWORK(J) = IWJ + 1
   60 CONTINUE
C
C     Reorder A and IROW.
C
      IFAIL = 1
      CALL M01EAF(A,1,NNZ,ICOL,IFAIL)
      IFAIL = 1
      CALL M01EBF(IROW,1,NNZ,ICOL,IFAIL)
C
C     Create ordered version of ICOL.
C
      DO 100 I = 1, N
         DO 80 J = ISTC(I), ISTC(I+1) - 1
            ICOL(J) = I
   80    CONTINUE
  100 CONTINUE
C
      RETURN
      END
