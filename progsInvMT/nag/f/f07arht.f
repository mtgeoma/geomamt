      SUBROUTINE F07ARH(M,N,A,LDA,PIV,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 16 REVISED. IER-1017 (JUN 1993).
C
C  Purpose
C  =======
C
C  F07ARH computes an LU factorization of a general m-by-n matrix A
C  using partial pivoting with row interchanges.
C
C  The factorization has the form
C     A = P * L * U
C  where P is a permutation matrix, L is lower triangular (lower
C  trapezoidal if m > n), and U is upper triangular with unit diagonal
C  elements (upper trapezoidal if m < n).
C
C  This is the Level 2 BLAS version of the right-looking algorithm.
C
C  Arguments
C  =========
C
C  M       (input) INTEGER
C          The number of rows of the matrix A.  M >= 0.
C
C  N       (input) INTEGER
C          The number of columns of the matrix A.  N >= 0.
C
C  A       (input/output) COMPLEX array, dimension (LDA,N)
C          On entry, the m by n matrix to be factored.
C          On exit, the factors L and U from the factorization
C          A = P*L*U; the unit diagonal elements of U are not stored.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,M).
C
C  PIV     (input/output) REAL array, dimension (M)
C          On entry, M scale factors for equilibrating the rows of A.
C          On exit, the pivot indices; for 1 <= i <= min(M,N), row i of
C          the matrix was interchanged with row PIV(i).
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C          > 0: if INFO = k, L(k,k) is exactly zero; the factorization
C               has not been completed.
C
C  This is a modified version of the LAPACK routine F07ARZ/ZGETF2, in
C  which the INTEGER array IPIV has been replaced by a REAL array PIV,
C  row-equilibration is used in the choice of pivot, U has unit diagonal
C  elements, and the routine exits immediately if singularity is
C  detected.
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ONE
      DOUBLE PRECISION  ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, M, N
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*)
      DOUBLE PRECISION  PIV(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  X, Y
      INTEGER           I, J, JP
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, ZGERU, ZSCAL, ZSWAP
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF (M.LT.0) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (LDA.LT.MAX(1,M)) THEN
         INFO = -4
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07ARH       ',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
C
      DO 40 J = 1, MIN(M,N)
C
C        Find pivot and test for singularity.
C
         JP = J
         X = ZERO
         DO 20 I = J, M
            Y = ABS(A(I,J))*PIV(I)
            IF (Y.GT.X) THEN
               JP = I
               X = Y
            END IF
   20    CONTINUE
         PIV(JP) = PIV(J)
         PIV(J) = JP
         IF (X.NE.ZERO) THEN
C
C           Apply interchange to columns 1:N.
C
            IF (JP.NE.J) CALL ZSWAP(N,A(J,1),LDA,A(JP,1),LDA)
C
C           Compute elements J+1:N of J-th row.
C
            IF (J.LT.N) CALL ZSCAL(N-J,ONE/A(J,J),A(J,J+1),LDA)
C
            IF (J.LT.MIN(M,N)) THEN
C
C              Update trailing submatrix.
C
               CALL ZGERU(M-J,N-J,-ONE,A(J+1,J),1,A(J,J+1),LDA,
     *                    A(J+1,J+1),LDA)
            END IF
C
         ELSE
C
C           If A( JP, J ) is zero, set INFO to indicate that a zero
C           pivot has been found.
C
            INFO = J
            RETURN
         END IF
   40 CONTINUE
      RETURN
C
C     End of F07ARH
C
      END
