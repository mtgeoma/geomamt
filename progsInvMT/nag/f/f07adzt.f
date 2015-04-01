      SUBROUTINE F07ADZ(M,N,A,LDA,IPIV,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 16 REVISED. IER-1012 (JUN 1993).
C     ENTRY             DGETF2(M,N,A,LDA,IPIV,INFO)
C
C  Purpose
C  =======
C
C  DGETF2 computes an LU factorization of a general m-by-n matrix A
C  using partial pivoting with row interchanges.
C
C  The factorization has the form
C     A = P * L * U
C  where P is a permutation matrix, L is lower triangular with unit
C  diagonal elements (lower trapezoidal if m > n), and U is upper
C  triangular (upper trapezoidal if m < n).
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
C  A       (input/output) REAL array, dimension (LDA,N)
C          On entry, the m by n matrix to be factored.
C          On exit, the factors L and U from the factorization
C          A = P*L*U; the unit diagonal elements of L are not stored.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,M).
C
C  IPIV    (output) INTEGER array, dimension (min(M,N))
C          The pivot indices; for 1 <= i <= min(M,N), row i of the
C          matrix was interchanged with row IPIV(i).
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
C               has been completed, but the factor U is exactly
C               singular, and division by zero will occur if it is used
C               to solve a system of equations.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*)
      INTEGER           IPIV(*)
C     .. Local Scalars ..
      INTEGER           J, JP
C     .. External Functions ..
      INTEGER           IDAMAX
      EXTERNAL          IDAMAX
C     .. External Subroutines ..
      EXTERNAL          DGER, DSCAL, DSWAP, F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
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
         CALL F06AAZ('F07ADZ/DGETF2',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
C
      DO 20 J = 1, MIN(M,N)
C
C        Find pivot and test for singularity.
C
         JP = J - 1 + IDAMAX(M-J+1,A(J,J),1)
         IPIV(J) = JP
         IF (A(JP,J).NE.ZERO) THEN
C
C           Apply interchange to columns 1:N.
C
            IF (JP.NE.J) CALL DSWAP(N,A(J,1),LDA,A(JP,1),LDA)
C
C           Compute elements J+1:M of J-th column.
C
            IF (J.LT.M) CALL DSCAL(M-J,ONE/A(J,J),A(J+1,J),1)
C
         ELSE IF (INFO.EQ.0) THEN
C
C           If A( JP, J ) is zero, set INFO to indicate that a zero
C           pivot has been found.
C
            INFO = J
         END IF
C
         IF (J.LT.MIN(M,N)) THEN
C
C           Update trailing submatrix.
C
            CALL DGER(M-J,N-J,-ONE,A(J+1,J),1,A(J,J+1),LDA,A(J+1,J+1),
     *                LDA)
         END IF
   20 CONTINUE
      RETURN
C
C     End of F07ADZ (DGETF2)
C
      END
