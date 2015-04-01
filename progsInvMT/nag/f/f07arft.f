      SUBROUTINE F07ARF(M,N,A,LDA,IPIV,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 16 REVISED. IER-1015 (JUN 1993).
C
C     .. Entry Points ..
      ENTRY             ZGETRF(M,N,A,LDA,IPIV,INFO)
C
C  Purpose
C  =======
C
C  ZGETRF computes an LU factorization of a general m-by-n matrix A
C  using partial pivoting with row interchanges.
C
C  The factorization has the form
C     A = P * L * U
C  where P is a permutation matrix, L is lower triangular with unit
C  diagonal elements (lower trapezoidal if m > n), and U is upper
C  triangular (upper trapezoidal if m < n).
C
C  This is the Level 3 BLAS version of the right-looking algorithm.
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
      COMPLEX*16        ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, M, N
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*)
      INTEGER           IPIV(*)
C     .. Local Scalars ..
      INTEGER           I, IINFO, J, JB, NB
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07ARY, F07ARZ, F07ZAZ, ZGEMM, ZTRSM
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
         CALL F06AAZ('F07ARF/ZGETRF',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
C
C     Determine the block size for this environment.
C
      CALL F07ZAZ(1,'F07ARF',NB,0)
      IF (NB.LE.1) THEN
C
C        Use unblocked code.
C
         CALL F07ARZ(M,N,A,LDA,IPIV,INFO)
      ELSE
C
C        Use blocked code.
C
         DO 40 J = 1, MIN(M,N), NB
            JB = MIN(MIN(M,N)-J+1,NB)
C
C           Factor diagonal and subdiagonal blocks and test for exact
C           singularity.
C
            CALL F07ARZ(M-J+1,JB,A(J,J),LDA,IPIV(J),IINFO)
C
C           Adjust INFO and the pivot indices.
C
            IF (INFO.EQ.0 .AND. IINFO.GT.0) INFO = IINFO + J - 1
            DO 20 I = J, MIN(M,J+JB-1)
               IPIV(I) = J - 1 + IPIV(I)
   20       CONTINUE
C
C           Apply interchanges to columns 1:J-1.
C
            CALL F07ARY(J-1,A,LDA,J,J+JB-1,IPIV,1)
C
            IF (J+JB.LE.N) THEN
C
C              Apply interchanges to columns J+JB:N.
C
               CALL F07ARY(N-J-JB+1,A(1,J+JB),LDA,J,J+JB-1,IPIV,1)
C
C              Compute block row of U.
C
               CALL ZTRSM('Left','Lower','No transpose','Unit',JB,
     *                    N-J-JB+1,ONE,A(J,J),LDA,A(J,J+JB),LDA)
               IF (J+JB.LE.M) THEN
C
C                 Update trailing submatrix.
C
                  CALL ZGEMM('No transpose','No transpose',M-J-JB+1,
     *                       N-J-JB+1,JB,-ONE,A(J+JB,J),LDA,A(J,J+JB),
     *                       LDA,ONE,A(J+JB,J+JB),LDA)
               END IF
            END IF
   40    CONTINUE
      END IF
      RETURN
C
C     End of F07ARF (ZGETRF)
C
      END
