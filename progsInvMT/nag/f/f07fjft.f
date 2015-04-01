      SUBROUTINE F07FJF(UPLO,N,A,LDA,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             DPOTRI(UPLO,N,A,LDA,INFO)
C
C  Purpose
C  =======
C
C  DPOTRI computes the inverse of a real symmetric positive definite
C  matrix A using the Cholesky factorization A = U'*U or A = L*L'
C  computed by F07FDF.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the factor stored in A is upper or lower
C          triangular.
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  A       (input/output) REAL array, dimension (LDA,N)
C          On entry, the triangular factor U or L from the Cholesky
C          factorization A = U'*U or A = L*L', as computed by F07FDF.
C          On exit, the upper or lower triangle of the (symmetric)
C          inverse of A, overwriting the input factor U or L.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C          > 0: if INFO = k, the (k,k) element of the factor U or L is
C               zero, and the inverse could not be computed.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*)
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07FJZ, F07TJF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF ( .NOT. (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
     *    .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l')) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -4
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07FJF/DPOTRI',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
C     Invert the triangular Cholesky factor U or L.
C
      CALL F07TJF(UPLO,'Non-unit',N,A,LDA,INFO)
      IF (INFO.GT.0) RETURN
C
C     Form inv(U)*inv(U)' or inv(L)'*inv(L).
C
      CALL F07FJZ(UPLO,N,A,LDA,INFO)
C
      RETURN
C
C     End of F07FJF (DPOTRI)
C
      END
