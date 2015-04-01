      SUBROUTINE F07GSF(UPLO,N,NRHS,AP,B,LDB,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             ZPPTRS(UPLO,N,NRHS,AP,B,LDB,INFO)
C
C  Purpose
C  =======
C
C  ZPPTRS solves a system of linear equations A*X = B with a Hermitian
C  positive definite matrix A in packed storage using the Cholesky
C  factorization A = U'*U or A = L*L' computed by F07GRF.
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
C  NRHS    (input) INTEGER
C          The number of right hand sides, i.e., the number of columns
C          of the matrix B.  NRHS >= 0.
C
C  AP      (input) COMPLEX array, dimension (N*(N+1)/2)
C          The triangular factor U or L from the Cholesky factorization
C          A = U'*U or A = L*L', packed columnwise in a linear array.
C          The j-th column of U or L is stored in the array AP as
C          follows:
C          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;
C          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n.
C
C  B       (input/output) COMPLEX array, dimension (LDB,NRHS)
C          On entry, the right hand side vectors B for the system of
C          linear equations.
C          On exit, the solution vectors, X.
C
C  LDB     (input) INTEGER
C          The leading dimension of the array B.  LDB >= max(1,N).
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Scalar Arguments ..
      INTEGER           INFO, LDB, N, NRHS
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        AP(*), B(LDB,*)
C     .. Local Scalars ..
      INTEGER           I
      LOGICAL           UPPER
C     .. External Subroutines ..
      EXTERNAL          ZTPSV, F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
      IF ( .NOT. UPPER .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l')) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (NRHS.LT.0) THEN
         INFO = -3
      ELSE IF (LDB.LT.MAX(1,N)) THEN
         INFO = -6
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07GSF/ZPPTRS',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0 .OR. NRHS.EQ.0) RETURN
C
      IF (UPPER) THEN
C
C        Solve A*X = B where A = U'*U.
C
         DO 20 I = 1, NRHS
C
C           Solve U'*X = B, overwriting B with X.
C
            CALL ZTPSV('Upper','Conjugate transpose','Non-unit',N,AP,
     *                 B(1,I),1)
C
C           Solve U*X = B, overwriting B with X.
C
            CALL ZTPSV('Upper','No transpose','Non-unit',N,AP,B(1,I),1)
   20    CONTINUE
      ELSE
C
C        Solve A*X = B where A = L*L'.
C
         DO 40 I = 1, NRHS
C
C           Solve L*Y = B, overwriting B with X.
C
            CALL ZTPSV('Lower','No transpose','Non-unit',N,AP,B(1,I),1)
C
C           Solve L'*X = Y, overwriting B with X.
C
            CALL ZTPSV('Lower','Conjugate transpose','Non-unit',N,AP,
     *                 B(1,I),1)
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F07GSF (ZPPTRS)
C
      END
