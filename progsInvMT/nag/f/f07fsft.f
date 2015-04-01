      SUBROUTINE F07FSF(UPLO,N,NRHS,A,LDA,B,LDB,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 16 REVISED. IER-1022 (JUN 1993).
C     .. Entry Points ..
      ENTRY             ZPOTRS(UPLO,N,NRHS,A,LDA,B,LDB,INFO)
C
C  Purpose
C  =======
C
C  ZPOTRS solves a system of linear equations A*X = B with a Hermitian
C  positive definite matrix A using the Cholesky factorization A = U'*U
C  or A = L*L' computed by F07FRF.
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
C  A       (input) COMPLEX array, dimension (LDA,N)
C          The triangular factor U or L from the Cholesky factorization
C          A = U'*U or A = L*L', as computed by F07FRF.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
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
C     .. Parameters ..
      COMPLEX*16        ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDB, N, NRHS
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), B(LDB,*)
C     .. Local Scalars ..
      LOGICAL           UPPER
C     .. External Subroutines ..
      EXTERNAL          ZTRSM, ZTRSV, F06AAZ
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
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -5
      ELSE IF (LDB.LT.MAX(1,N)) THEN
         INFO = -7
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07FSF/ZPOTRS',-INFO)
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
C        Solve U'*Y = B, overwriting B with Y, and then solve U*X = Y,
C        overwriting Y with X.
C
         IF (NRHS.EQ.1) THEN
            CALL ZTRSV('Upper','Conjugate transpose','Non-unit',N,A,LDA,
     *                 B,1)
            CALL ZTRSV('Upper','No transpose','Non-unit',N,A,LDA,B,1)
         ELSE
            CALL ZTRSM('Left','Upper','Conjugate transpose','Non-unit',
     *                 N,NRHS,ONE,A,LDA,B,LDB)
            CALL ZTRSM('Left','Upper','No transpose','Non-unit',N,NRHS,
     *                 ONE,A,LDA,B,LDB)
         END IF
      ELSE
C
C        Solve A*X = B where A = L*L'.
C
C        Solve L*Y = B, overwriting B with Y, and then solve L'*X = Y,
C        overwriting Y with X.
C
         IF (NRHS.EQ.1) THEN
            CALL ZTRSV('Lower','No transpose','Non-unit',N,A,LDA,B,1)
            CALL ZTRSV('Lower','Conjugate transpose','Non-unit',N,A,LDA,
     *                 B,1)
         ELSE
            CALL ZTRSM('Left','Lower','No transpose','Non-unit',N,NRHS,
     *                 ONE,A,LDA,B,LDB)
            CALL ZTRSM('Left','Lower','Conjugate transpose','Non-unit',
     *                 N,NRHS,ONE,A,LDA,B,LDB)
         END IF
      END IF
C
      RETURN
C
C     End of F07FSF (ZPOTRS)
C
      END
