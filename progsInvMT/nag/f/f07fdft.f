      SUBROUTINE F07FDF(UPLO,N,A,LDA,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             DPOTRF(UPLO,N,A,LDA,INFO)
C
C  Purpose
C  =======
C
C  DPOTRF computes the Cholesky factorization of a real symmetric
C  positive definite matrix A.
C
C  The factorization has the form
C     A = U' * U ,  if UPLO = 'U', or
C     A = L  * L',  if UPLO = 'L',
C  where U is an upper triangular matrix and L is lower triangular.
C
C  This is the block version of the algorithm, calling Level 3 BLAS.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the upper or lower triangular part of the
C          symmetric matrix A is stored.
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  A       (input/output) REAL array, dimension (LDA,N)
C          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
C          n by n upper triangular part of A contains the upper
C          triangular part of the matrix A, and the strictly lower
C          triangular part of A is not referenced.  If UPLO = 'L', the
C          leading n by n lower triangular part of A contains the lower
C          triangular part of the matrix A, and the strictly upper
C          triangular part of A is not referenced.
C
C          On exit, if INFO = 0, the factor U or L from the Cholesky
C          factorization A = U'*U or A = L*L'.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C          > 0: if INFO = k, the leading minor of order k is not
C               positive definite, and the factorization could not be
C               completed.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*)
C     .. Local Scalars ..
      INTEGER           J, JB, NB
      LOGICAL           UPPER
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07FDZ, F07ZAZ, DGEMM, DSYRK, DTRSM
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
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
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -4
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07FDF/DPOTRF',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
C     Determine the block size for this environment.
C
      CALL F07ZAZ(1,'F07FDF',NB,0)
      IF (NB.LE.1) THEN
C
C        Use unblocked code.
C
         CALL F07FDZ(UPLO,N,A,LDA,INFO)
      ELSE
C
C        Use blocked code.
C
         IF (UPPER) THEN
C
C           Compute the Cholesky factorization A = U'*U.
C
            DO 20 J = 1, N, NB
C
C              Update and factorize the current diagonal block and test
C              for non-positive-definiteness.
C
               JB = MIN(NB,N-J+1)
               CALL DSYRK('Upper','Transpose',JB,J-1,-ONE,A(1,J),LDA,
     *                    ONE,A(J,J),LDA)
               CALL F07FDZ('Upper',JB,A(J,J),LDA,INFO)
               IF (INFO.NE.0) GO TO 60
               IF (J+JB.LE.N) THEN
C
C                 Compute the current block row.
C
                  CALL DGEMM('Transpose','No transpose',JB,N-J-JB+1,J-1,
     *                       -ONE,A(1,J),LDA,A(1,J+JB),LDA,ONE,A(J,J+JB)
     *                       ,LDA)
                  CALL DTRSM('Left','Upper','Transpose','Non-unit',JB,
     *                       N-J-JB+1,ONE,A(J,J),LDA,A(J,J+JB),LDA)
               END IF
   20       CONTINUE
C
         ELSE
C
C           Compute the Cholesky factorization A = L*L'.
C
            DO 40 J = 1, N, NB
C
C              Update and factorize the current diagonal block and test
C              for non-positive-definiteness.
C
               JB = MIN(NB,N-J+1)
               CALL DSYRK('Lower','No transpose',JB,J-1,-ONE,A(J,1),LDA,
     *                    ONE,A(J,J),LDA)
               CALL F07FDZ('Lower',JB,A(J,J),LDA,INFO)
               IF (INFO.NE.0) GO TO 60
               IF (J+JB.LE.N) THEN
C
C                 Compute the current block column.
C
                  CALL DGEMM('No transpose','Transpose',N-J-JB+1,JB,J-1,
     *                       -ONE,A(J+JB,1),LDA,A(J,1),LDA,ONE,A(J+JB,J)
     *                       ,LDA)
                  CALL DTRSM('Right','Lower','Transpose','Non-unit',
     *                       N-J-JB+1,JB,ONE,A(J,J),LDA,A(J+JB,J),LDA)
               END IF
   40       CONTINUE
         END IF
      END IF
      GO TO 80
C
   60 CONTINUE
      INFO = INFO + J - 1
C
   80 CONTINUE
      RETURN
C
C     End of F07FDF (DPOTRF)
C
      END
