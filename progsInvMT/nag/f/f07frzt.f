      SUBROUTINE F07FRZ(UPLO,N,A,LDA,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             ZPOTF2(UPLO,N,A,LDA,INFO)
C
C  Purpose
C  =======
C
C  ZPOTF2 computes the Cholesky factorization of a complex Hermitian
C  positive definite matrix A.
C
C  The factorization has the form
C     A = U' * U ,  if UPLO = 'U', or
C     A = L  * L',  if UPLO = 'L',
C  where U is an upper triangular matrix and L is lower triangular.
C
C  This is the unblocked version of the algorithm, calling Level 2 BLAS.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the upper or lower triangular part of the
C          Hermitian matrix A is stored.
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  A       (input/output) COMPLEX array, dimension (LDA,N)
C          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
C          n by n upper triangular part of A contains the upper
C          triangular part of the matrix A, and the strictly lower
C          triangular part of A is not referenced.  If UPLO = 'L', the
C          leading n by n lower triangular part of A contains the lower
C          triangular part of the matrix A, and the strictly upper
C          triangular part of A is not referenced.
C
C          On exit, if INFO = 0, the factor U or L from the Cholesky
C          factorization A = U'*U  or A = L*L'.
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
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
      COMPLEX*16        CONE
      PARAMETER         (CONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AJJ
      INTEGER           J
      LOGICAL           UPPER
C     .. External Functions ..
      COMPLEX*16        ZDOTC
      EXTERNAL          ZDOTC
C     .. External Subroutines ..
      EXTERNAL          ZGEMV, ZDSCAL, F06AAZ, F07FRY
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, DBLE, SQRT
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
         CALL F06AAZ('F07FRZ/ZPOTF2',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
      IF (UPPER) THEN
C
C        Compute the Cholesky factorization A = U'*U.
C
         DO 20 J = 1, N
C
C           Compute U(J,J) and test for non-positive-definiteness.
C
            AJJ = DBLE(A(J,J)) - ZDOTC(J-1,A(1,J),1,A(1,J),1)
            IF (AJJ.LE.ZERO) THEN
               A(J,J) = AJJ
               GO TO 60
            END IF
            AJJ = SQRT(AJJ)
            A(J,J) = AJJ
C
C           Compute elements J+1:N of row J.
C
            IF (J.LT.N) THEN
               CALL F07FRY(J-1,A(1,J),1)
               CALL ZGEMV('Transpose',J-1,N-J,-CONE,A(1,J+1),LDA,A(1,J),
     *                    1,CONE,A(J,J+1),LDA)
               CALL F07FRY(J-1,A(1,J),1)
               CALL ZDSCAL(N-J,ONE/AJJ,A(J,J+1),LDA)
            END IF
   20    CONTINUE
      ELSE
C
C        Compute the Cholesky factorization A = L*L'.
C
         DO 40 J = 1, N
C
C           Compute L(J,J) and test for non-positive-definiteness.
C
            AJJ = DBLE(A(J,J)) - ZDOTC(J-1,A(J,1),LDA,A(J,1),LDA)
            IF (AJJ.LE.ZERO) THEN
               A(J,J) = AJJ
               GO TO 60
            END IF
            AJJ = SQRT(AJJ)
            A(J,J) = AJJ
C
C           Compute elements J+1:N of column J.
C
            IF (J.LT.N) THEN
               CALL F07FRY(J-1,A(J,1),LDA)
               CALL ZGEMV('No transpose',N-J,J-1,-CONE,A(J+1,1),LDA,
     *                    A(J,1),LDA,CONE,A(J+1,J),1)
               CALL F07FRY(J-1,A(J,1),LDA)
               CALL ZDSCAL(N-J,ONE/AJJ,A(J+1,J),1)
            END IF
   40    CONTINUE
      END IF
      GO TO 80
C
   60 CONTINUE
      INFO = J
C
   80 CONTINUE
      RETURN
C
C     End of F07FRZ (ZPOTF2)
C
      END
