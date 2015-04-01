      SUBROUTINE F07FWZ(UPLO,N,A,LDA,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             ZLAUUM(UPLO,N,A,LDA,INFO)
C
C  Purpose
C  =======
C
C  ZLAUUM computes the product U * U' or L' * L, where the triangular
C  factor U or L is stored in the upper or lower triangular part of
C  the array A.
C
C  If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
C  overwriting the factor U in A.
C  If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
C  overwriting the factor L in A.
C
C  This is the blocked form of the algorithm, calling Level 3 BLAS.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the triangular factor stored in the array A
C          is upper or lower triangular:
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  N       (input) INTEGER
C          The order of the triangular factor U or L.  N >= 0.
C
C  A       (input/output) COMPLEX array, dimension (LDA,N)
C          On entry, the triangular factor U or L.
C          On exit, if UPLO = 'U', the upper triangle of A is
C          overwritten with the upper triangle of the product U * U';
C          if UPLO = 'L', the lower triangle of A is overwritten with
C          the lower triangle of the product L' * L.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
      COMPLEX*16        CONE
      PARAMETER         (CONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*)
C     .. Local Scalars ..
      INTEGER           I, IB, NB
      LOGICAL           UPPER
C     .. External Subroutines ..
      EXTERNAL          ZGEMM, ZHERK, ZTRMM, F06AAZ, F07FWY, F07ZAZ
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
         CALL F06AAZ('F07FWZ/ZLAUUM',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
C     Determine the block size for this environment.
C
      CALL F07ZAZ(1,'F07FWZ',NB,0)
      IF (NB.LE.1) THEN
         CALL F07FWY(UPLO,N,A,LDA,INFO)
      ELSE
         IF (UPPER) THEN
C
C           Compute the product U * U'.
C
            DO 20 I = 1, N, NB
               IB = MIN(NB,N-I+1)
               CALL ZTRMM('Right','Upper','Conjugate transpose',
     *                    'Non-unit',I-1,IB,CONE,A(I,I),LDA,A(1,I),LDA)
               CALL F07FWY('Upper',IB,A(I,I),LDA,INFO)
               IF (I+IB.LE.N) THEN
                  CALL ZGEMM('No transpose','Conjugate transpose',I-1,
     *                       IB,N-I-IB+1,CONE,A(1,I+IB),LDA,A(I,I+IB),
     *                       LDA,CONE,A(1,I),LDA)
                  CALL ZHERK('Upper','No transpose',IB,N-I-IB+1,ONE,
     *                       A(I,I+IB),LDA,ONE,A(I,I),LDA)
               END IF
   20       CONTINUE
         ELSE
C
C           Compute the product L' * L.
C
            DO 40 I = 1, N, NB
               IB = MIN(NB,N-I+1)
               CALL ZTRMM('Left','Lower','Conjugate transpose',
     *                    'Non-unit',IB,I-1,CONE,A(I,I),LDA,A(I,1),LDA)
               CALL F07FWY('Lower',IB,A(I,I),LDA,INFO)
               IF (I+IB.LE.N) THEN
                  CALL ZGEMM('Conjugate transpose','No transpose',IB,
     *                       I-1,N-I-IB+1,CONE,A(I+IB,I),LDA,A(I+IB,1),
     *                       LDA,CONE,A(I,1),LDA)
                  CALL ZHERK('Lower','Conjugate transpose',IB,N-I-IB+1,
     *                       ONE,A(I+IB,I),LDA,ONE,A(I,I),LDA)
               END IF
   40       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F07FWZ (ZLAUUM)
C
      END
