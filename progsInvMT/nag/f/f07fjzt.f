      SUBROUTINE F07FJZ(UPLO,N,A,LDA,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             DLAUUM(UPLO,N,A,LDA,INFO)
C
C  Purpose
C  =======
C
C  DLAUUM computes the product U * U' or L' * L, where the triangular
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
C  A       (input/output) REAL array, dimension (LDA,N)
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
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*)
C     .. Local Scalars ..
      INTEGER           I, IB, NB
      LOGICAL           UPPER
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07FJY, F07ZAZ, DGEMM, DSYRK, DTRMM
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
         CALL F06AAZ('F07FJZ/DLAUUM',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
C     Determine the block size for this environment.
C
      CALL F07ZAZ(1,'F07FJZ',NB,0)
      IF (NB.LE.1) THEN
         CALL F07FJY(UPLO,N,A,LDA,INFO)
      ELSE
         IF (UPPER) THEN
C
C           Compute the product U * U'.
C
            DO 20 I = 1, N, NB
               IB = MIN(NB,N-I+1)
               CALL DTRMM('Right','Upper','Transpose','Non-unit',I-1,IB,
     *                    ONE,A(I,I),LDA,A(1,I),LDA)
               CALL F07FJY('Upper',IB,A(I,I),LDA,INFO)
               IF (I+IB.LE.N) THEN
                  CALL DGEMM('No transpose','Transpose',I-1,IB,N-I-IB+1,
     *                       ONE,A(1,I+IB),LDA,A(I,I+IB),LDA,ONE,A(1,I),
     *                       LDA)
                  CALL DSYRK('Upper','No transpose',IB,N-I-IB+1,ONE,
     *                       A(I,I+IB),LDA,ONE,A(I,I),LDA)
               END IF
   20       CONTINUE
         ELSE
C
C           Compute the product L' * L.
C
            DO 40 I = 1, N, NB
               IB = MIN(NB,N-I+1)
               CALL DTRMM('Left','Lower','Transpose','Non-unit',IB,I-1,
     *                    ONE,A(I,I),LDA,A(I,1),LDA)
               CALL F07FJY('Lower',IB,A(I,I),LDA,INFO)
               IF (I+IB.LE.N) THEN
                  CALL DGEMM('Transpose','No transpose',IB,I-1,N-I-IB+1,
     *                       ONE,A(I+IB,I),LDA,A(I+IB,1),LDA,ONE,A(I,1),
     *                       LDA)
                  CALL DSYRK('Lower','Transpose',IB,N-I-IB+1,ONE,
     *                       A(I+IB,I),LDA,ONE,A(I,I),LDA)
               END IF
   40       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F07FJZ (DLAUUM)
C
      END
