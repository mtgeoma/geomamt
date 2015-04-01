      SUBROUTINE F07FWY(UPLO,N,A,LDA,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             ZLAUU2(UPLO,N,A,LDA,INFO)
C
C  Purpose
C  =======
C
C  ZLAUU2 computes the product U * U' or L' * L, where the triangular
C  factor U or L is stored in the upper or lower triangular part of
C  the array A.
C
C  If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
C  overwriting the factor U in A.
C  If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
C  overwriting the factor L in A.
C
C  This is the unblocked form of the algorithm, calling Level 2 BLAS.
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
      COMPLEX*16        ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AII
      INTEGER           I
      LOGICAL           UPPER
C     .. External Functions ..
      COMPLEX*16        ZDOTC
      EXTERNAL          ZDOTC
C     .. External Subroutines ..
      EXTERNAL          ZGEMV, ZDSCAL, F06AAZ, F07FRY
C     .. Intrinsic Functions ..
      INTRINSIC         DCMPLX, MAX, DBLE
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
         CALL F06AAZ('F07FWY/ZLAUU2',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
      IF (UPPER) THEN
C
C        Compute the product U * U'.
C
         DO 20 I = 1, N
            AII = A(I,I)
            IF (I.LT.N) THEN
               A(I,I) = AII*AII + DBLE(ZDOTC(N-I,A(I,I+1),LDA,A(I,I+1)
     *                  ,LDA))
               CALL F07FRY(N-I,A(I,I+1),LDA)
               CALL ZGEMV('No transpose',I-1,N-I,ONE,A(1,I+1),LDA,
     *                    A(I,I+1),LDA,DCMPLX(AII),A(1,I),1)
               CALL F07FRY(N-I,A(I,I+1),LDA)
            ELSE
               CALL ZDSCAL(I,AII,A(1,I),1)
            END IF
   20    CONTINUE
C
      ELSE
C
C        Compute the product L' * L.
C
         DO 40 I = 1, N
            AII = A(I,I)
            IF (I.LT.N) THEN
               A(I,I) = AII*AII + DBLE(ZDOTC(N-I,A(I+1,I),1,A(I+1,I),1))
               CALL F07FRY(I-1,A(I,1),LDA)
               CALL ZGEMV('Conjugate transpose',N-I,I-1,ONE,A(I+1,1),
     *                    LDA,A(I+1,I),1,DCMPLX(AII),A(I,1),LDA)
               CALL F07FRY(I-1,A(I,1),LDA)
            ELSE
               CALL ZDSCAL(I,AII,A(I,1),LDA)
            END IF
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F07FWY (ZLAUU2)
C
      END
