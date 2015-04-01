      DOUBLE PRECISION FUNCTION F06RCF(NORM,UPLO,N,A,LDA,WORK)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     REAL                 DLANSY
C     ENTRY                DLANSY(NORM,UPLO,N,A,LDA,WORK)
C
C  Purpose
C  =======
C
C  DLANSY  returns the value of the one norm,  or the Frobenius norm, or
C  the  infinity norm,  or the  element of  largest absolute value  of a
C  real symmetric matrix A.
C
C  Description
C  ===========
C
C  DLANSY returns the value
C
C     DLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm'
C              (
C              ( norm1(A),         NORM = '1', 'O' or 'o'
C              (
C              ( normI(A),         NORM = 'I' or 'i'
C              (
C              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
C
C  where  norm1  denotes the  one norm of a matrix (maximum column sum),
C  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
C  normF  denotes the  Frobenius norm of a matrix (square root of sum of
C  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
C
C  Arguments
C  =========
C
C  NORM    (input) CHARACTER*1
C          Specifies the value to be returned in DLANSY as described
C          above.
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the upper or lower triangular part of the
C          symmetric matrix A is to be referenced.
C          = 'U':  Upper triangular part of A is referenced
C          = 'L':  Lower triangular part of A is referenced
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.  When N = 0, DLANSY is
C          set to zero.
C
C  A       (input) REAL array, dimension (LDA,N)
C          The symmetric matrix A.  If UPLO = 'U', the leading n by n
C          upper triangular part of A contains the upper triangular part
C          of the matrix A, and the strictly lower triangular part of A
C          is not referenced.  If UPLO = 'L', the leading n by n lower
C          triangular part of A contains the lower triangular part of
C          the matrix A, and the strictly upper triangular part of A is
C          not referenced.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(N,1).
C
C  WORK    (workspace) REAL array, dimension (LWORK),
C          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
C          WORK is not referenced.
C
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C     .. Parameters ..
      DOUBLE PRECISION                 ONE, ZERO
      PARAMETER                        (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER                          LDA, N
      CHARACTER                        NORM, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION                 A(LDA,*), WORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION                 ABSA, SCALE, SUM, VALUE
      INTEGER                          I, J
C     .. External Subroutines ..
      EXTERNAL                         F06FJF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, MAX, SQRT
C     .. Executable Statements ..
C
      IF (N.EQ.0) THEN
         VALUE = ZERO
      ELSE IF ((NORM.EQ.'M' .OR. NORM.EQ.'m')) THEN
C
C        Find max(abs(A(i,j))).
C
         VALUE = ZERO
         IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
            DO 40 J = 1, N
               DO 20 I = 1, J
                  VALUE = MAX(VALUE,ABS(A(I,J)))
   20          CONTINUE
   40       CONTINUE
         ELSE
            DO 80 J = 1, N
               DO 60 I = J, N
                  VALUE = MAX(VALUE,ABS(A(I,J)))
   60          CONTINUE
   80       CONTINUE
         END IF
      ELSE IF (((NORM.EQ.'I' .OR. NORM.EQ.'i'))
     *         .OR. ((NORM.EQ.'O' .OR. NORM.EQ.'o')) .OR. (NORM.EQ.'1'))
     *         THEN
C
C        Find normI(A) ( = norm1(A), since A is symmetric).
C
         VALUE = ZERO
         IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
            DO 120 J = 1, N
               SUM = ZERO
               DO 100 I = 1, J - 1
                  ABSA = ABS(A(I,J))
                  SUM = SUM + ABSA
                  WORK(I) = WORK(I) + ABSA
  100          CONTINUE
               WORK(J) = SUM + ABS(A(J,J))
  120       CONTINUE
            DO 140 I = 1, N
               VALUE = MAX(VALUE,WORK(I))
  140       CONTINUE
         ELSE
            DO 160 I = 1, N
               WORK(I) = ZERO
  160       CONTINUE
            DO 200 J = 1, N
               SUM = WORK(J) + ABS(A(J,J))
               DO 180 I = J + 1, N
                  ABSA = ABS(A(I,J))
                  SUM = SUM + ABSA
                  WORK(I) = WORK(I) + ABSA
  180          CONTINUE
               VALUE = MAX(VALUE,SUM)
  200       CONTINUE
         END IF
      ELSE IF (((NORM.EQ.'F' .OR. NORM.EQ.'f'))
     *         .OR. ((NORM.EQ.'E' .OR. NORM.EQ.'e'))) THEN
C
C        Find normF(A).
C
         SCALE = ZERO
         SUM = ONE
         IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
            DO 220 J = 2, N
               CALL F06FJF(J-1,A(1,J),1,SCALE,SUM)
  220       CONTINUE
         ELSE
            DO 240 J = 1, N - 1
               CALL F06FJF(N-J,A(J+1,J),1,SCALE,SUM)
  240       CONTINUE
         END IF
         SUM = 2*SUM
         CALL F06FJF(N,A,LDA+1,SCALE,SUM)
         VALUE = SCALE*SQRT(SUM)
      END IF
C
      F06RCF = VALUE
      RETURN
C
C     End of F06RCF (DLANSY)
C
      END
