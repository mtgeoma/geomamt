      DOUBLE PRECISION FUNCTION F06UEF(NORM,UPLO,N,K,AB,LDAB,WORK)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     REAL                 ZLANHB
C     ENTRY                ZLANHB(NORM,UPLO,N,K,AB,LDAB,WORK)
C
C  Purpose
C  =======
C
C  ZLANHB  returns the value of the one norm,  or the Frobenius norm, or
C  the  infinity norm,  or the element of  largest absolute value  of an
C  n by n hermitian band matrix A,  with k super-diagonals.
C
C  Description
C  ===========
C
C  ZLANHB returns the value
C
C     ZLANHB = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
C          Specifies the value to be returned in ZLANHB as described
C          above.
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the upper or lower triangular part of the
C          band matrix A is supplied.
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.  When N = 0, ZLANHB is
C          set to zero.
C
C  K       (input) INTEGER
C          The number of super-diagonals or sub-diagonals of the
C          band matrix A.  K >= 0.
C
C  AB      (input) COMPLEX array, dimension (LDAB,N)
C          The upper or lower triangle of the hermitian band matrix A,
C          stored in the first K+1 rows of AB.  The j-th column of A is
C          stored in the j-th column of the array AB as follows:
C          if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for max(1,j-k)<=i<=j;
C          if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=min(n,j+k).
C          Note that the imaginary parts of the diagonal elements need
C          not be set and are assumed to be zero.
C
C  LDAB    (input) INTEGER
C          The leading dimension of the array AB.  LDAB >= K+1.
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
      INTEGER                          K, LDAB, N
      CHARACTER                        NORM, UPLO
C     .. Array Arguments ..
      COMPLEX*16                       AB(LDAB,*)
      DOUBLE PRECISION                 WORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION                 ABSA, SCALE, SUM, VALUE
      INTEGER                          I, J, L
C     .. External Subroutines ..
      EXTERNAL                         F06KJF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, MAX, MIN, DBLE, SQRT
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
               DO 20 I = MAX(K+2-J,1), K
                  VALUE = MAX(VALUE,ABS(AB(I,J)))
   20          CONTINUE
               VALUE = MAX(VALUE,ABS(DBLE(AB(K+1,J))))
   40       CONTINUE
         ELSE
            DO 80 J = 1, N
               VALUE = MAX(VALUE,ABS(DBLE(AB(1,J))))
               DO 60 I = 2, MIN(N+1-J,K+1)
                  VALUE = MAX(VALUE,ABS(AB(I,J)))
   60          CONTINUE
   80       CONTINUE
         END IF
      ELSE IF (((NORM.EQ.'I' .OR. NORM.EQ.'i'))
     *         .OR. ((NORM.EQ.'O' .OR. NORM.EQ.'o')) .OR. (NORM.EQ.'1'))
     *         THEN
C
C        Find normI(A) ( = norm1(A), since A is hermitian).
C
         VALUE = ZERO
         IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
            DO 120 J = 1, N
               SUM = ZERO
               L = K + 1 - J
               DO 100 I = MAX(1,J-K), J - 1
                  ABSA = ABS(AB(L+I,J))
                  SUM = SUM + ABSA
                  WORK(I) = WORK(I) + ABSA
  100          CONTINUE
               WORK(J) = SUM + ABS(DBLE(AB(K+1,J)))
  120       CONTINUE
            DO 140 I = 1, N
               VALUE = MAX(VALUE,WORK(I))
  140       CONTINUE
         ELSE
            DO 160 I = 1, N
               WORK(I) = ZERO
  160       CONTINUE
            DO 200 J = 1, N
               SUM = WORK(J) + ABS(DBLE(AB(1,J)))
               L = 1 - J
               DO 180 I = J + 1, MIN(N,J+K)
                  ABSA = ABS(AB(L+I,J))
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
         IF (K.GT.0) THEN
            IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
               DO 220 J = 2, N
                  CALL F06KJF(MIN(J-1,K),AB(MAX(K+2-J,1),J),1,SCALE,SUM)
  220          CONTINUE
               L = K + 1
            ELSE
               DO 240 J = 1, N - 1
                  CALL F06KJF(MIN(N-J,K),AB(2,J),1,SCALE,SUM)
  240          CONTINUE
               L = 1
            END IF
            SUM = 2*SUM
         ELSE
            L = 1
         END IF
         DO 260 J = 1, N
            IF (DBLE(AB(L,J)).NE.ZERO) THEN
               ABSA = ABS(DBLE(AB(L,J)))
               IF (SCALE.LT.ABSA) THEN
                  SUM = ONE + SUM*(SCALE/ABSA)**2
                  SCALE = ABSA
               ELSE
                  SUM = SUM + (ABSA/SCALE)**2
               END IF
            END IF
  260    CONTINUE
         VALUE = SCALE*SQRT(SUM)
      END IF
C
      F06UEF = VALUE
      RETURN
C
C     End of F06UEF (ZLANHB)
C
      END
