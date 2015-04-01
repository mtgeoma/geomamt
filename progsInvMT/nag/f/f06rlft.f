      DOUBLE PRECISION FUNCTION F06RLF(NORM,UPLO,DIAG,N,K,AB,LDAB,WORK)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     REAL                 DLANTB
C     ENTRY                DLANTB(NORM,UPLO,DIAG,N,K,AB,LDAB,WORK)
C
C  Purpose
C  =======
C
C  DLANTB  returns the value of the one norm,  or the Frobenius norm, or
C  the  infinity norm,  or the element of  largest absolute value  of an
C  n by n triangular band matrix A,  with ( k + 1 ) diagonals.
C
C  Description
C  ===========
C
C  DLANTB returns the value
C
C     DLANTB = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
C          Specifies the value to be returned in DLANTB as described
C          above.
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the matrix A is upper or lower triangular.
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  DIAG    (input) CHARACTER*1
C          Specifies whether or not the matrix A is unit triangular.
C          = 'N':  Non-unit triangular
C          = 'U':  Unit triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.  When N = 0, DLANTB is
C          set to zero.
C
C  K       (input) INTEGER
C          The number of super-diagonals of the matrix A if UPLO = 'U',
C          or the number of sub-diagonals of the matrix A if UPLO = 'L'.
C          K >= 0.
C
C  AB      (input) REAL array, dimension (LDAB,N)
C          The upper or lower triangular band matrix A, stored in the
C          first k+1 rows of AB.  The j-th column of A is stored
C          in the j-th column of the array AB as follows:
C          if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for max(1,j-k)<=i<=j;
C          if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=min(n,j+k).
C          Note that when DIAG = 'U', the elements of the array AB
C          corresponding to the diagonal elements of the matrix A are
C          not referenced, but are assumed to be one.
C
C  LDAB    (input) INTEGER
C          The leading dimension of the array AB.  LDAB >= K+1.
C
C  WORK    (workspace) REAL array, dimension (LWORK),
C          where LWORK >= N when NORM = 'I'; otherwise, WORK is not
C          referenced.
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
      CHARACTER                        DIAG, NORM, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION                 AB(LDAB,*), WORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION                 SCALE, SUM, VALUE
      INTEGER                          I, J, L
      LOGICAL                          UDIAG
C     .. External Subroutines ..
      EXTERNAL                         F06FJF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, MAX, MIN, SQRT
C     .. Executable Statements ..
C
      IF (N.EQ.0) THEN
         VALUE = ZERO
      ELSE IF ((NORM.EQ.'M' .OR. NORM.EQ.'m')) THEN
C
C        Find max(abs(A(i,j))).
C
         IF ((DIAG.EQ.'U' .OR. DIAG.EQ.'u')) THEN
            VALUE = ONE
            IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
               DO 40 J = 1, N
                  DO 20 I = MAX(K+2-J,1), K
                     VALUE = MAX(VALUE,ABS(AB(I,J)))
   20             CONTINUE
   40          CONTINUE
            ELSE
               DO 80 J = 1, N
                  DO 60 I = 2, MIN(N+1-J,K+1)
                     VALUE = MAX(VALUE,ABS(AB(I,J)))
   60             CONTINUE
   80          CONTINUE
            END IF
         ELSE
            VALUE = ZERO
            IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
               DO 120 J = 1, N
                  DO 100 I = MAX(K+2-J,1), K + 1
                     VALUE = MAX(VALUE,ABS(AB(I,J)))
  100             CONTINUE
  120          CONTINUE
            ELSE
               DO 160 J = 1, N
                  DO 140 I = 1, MIN(N+1-J,K+1)
                     VALUE = MAX(VALUE,ABS(AB(I,J)))
  140             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE IF (((NORM.EQ.'O' .OR. NORM.EQ.'o')) .OR. (NORM.EQ.'1')) THEN
C
C        Find norm1(A).
C
         VALUE = ZERO
         UDIAG = (DIAG.EQ.'U' .OR. DIAG.EQ.'u')
         IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
            DO 220 J = 1, N
               IF (UDIAG) THEN
                  SUM = ONE
                  DO 180 I = MAX(K+2-J,1), K
                     SUM = SUM + ABS(AB(I,J))
  180             CONTINUE
               ELSE
                  SUM = ZERO
                  DO 200 I = MAX(K+2-J,1), K + 1
                     SUM = SUM + ABS(AB(I,J))
  200             CONTINUE
               END IF
               VALUE = MAX(VALUE,SUM)
  220       CONTINUE
         ELSE
            DO 280 J = 1, N
               IF (UDIAG) THEN
                  SUM = ONE
                  DO 240 I = 2, MIN(N+1-J,K+1)
                     SUM = SUM + ABS(AB(I,J))
  240             CONTINUE
               ELSE
                  SUM = ZERO
                  DO 260 I = 1, MIN(N+1-J,K+1)
                     SUM = SUM + ABS(AB(I,J))
  260             CONTINUE
               END IF
               VALUE = MAX(VALUE,SUM)
  280       CONTINUE
         END IF
      ELSE IF ((NORM.EQ.'I' .OR. NORM.EQ.'i')) THEN
C
C        Find normI(A).
C
         VALUE = ZERO
         IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
            IF ((DIAG.EQ.'U' .OR. DIAG.EQ.'u')) THEN
               DO 300 I = 1, N
                  WORK(I) = ONE
  300          CONTINUE
               DO 340 J = 1, N
                  L = K + 1 - J
                  DO 320 I = MAX(1,J-K), J - 1
                     WORK(I) = WORK(I) + ABS(AB(L+I,J))
  320             CONTINUE
  340          CONTINUE
            ELSE
               DO 360 I = 1, N
                  WORK(I) = ZERO
  360          CONTINUE
               DO 400 J = 1, N
                  L = K + 1 - J
                  DO 380 I = MAX(1,J-K), J
                     WORK(I) = WORK(I) + ABS(AB(L+I,J))
  380             CONTINUE
  400          CONTINUE
            END IF
         ELSE
            IF ((DIAG.EQ.'U' .OR. DIAG.EQ.'u')) THEN
               DO 420 I = 1, N
                  WORK(I) = ONE
  420          CONTINUE
               DO 460 J = 1, N
                  L = 1 - J
                  DO 440 I = J + 1, MIN(N,J+K)
                     WORK(I) = WORK(I) + ABS(AB(L+I,J))
  440             CONTINUE
  460          CONTINUE
            ELSE
               DO 480 I = 1, N
                  WORK(I) = ZERO
  480          CONTINUE
               DO 520 J = 1, N
                  L = 1 - J
                  DO 500 I = J, MIN(N,J+K)
                     WORK(I) = WORK(I) + ABS(AB(L+I,J))
  500             CONTINUE
  520          CONTINUE
            END IF
         END IF
         DO 540 I = 1, N
            VALUE = MAX(VALUE,WORK(I))
  540    CONTINUE
      ELSE IF (((NORM.EQ.'F' .OR. NORM.EQ.'f'))
     *         .OR. ((NORM.EQ.'E' .OR. NORM.EQ.'e'))) THEN
C
C        Find normF(A).
C
         IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
            IF ((DIAG.EQ.'U' .OR. DIAG.EQ.'u')) THEN
               SCALE = ONE
               SUM = N
               IF (K.GT.0) THEN
                  DO 560 J = 2, N
                     CALL F06FJF(MIN(J-1,K),AB(MAX(K+2-J,1),J),1,SCALE,
     *                           SUM)
  560             CONTINUE
               END IF
            ELSE
               SCALE = ZERO
               SUM = ONE
               DO 580 J = 1, N
                  CALL F06FJF(MIN(J,K+1),AB(MAX(K+2-J,1),J),1,SCALE,SUM)
  580          CONTINUE
            END IF
         ELSE
            IF ((DIAG.EQ.'U' .OR. DIAG.EQ.'u')) THEN
               SCALE = ONE
               SUM = N
               IF (K.GT.0) THEN
                  DO 600 J = 1, N - 1
                     CALL F06FJF(MIN(N-J,K),AB(2,J),1,SCALE,SUM)
  600             CONTINUE
               END IF
            ELSE
               SCALE = ZERO
               SUM = ONE
               DO 620 J = 1, N
                  CALL F06FJF(MIN(N-J+1,K+1),AB(1,J),1,SCALE,SUM)
  620          CONTINUE
            END IF
         END IF
         VALUE = SCALE*SQRT(SUM)
      END IF
C
      F06RLF = VALUE
      RETURN
C
C     End of F06RLF (DLANTB)
C
      END
