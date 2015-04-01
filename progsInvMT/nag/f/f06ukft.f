      DOUBLE PRECISION FUNCTION F06UKF(NORM,UPLO,DIAG,N,AP,WORK)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     REAL                 ZLANTP
C     ENTRY                ZLANTP(NORM,UPLO,DIAG,N,AP,WORK)
C
C  Purpose
C  =======
C
C  ZLANTP  returns the value of the one norm,  or the Frobenius norm, or
C  the  infinity norm,  or the  element of  largest absolute value  of a
C  triangular matrix A, supplied in packed form.
C
C  Description
C  ===========
C
C  ZLANTP returns the value
C
C     ZLANTP = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
C          Specifies the value to be returned in ZLANTP as described
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
C          The order of the matrix A.  N >= 0.  When N = 0, ZLANTP is
C          set to zero.
C
C  AP      (input) COMPLEX array, dimension (N*(N+1)/2)
C          The upper or lower triangular matrix A, packed columnwise in
C          a linear array.  The j-th column of A is stored in the array
C          AP as follows:
C          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
C          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
C          Note that when DIAG = 'U', the elements of the array AP
C          corresponding to the diagonal elements of the matrix A are
C          not referenced, but are assumed to be one.
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
      INTEGER                          N
      CHARACTER                        DIAG, NORM, UPLO
C     .. Array Arguments ..
      COMPLEX*16                       AP(*)
      DOUBLE PRECISION                 WORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION                 SCALE, SUM, VALUE
      INTEGER                          I, J, K
      LOGICAL                          UDIAG
C     .. External Subroutines ..
      EXTERNAL                         F06KJF
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
         K = 1
         IF ((DIAG.EQ.'U' .OR. DIAG.EQ.'u')) THEN
            VALUE = ONE
            IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
               DO 40 J = 1, N
                  DO 20 I = K, K + J - 2
                     VALUE = MAX(VALUE,ABS(AP(I)))
   20             CONTINUE
                  K = K + J
   40          CONTINUE
            ELSE
               DO 80 J = 1, N
                  DO 60 I = K + 1, K + N - J
                     VALUE = MAX(VALUE,ABS(AP(I)))
   60             CONTINUE
                  K = K + N - J + 1
   80          CONTINUE
            END IF
         ELSE
            VALUE = ZERO
            IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
               DO 120 J = 1, N
                  DO 100 I = K, K + J - 1
                     VALUE = MAX(VALUE,ABS(AP(I)))
  100             CONTINUE
                  K = K + J
  120          CONTINUE
            ELSE
               DO 160 J = 1, N
                  DO 140 I = K, K + N - J
                     VALUE = MAX(VALUE,ABS(AP(I)))
  140             CONTINUE
                  K = K + N - J + 1
  160          CONTINUE
            END IF
         END IF
      ELSE IF (((NORM.EQ.'O' .OR. NORM.EQ.'o')) .OR. (NORM.EQ.'1')) THEN
C
C        Find norm1(A).
C
         VALUE = ZERO
         K = 1
         UDIAG = (DIAG.EQ.'U' .OR. DIAG.EQ.'u')
         IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
            DO 220 J = 1, N
               IF (UDIAG) THEN
                  SUM = ONE
                  DO 180 I = K, K + J - 2
                     SUM = SUM + ABS(AP(I))
  180             CONTINUE
               ELSE
                  SUM = ZERO
                  DO 200 I = K, K + J - 1
                     SUM = SUM + ABS(AP(I))
  200             CONTINUE
               END IF
               K = K + J
               VALUE = MAX(VALUE,SUM)
  220       CONTINUE
         ELSE
            DO 280 J = 1, N
               IF (UDIAG) THEN
                  SUM = ONE
                  DO 240 I = K + 1, K + N - J
                     SUM = SUM + ABS(AP(I))
  240             CONTINUE
               ELSE
                  SUM = ZERO
                  DO 260 I = K, K + N - J
                     SUM = SUM + ABS(AP(I))
  260             CONTINUE
               END IF
               K = K + N - J + 1
               VALUE = MAX(VALUE,SUM)
  280       CONTINUE
         END IF
      ELSE IF ((NORM.EQ.'I' .OR. NORM.EQ.'i')) THEN
C
C        Find normI(A).
C
         K = 1
         IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
            IF ((DIAG.EQ.'U' .OR. DIAG.EQ.'u')) THEN
               DO 300 I = 1, N
                  WORK(I) = ONE
  300          CONTINUE
               DO 340 J = 1, N
                  DO 320 I = 1, J - 1
                     WORK(I) = WORK(I) + ABS(AP(K))
                     K = K + 1
  320             CONTINUE
                  K = K + 1
  340          CONTINUE
            ELSE
               DO 360 I = 1, N
                  WORK(I) = ZERO
  360          CONTINUE
               DO 400 J = 1, N
                  DO 380 I = 1, J
                     WORK(I) = WORK(I) + ABS(AP(K))
                     K = K + 1
  380             CONTINUE
  400          CONTINUE
            END IF
         ELSE
            IF ((DIAG.EQ.'U' .OR. DIAG.EQ.'u')) THEN
               DO 420 I = 1, N
                  WORK(I) = ONE
  420          CONTINUE
               DO 460 J = 1, N
                  K = K + 1
                  DO 440 I = J + 1, N
                     WORK(I) = WORK(I) + ABS(AP(K))
                     K = K + 1
  440             CONTINUE
  460          CONTINUE
            ELSE
               DO 480 I = 1, N
                  WORK(I) = ZERO
  480          CONTINUE
               DO 520 J = 1, N
                  DO 500 I = J, N
                     WORK(I) = WORK(I) + ABS(AP(K))
                     K = K + 1
  500             CONTINUE
  520          CONTINUE
            END IF
         END IF
         VALUE = ZERO
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
               K = 2
               DO 560 J = 2, N
                  CALL F06KJF(J-1,AP(K),1,SCALE,SUM)
                  K = K + J
  560          CONTINUE
            ELSE
               SCALE = ZERO
               SUM = ONE
               K = 1
               DO 580 J = 1, N
                  CALL F06KJF(J,AP(K),1,SCALE,SUM)
                  K = K + J
  580          CONTINUE
            END IF
         ELSE
            IF ((DIAG.EQ.'U' .OR. DIAG.EQ.'u')) THEN
               SCALE = ONE
               SUM = N
               K = 2
               DO 600 J = 1, N - 1
                  CALL F06KJF(N-J,AP(K),1,SCALE,SUM)
                  K = K + N - J + 1
  600          CONTINUE
            ELSE
               SCALE = ZERO
               SUM = ONE
               K = 1
               DO 620 J = 1, N
                  CALL F06KJF(N-J+1,AP(K),1,SCALE,SUM)
                  K = K + N - J + 1
  620          CONTINUE
            END IF
         END IF
         VALUE = SCALE*SQRT(SUM)
      END IF
C
      F06UKF = VALUE
      RETURN
C
C     End of F06UKF (ZLANTP)
C
      END
