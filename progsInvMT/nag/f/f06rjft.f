      DOUBLE PRECISION FUNCTION F06RJF(NORM,UPLO,DIAG,M,N,A,LDA,WORK)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     REAL                 DLANTR
C     ENTRY                DLANTR(NORM,UPLO,DIAG,M,N,A,LDA,WORK)
C
C  Purpose
C  =======
C
C  DLANTR  returns the value of the one norm,  or the Frobenius norm, or
C  the  infinity norm,  or the  element of  largest absolute value  of a
C  trapezoidal or triangular matrix A.
C
C  Description
C  ===========
C
C  DLANTR returns the value
C
C     DLANTR = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
C          Specifies the value to be returned in DLANTR as described
C          above.
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the matrix A is upper or lower trapezoidal.
C          = 'U':  Upper trapezoidal
C          = 'L':  Lower trapezoidal
C          Note that A is triangular instead of trapezoidal if M = N.
C
C  DIAG    (input) CHARACTER*1
C          Specifies whether or not the matrix A has unit diagonal.
C          = 'N':  Non-unit diagonal
C          = 'U':  Unit diagonal
C
C  M       (input) INTEGER
C          The number of rows of the matrix A.  M >= 0, and if
C          UPLO = 'U', M <= N.  When M = 0, DLANTR is set to zero.
C
C  N       (input) INTEGER
C          The number of columns of the matrix A.  N >= 0, and if
C          UPLO = 'L', N <= M.  When N = 0, DLANTR is set to zero.
C
C  A       (input) REAL array, dimension (LDA,N)
C          The trapezoidal matrix A (A is triangular if M = N).
C          If UPLO = 'U', the leading m by n upper trapezoidal part of
C          the array A contains the upper trapezoidal matrix, and the
C          strictly lower triangular part of A is not referenced.
C          If UPLO = 'L', the leading m by n lower trapezoidal part of
C          the array A contains the lower trapezoidal matrix, and the
C          strictly upper triangular part of A is not referenced.  Note
C          that when DIAG = 'U', the diagonal elements of A are not
C          referenced and are assumed to be one.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(M,1).
C
C  WORK    (workspace) REAL array, dimension (LWORK),
C          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
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
      INTEGER                          LDA, M, N
      CHARACTER                        DIAG, NORM, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION                 A(LDA,*), WORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION                 SCALE, SUM, VALUE
      INTEGER                          I, J
      LOGICAL                          UDIAG
C     .. External Subroutines ..
      EXTERNAL                         F06FJF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, MAX, MIN, SQRT
C     .. Executable Statements ..
C
      IF (MIN(M,N).EQ.0) THEN
         VALUE = ZERO
      ELSE IF ((NORM.EQ.'M' .OR. NORM.EQ.'m')) THEN
C
C        Find max(abs(A(i,j))).
C
         IF ((DIAG.EQ.'U' .OR. DIAG.EQ.'u')) THEN
            VALUE = ONE
            IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
               DO 40 J = 1, N
                  DO 20 I = 1, MIN(M,J-1)
                     VALUE = MAX(VALUE,ABS(A(I,J)))
   20             CONTINUE
   40          CONTINUE
            ELSE
               DO 80 J = 1, N
                  DO 60 I = J + 1, M
                     VALUE = MAX(VALUE,ABS(A(I,J)))
   60             CONTINUE
   80          CONTINUE
            END IF
         ELSE
            VALUE = ZERO
            IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
               DO 120 J = 1, N
                  DO 100 I = 1, MIN(M,J)
                     VALUE = MAX(VALUE,ABS(A(I,J)))
  100             CONTINUE
  120          CONTINUE
            ELSE
               DO 160 J = 1, N
                  DO 140 I = J, M
                     VALUE = MAX(VALUE,ABS(A(I,J)))
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
               IF ((UDIAG) .AND. (J.LE.M)) THEN
                  SUM = ONE
                  DO 180 I = 1, J - 1
                     SUM = SUM + ABS(A(I,J))
  180             CONTINUE
               ELSE
                  SUM = ZERO
                  DO 200 I = 1, MIN(M,J)
                     SUM = SUM + ABS(A(I,J))
  200             CONTINUE
               END IF
               VALUE = MAX(VALUE,SUM)
  220       CONTINUE
         ELSE
            DO 280 J = 1, N
               IF (UDIAG) THEN
                  SUM = ONE
                  DO 240 I = J + 1, M
                     SUM = SUM + ABS(A(I,J))
  240             CONTINUE
               ELSE
                  SUM = ZERO
                  DO 260 I = J, M
                     SUM = SUM + ABS(A(I,J))
  260             CONTINUE
               END IF
               VALUE = MAX(VALUE,SUM)
  280       CONTINUE
         END IF
      ELSE IF ((NORM.EQ.'I' .OR. NORM.EQ.'i')) THEN
C
C        Find normI(A).
C
         IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
            IF ((DIAG.EQ.'U' .OR. DIAG.EQ.'u')) THEN
               DO 300 I = 1, M
                  WORK(I) = ONE
  300          CONTINUE
               DO 340 J = 1, N
                  DO 320 I = 1, MIN(M,J-1)
                     WORK(I) = WORK(I) + ABS(A(I,J))
  320             CONTINUE
  340          CONTINUE
            ELSE
               DO 360 I = 1, M
                  WORK(I) = ZERO
  360          CONTINUE
               DO 400 J = 1, N
                  DO 380 I = 1, MIN(M,J)
                     WORK(I) = WORK(I) + ABS(A(I,J))
  380             CONTINUE
  400          CONTINUE
            END IF
         ELSE
            IF ((DIAG.EQ.'U' .OR. DIAG.EQ.'u')) THEN
               DO 420 I = 1, N
                  WORK(I) = ONE
  420          CONTINUE
               DO 440 I = N + 1, M
                  WORK(I) = ZERO
  440          CONTINUE
               DO 480 J = 1, N
                  DO 460 I = J + 1, M
                     WORK(I) = WORK(I) + ABS(A(I,J))
  460             CONTINUE
  480          CONTINUE
            ELSE
               DO 500 I = 1, M
                  WORK(I) = ZERO
  500          CONTINUE
               DO 540 J = 1, N
                  DO 520 I = J, M
                     WORK(I) = WORK(I) + ABS(A(I,J))
  520             CONTINUE
  540          CONTINUE
            END IF
         END IF
         VALUE = ZERO
         DO 560 I = 1, M
            VALUE = MAX(VALUE,WORK(I))
  560    CONTINUE
      ELSE IF (((NORM.EQ.'F' .OR. NORM.EQ.'f'))
     *         .OR. ((NORM.EQ.'E' .OR. NORM.EQ.'e'))) THEN
C
C        Find normF(A).
C
         IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
            IF ((DIAG.EQ.'U' .OR. DIAG.EQ.'u')) THEN
               SCALE = ONE
               SUM = MIN(M,N)
               DO 580 J = 2, N
                  CALL F06FJF(MIN(M,J-1),A(1,J),1,SCALE,SUM)
  580          CONTINUE
            ELSE
               SCALE = ZERO
               SUM = ONE
               DO 600 J = 1, N
                  CALL F06FJF(MIN(M,J),A(1,J),1,SCALE,SUM)
  600          CONTINUE
            END IF
         ELSE
            IF ((DIAG.EQ.'U' .OR. DIAG.EQ.'u')) THEN
               SCALE = ONE
               SUM = MIN(M,N)
               DO 620 J = 1, N
                  CALL F06FJF(M-J,A(MIN(M,J+1),J),1,SCALE,SUM)
  620          CONTINUE
            ELSE
               SCALE = ZERO
               SUM = ONE
               DO 640 J = 1, N
                  CALL F06FJF(M-J+1,A(J,J),1,SCALE,SUM)
  640          CONTINUE
            END IF
         END IF
         VALUE = SCALE*SQRT(SUM)
      END IF
C
      F06RJF = VALUE
      RETURN
C
C     End of F06RJF (DLANTR)
C
      END
