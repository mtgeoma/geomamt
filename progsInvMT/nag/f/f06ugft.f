      DOUBLE PRECISION FUNCTION F06UGF(NORM,UPLO,N,AP,WORK)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     REAL                 ZLANSP
C     ENTRY                ZLANSP(NORM,UPLO,N,AP,WORK)
C
C  Purpose
C  =======
C
C  ZLANSP  returns the value of the one norm,  or the Frobenius norm, or
C  the  infinity norm,  or the  element of  largest absolute value  of a
C  complex symmetric matrix A,  supplied in packed form.
C
C  Description
C  ===========
C
C  ZLANSP returns the value
C
C     ZLANSP = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
C          Specifies the value to be returned in ZLANSP as described
C          above.
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the upper or lower triangular part of the
C          symmetric matrix A is supplied.
C          = 'U':  Upper triangular part of A is supplied
C          = 'L':  Lower triangular part of A is supplied
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.  When N = 0, ZLANSP is
C          set to zero.
C
C  AP      (input) COMPLEX array, dimension (N*(N+1)/2)
C          The upper or lower triangle of the symmetric matrix A, packed
C          columnwise in a linear array.  The j-th column of A is stored
C          in the array AP as follows:
C          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
C          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
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
      INTEGER                          N
      CHARACTER                        NORM, UPLO
C     .. Array Arguments ..
      COMPLEX*16                       AP(*)
      DOUBLE PRECISION                 WORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION                 ABSA, SCALE, SUM, VALUE
      INTEGER                          I, J, K
C     .. External Subroutines ..
      EXTERNAL                         F06KJF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, DIMAG, MAX, DBLE, SQRT
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
            K = 1
            DO 40 J = 1, N
               DO 20 I = K, K + J - 1
                  VALUE = MAX(VALUE,ABS(AP(I)))
   20          CONTINUE
               K = K + J
   40       CONTINUE
         ELSE
            K = 1
            DO 80 J = 1, N
               DO 60 I = K, K + N - J
                  VALUE = MAX(VALUE,ABS(AP(I)))
   60          CONTINUE
               K = K + N - J + 1
   80       CONTINUE
         END IF
      ELSE IF (((NORM.EQ.'I' .OR. NORM.EQ.'i'))
     *         .OR. ((NORM.EQ.'O' .OR. NORM.EQ.'o')) .OR. (NORM.EQ.'1'))
     *         THEN
C
C        Find normI(A) ( = norm1(A), since A is symmetric).
C
         VALUE = ZERO
         K = 1
         IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
            DO 120 J = 1, N
               SUM = ZERO
               DO 100 I = 1, J - 1
                  ABSA = ABS(AP(K))
                  SUM = SUM + ABSA
                  WORK(I) = WORK(I) + ABSA
                  K = K + 1
  100          CONTINUE
               WORK(J) = SUM + ABS(AP(K))
               K = K + 1
  120       CONTINUE
            DO 140 I = 1, N
               VALUE = MAX(VALUE,WORK(I))
  140       CONTINUE
         ELSE
            DO 160 I = 1, N
               WORK(I) = ZERO
  160       CONTINUE
            DO 200 J = 1, N
               SUM = WORK(J) + ABS(AP(K))
               K = K + 1
               DO 180 I = J + 1, N
                  ABSA = ABS(AP(K))
                  SUM = SUM + ABSA
                  WORK(I) = WORK(I) + ABSA
                  K = K + 1
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
         K = 2
         IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
            DO 220 J = 2, N
               CALL F06KJF(J-1,AP(K),1,SCALE,SUM)
               K = K + J
  220       CONTINUE
         ELSE
            DO 240 J = 1, N - 1
               CALL F06KJF(N-J,AP(K),1,SCALE,SUM)
               K = K + N - J + 1
  240       CONTINUE
         END IF
         SUM = 2*SUM
         K = 1
         DO 260 I = 1, N
            IF (DBLE(AP(K)).NE.ZERO) THEN
               ABSA = ABS(DBLE(AP(K)))
               IF (SCALE.LT.ABSA) THEN
                  SUM = ONE + SUM*(SCALE/ABSA)**2
                  SCALE = ABSA
               ELSE
                  SUM = SUM + (ABSA/SCALE)**2
               END IF
            END IF
            IF (DIMAG(AP(K)).NE.ZERO) THEN
               ABSA = ABS(DIMAG(AP(K)))
               IF (SCALE.LT.ABSA) THEN
                  SUM = ONE + SUM*(SCALE/ABSA)**2
                  SCALE = ABSA
               ELSE
                  SUM = SUM + (ABSA/SCALE)**2
               END IF
            END IF
            IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
               K = K + I + 1
            ELSE
               K = K + N - I + 1
            END IF
  260    CONTINUE
         VALUE = SCALE*SQRT(SUM)
      END IF
C
      F06UGF = VALUE
      RETURN
C
C     End of F06UGF (ZLANSP)
C
      END
