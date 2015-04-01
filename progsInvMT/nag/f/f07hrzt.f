      SUBROUTINE F07HRZ(UPLO,N,KD,AB,LDAB,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             ZPBTF2(UPLO,N,KD,AB,LDAB,INFO)
C
C  Purpose
C  =======
C
C  ZPBTF2 computes the Cholesky factorization of a complex Hermitian
C  positive definite band matrix A.
C
C  The factorization has the form
C     A = U' * U ,  if UPLO = 'U', or
C     A = L  * L',  if UPLO = 'L',
C  where U is an upper triangular matrix, U' is the conjugate transpose
C  of U, and L is lower triangular.
C
C  This is the unblocked version of the algorithm, calling Level 2 BLAS.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the upper or lower triangular part of the
C          Hermitian matrix A is stored:
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  KD      (input) INTEGER
C          The number of super-diagonals of the matrix A if UPLO = 'U',
C          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.
C
C  AB      (input/output) COMPLEX array, dimension (LDAB,N)
C          On entry, the upper or lower triangle of the Hermitian band
C          matrix A, stored in the first KD+1 rows of the array.  The
C          j-th column of A is stored in the j-th column of the array AB
C          as follows:
C          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
C          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
C
C          On exit, if INFO = 0, the triangular factor U or L from the
C          Cholesky factorization A = U'*U or A = L*L' of the band
C          matrix A, in the same storage format as A.
C
C  LDAB    (input) INTEGER
C          The leading dimension of the array AB.  LDAB >= KD+1.
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C          > 0: if INFO = k, the leading minor of order k is not
C               positive definite, and the factorization could not be
C               completed.
C
C  Further Details
C  ===============
C
C  The band storage scheme is illustrated by the following example, when
C  N = 6, KD = 2, and UPLO = 'U':
C
C  On entry:                       On exit:
C
C      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
C      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
C     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
C
C  Similarly, if UPLO = 'L' the format of A is as follows:
C
C  On entry:                       On exit:
C
C     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
C     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
C     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
C
C  Array elements marked * are not used by the routine.
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
C     .. Scalar Arguments ..
      INTEGER           INFO, KD, LDAB, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        AB(LDAB,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AJJ
      INTEGER           J, KLD, KN
      LOGICAL           UPPER
C     .. External Functions ..
      COMPLEX*16        ZDOTC
      EXTERNAL          ZDOTC
C     .. External Subroutines ..
      EXTERNAL          ZHER, ZDSCAL, ZTRSV, F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN, DBLE, SQRT
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
      ELSE IF (KD.LT.0) THEN
         INFO = -3
      ELSE IF (LDAB.LT.KD+1) THEN
         INFO = -5
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07HRZ/ZPBTF2',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
      KLD = MAX(1,LDAB-1)
C
      IF (UPPER) THEN
C
C        Compute the Cholesky factorization A = U'*U.
C
         DO 20 J = 1, N
C
C           Compute elements J-KN:J-1 of column J.
C
            KN = MIN(J-1,KD)
            CALL ZTRSV('Upper','Conjugate transpose','Non-unit',KN,
     *                 AB(KD+1,J-KN),KLD,AB(KD+1-KN,J),1)
C
C           Compute U(J,J) and test for non-positive-definiteness.
C
            AJJ = DBLE(AB(KD+1,J)) - ZDOTC(KN,AB(KD+1-KN,J),1,
     *            AB(KD+1-KN,J),1)
            IF (AJJ.LE.ZERO) THEN
               AB(KD+1,J) = AJJ
               GO TO 60
            END IF
            AB(KD+1,J) = SQRT(AJJ)
   20    CONTINUE
      ELSE
C
C        Compute the Cholesky factorization A = L*L'.
C
         DO 40 J = 1, N
C
C           Compute L(J,J) and test for non-positive-definiteness.
C
            AJJ = DBLE(AB(1,J))
            IF (AJJ.LE.ZERO) THEN
               AB(1,J) = AJJ
               GO TO 60
            END IF
            AJJ = SQRT(AJJ)
            AB(1,J) = AJJ
C
C           Compute elements J+1:J+KN of column J and update the
C           trailing submatrix within the band.
C
            KN = MIN(KD,N-J)
            IF (KN.GT.0) THEN
               CALL ZDSCAL(KN,ONE/AJJ,AB(2,J),1)
               CALL ZHER('Lower',KN,-ONE,AB(2,J),1,AB(1,J+1),KLD)
            END IF
   40    CONTINUE
      END IF
      RETURN
C
   60 CONTINUE
      INFO = J
      RETURN
C
C     End of F07HRZ (ZPBTF2)
C
      END
