      SUBROUTINE F07GRF(UPLO,N,AP,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             ZPPTRF(UPLO,N,AP,INFO)
C
C  Purpose
C  =======
C
C  ZPPTRF computes the Cholesky factorization of a complex Hermitian
C  positive definite matrix stored in packed format.
C
C  The factorization has the form
C     A = U' * U ,  if UPLO = 'U', or
C     A = L  * L',  if UPLO = 'L',
C  where U is an upper triangular matrix and L is lower triangular.
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
C  AP      (input/output) COMPLEX array, dimension (N*(N+1)/2)
C          On entry, the upper or lower triangle of the Hermitian matrix
C          A, packed columnwise in a linear array.  The j-th column of A
C          is stored in the array AP as follows:
C          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
C          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
C          See below for further details.
C
C          On exit, if INFO = 0, the triangular factor U or L from the
C          Cholesky factorization A = U'*U or A = L*L', in the same
C          storage format as A.
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C          > 0: if INFO = k, the leading minor of order k is not
C               positive definite, and the factorization could not be
C               completed.
C
C  Further Details
C  ======= =======
C
C  The packed storage scheme is illustrated by the following example
C  when N = 4, UPLO = 'U':
C
C  Two-dimensional storage of the Hermitian matrix A:
C
C     a11 a12 a13 a14
C         a22 a23 a24
C             a33 a34     (aij = aji)
C                 a44
C
C  Packed storage of the upper triangle of A:
C
C  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        AP(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AJJ
      INTEGER           J, JC, JJ
      LOGICAL           UPPER
C     .. External Functions ..
      COMPLEX*16        ZDOTC
      EXTERNAL          ZDOTC
C     .. External Subroutines ..
      EXTERNAL          ZHPR, ZDSCAL, ZTPSV, F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
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
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07GRF/ZPPTRF',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
      IF (UPPER) THEN
C
C        Compute the Cholesky factorization A = U'*U.
C
         JJ = 0
         DO 20 J = 1, N
            JC = JJ + 1
            JJ = JJ + J
C
C           Compute elements 1:J-1 of column J.
C
            IF (J.GT.1) CALL ZTPSV('Upper','Conjugate transpose',
     *                             'Non-unit',J-1,AP,AP(JC),1)
C
C           Compute U(J,J) and test for non-positive-definiteness.
C
            AJJ = DBLE(AP(JJ)) - ZDOTC(J-1,AP(JC),1,AP(JC),1)
            IF (AJJ.LE.ZERO) THEN
               AP(JJ) = AJJ
               GO TO 60
            END IF
            AP(JJ) = SQRT(AJJ)
   20    CONTINUE
      ELSE
C
C        Compute the Cholesky factorization A = L*L'.
C
         JJ = 1
         DO 40 J = 1, N
C
C           Compute L(J,J) and test for non-positive-definiteness.
C
            AJJ = DBLE(AP(JJ))
            IF (AJJ.LE.ZERO) THEN
               AP(JJ) = AJJ
               GO TO 60
            END IF
            AJJ = SQRT(AJJ)
            AP(JJ) = AJJ
C
C           Compute elements J+1:N of column J and update the trailing
C           submatrix.
C
            IF (J.LT.N) THEN
               CALL ZDSCAL(N-J,ONE/AJJ,AP(JJ+1),1)
               CALL ZHPR('Lower',N-J,-ONE,AP(JJ+1),1,AP(JJ+N-J+1))
               JJ = JJ + N - J + 1
            END IF
   40    CONTINUE
      END IF
      GO TO 80
C
   60 CONTINUE
      INFO = J
C
   80 CONTINUE
      RETURN
C
C     End of F07GRF (ZPPTRF)
C
      END
