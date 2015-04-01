      SUBROUTINE F07GWF(UPLO,N,AP,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             ZPPTRI(UPLO,N,AP,INFO)
C
C  Purpose
C  =======
C
C  ZPPTRI computes the inverse of a complex Hermitian positive definite
C  matrix A using the Cholesky factorization A = U'*U or A = L*L'
C  computed by F07GRF.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the factor stored in AP is upper or lower
C          triangular.
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  AP      (input) COMPLEX array, dimension (N*(N+1)/2)
C          On entry, the triangular factor U or L from the Cholesky
C          factorization A = U'*U or A = L*L', packed columnwise as a
C          linear array.  The j-th column of U or L is stored in the
C          array AP as follows:
C          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;
C          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n.
C
C          On exit, the upper or lower triangle of the (Hermitian)
C          inverse of A, overwriting the input factor U or L.
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C          > 0: if INFO = k, the (k,k) element of the factor U or L is
C               zero, and the inverse could not be computed.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        AP(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AJJ
      INTEGER           J, JC, JJ, JJN
      LOGICAL           UPPER
C     .. External Functions ..
      COMPLEX*16        ZDOTC
      EXTERNAL          ZDOTC
C     .. External Subroutines ..
      EXTERNAL          ZHPR, ZDSCAL, ZTPMV, F06AAZ, F07UWF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
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
         CALL F06AAZ('F07GWF/ZPPTRI',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
C     Invert the triangular Cholesky factor U or L.
C
      CALL F07UWF(UPLO,'Non-unit',N,AP,INFO)
      IF (INFO.GT.0) RETURN
      IF (UPPER) THEN
C
C        Compute the product inv(U) * inv(U)'.
C
         JJ = 0
         DO 20 J = 1, N
            JC = JJ + 1
            JJ = JJ + J
            IF (J.GT.1) CALL ZHPR('Upper',J-1,ONE,AP(JC),1,AP)
            AJJ = AP(JJ)
            CALL ZDSCAL(J,AJJ,AP(JC),1)
   20    CONTINUE
C
      ELSE
C
C        Compute the product inv(L)' * inv(L).
C
         JJ = 1
         DO 40 J = 1, N
            JJN = JJ + N - J + 1
            AP(JJ) = DBLE(ZDOTC(N-J+1,AP(JJ),1,AP(JJ),1))
            IF (J.LT.N) CALL ZTPMV('Lower','Conjugate transpose',
     *                             'Non-unit',N-J,AP(JJN),AP(JJ+1),1)
            JJ = JJN
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F07GWF (ZPPTRI)
C
      END
